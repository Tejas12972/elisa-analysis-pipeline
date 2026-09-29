"""Plate-reader export and plate-layout parsing.

A sandwich ELISA run produces two independent pieces of information:

1. **The plate export**: one optical density (OD) per well, written by the
   plate reader. Instruments disagree on the format. Most write an 8 x 12
   grid that mirrors the physical plate (rows A-H, columns 1-12), often with
   instrument metadata above it. Others write a "long" table with one
   ``(well, OD)`` pair per line.
2. **The plate layout**: what the scientist *put* in each well. That is the
   well type (standard, sample, blank, control), the sample identity, the
   nominal concentration for standards, the dilution factor, and the
   treatment group.

This module turns both into one tidy DataFrame (one row per well) and
subtracts the mean blank OD. Every later stage builds on that frame.

Validation philosophy
---------------------
Errors in these files cause silent scientific mistakes, not crashes. For
example, a mistyped standard concentration shifts every back-calculated
result on the plate. So the parsers are strict. They collect *every*
problem they find and raise it in one exception, which saves the user from
fixing errors one at a time. Only conditions that are unusual but
legitimate trigger a :class:`PlateWarning`, such as slightly negative raw
ODs from reader-side blanking.
"""

from __future__ import annotations

import csv
import io as _stdio
import re
import warnings
from collections.abc import Iterable
from pathlib import Path

import numpy as np
import pandas as pd

ROWS: tuple[str, ...] = tuple("ABCDEFGH")
COLUMNS: tuple[int, ...] = tuple(range(1, 13))
ALL_WELLS: tuple[str, ...] = tuple(f"{r}{c}" for r in ROWS for c in COLUMNS)

WELL_TYPES: tuple[str, ...] = ("standard", "sample", "blank", "control")

# Readers print a token instead of a number when a well saturates the detector.
# Such a well has a real signal that is simply too high to report. It is
# neither missing nor garbage, so it is kept as NaN with status "overflow" and
# the QC layer can report it.
OVERFLOW_TOKENS: frozenset[str] = frozenset(
    {"OVRFLW", "OVERFLOW", "OVER", "SAT", "SATURATED", "HIGH", "HI"}
)

_EXCEL_SUFFIXES = {".xlsx", ".xlsm", ".xls"}
_WELL_RE = re.compile(r"^\s*([A-Ha-h])\s*0*(1[0-2]|[1-9])\s*$")

_WELL_COL_ALIASES = ("well", "well_id", "wellid", "well_position", "position", "pos")
_OD_COL_RE = re.compile(r"^(od|abs|absorbance|value|signal|raw|reading)", re.IGNORECASE)

_LAYOUT_ALIASES: dict[str, tuple[str, ...]] = {
    "well": _WELL_COL_ALIASES,
    "type": ("type", "well_type", "role", "content"),
    "sample_id": ("sample_id", "sampleid", "sample", "id", "name"),
    "concentration": (
        "concentration",
        "conc",
        "nominal_concentration",
        "nominal_conc",
        "std_conc",
    ),
    "dilution": ("dilution", "dilution_factor", "df", "dil"),
    "group": ("group", "treatment", "treatment_group", "condition"),
}
_TYPE_ALIASES: dict[str, str] = {
    "standard": "standard",
    "std": "standard",
    "s": "standard",
    "calibrator": "standard",
    "sample": "sample",
    "unknown": "sample",
    "unk": "sample",
    "u": "sample",
    "blank": "blank",
    "blk": "blank",
    "b": "blank",
    "zero": "blank",
    "control": "control",
    "ctrl": "control",
    "qc": "control",
    "c": "control",
}


class InputFileError(ValueError):
    """Raised when a plate or layout file cannot be used as-is.

    Attributes:
        problems: Every individual problem found, one message each. This way
            a UI can list them all instead of stopping at the first.
    """

    def __init__(self, source: str | Path, problems: Iterable[str]):
        self.source = str(source)
        self.problems = list(problems)
        bullet_list = "\n".join(f"  - {p}" for p in self.problems)
        super().__init__(f"{self.source}: {len(self.problems)} problem(s)\n{bullet_list}")


class PlateFormatError(InputFileError):
    """The plate export is unreadable or contains invalid OD values."""


class LayoutError(InputFileError):
    """The plate layout is missing information or is internally inconsistent."""


class MergeError(InputFileError):
    """The plate and layout disagree, e.g. the layout uses a well with no reading."""


class PlateWarning(UserWarning):
    """A condition that is legitimate but worth a human look."""


# --------------------------------------------------------------------------- #
# Well-name helpers
# --------------------------------------------------------------------------- #
def normalize_well(name: object) -> str | None:
    """Return the canonical well name (``"A1"``), or None if *name* is not a well.

    Accepts the common spellings: ``"A1"``, ``"a01"``, and ``" A 1 "``.
    """
    if not isinstance(name, str):
        return None
    m = _WELL_RE.match(name)
    return f"{m.group(1).upper()}{int(m.group(2))}" if m else None


def well_row(well: str) -> str:
    """Row letter of a canonical well name (``"B7"`` -> ``"B"``)."""
    return well[0]


def well_column(well: str) -> int:
    """Column number of a canonical well name (``"B7"`` -> ``7``)."""
    return int(well[1:])


# --------------------------------------------------------------------------- #
# Low-level file reading
# --------------------------------------------------------------------------- #
def _is_excel(path: Path) -> bool:
    return path.suffix.lower() in _EXCEL_SUFFIXES


def _read_text(path: Path) -> str:
    # Instrument software on Windows frequently writes cp1252 (e.g. "µg/mL").
    for encoding in ("utf-8-sig", "cp1252"):
        try:
            return path.read_text(encoding=encoding)
        except UnicodeDecodeError:
            continue
    raise PlateFormatError(path, ["file is not valid UTF-8 or Windows-1252 text"])


def _read_raw_cells(path: Path, sheet: str | int | None) -> pd.DataFrame:
    """Read a file as a header-less grid of raw cells (object dtype).

    Plate exports are *ragged*: a few metadata lines, then a table, sometimes
    several tables. So CSVs go through :mod:`csv` rather than
    :func:`pandas.read_csv`, which insists on a rectangular file. The
    delimiter (comma, tab, or semicolon) is sniffed.
    """
    if not path.exists():
        raise FileNotFoundError(path)
    if _is_excel(path):
        return pd.read_excel(
            path, sheet_name=sheet if sheet is not None else 0, header=None, dtype=object
        )

    text = _read_text(path)
    rows = list(csv.reader(_stdio.StringIO(text), delimiter=_guess_delimiter(text)))
    width = max((len(r) for r in rows), default=0)
    padded = [r + [""] * (width - len(r)) for r in rows]
    return pd.DataFrame(padded, dtype=object)


def _guess_delimiter(text: str) -> str:
    """Pick the delimiter under which the most lines share one field count (>= 2).

    :class:`csv.Sniffer` is unreliable on ragged exports, because a metadata
    line with no delimiter defeats it. It also picks ',' for semicolon files
    that use decimal commas. A table is the set of lines that split into the
    same number of fields, so the best delimiter is the one that produces
    the largest such set.
    """
    lines = [ln for ln in text.splitlines()[:200] if ln.strip()]
    best, best_score = ",", 0
    for delim in (",", "\t", ";"):
        counts: dict[int, int] = {}
        for fields in csv.reader(lines, delimiter=delim):
            if len(fields) >= 2:
                counts[len(fields)] = counts.get(len(fields), 0) + 1
        score = max(counts.values(), default=0)
        if score > best_score:
            best, best_score = delim, score
    return best


def _cell_str(value: object) -> str:
    """Stringify a raw cell. NaN becomes '' and Excel's 1.0 becomes '1'."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ""
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value).strip()


def _parse_od(value: object) -> tuple[float, str] | None:
    """Parse one OD cell into ``(od, status)``, or None if unparseable.

    Status is ``"ok"``, ``"missing"`` (empty cell, i.e. well not read), or
    ``"overflow"`` (detector saturated).
    """
    text = _cell_str(value)
    if text == "":
        return np.nan, "missing"
    if text.upper() in OVERFLOW_TOKENS or text.startswith(">"):
        return np.nan, "overflow"
    # Semicolon-delimited exports from European-locale readers write "0,734".
    # An OD never needs a thousands separator, so a lone comma is a decimal point.
    if text.count(",") == 1 and "." not in text:
        text = text.replace(",", ".")
    try:
        od = float(text)
    except ValueError:
        return None
    if not np.isfinite(od):
        return None
    return od, "ok"


# --------------------------------------------------------------------------- #
# Plate exports
# --------------------------------------------------------------------------- #
def _find_grid_blocks(cells: pd.DataFrame) -> list[tuple[int, int]]:
    """Locate 8 x 12 blocks. Returns ``(header_row, first_value_col)`` pairs.

    A block is a row whose cells read 1, 2, ..., 12 left to right, followed
    by eight rows whose cell just left of the "1" column reads A ... H.
    """
    strs = cells.map(_cell_str)
    n_rows, n_cols = strs.shape
    blocks: list[tuple[int, int]] = []
    expected_header = [str(c) for c in COLUMNS]
    for i in range(n_rows - len(ROWS)):
        row = strs.iloc[i].tolist()
        for j in range(1, n_cols - len(COLUMNS) + 1):
            if row[j : j + len(COLUMNS)] != expected_header:
                continue
            labels = [strs.iat[i + 1 + k, j - 1].upper() for k in range(len(ROWS))]
            if labels == list(ROWS):
                blocks.append((i, j))
    return blocks


def _grid_to_long(cells: pd.DataFrame, top: int, left: int) -> list[tuple[str, object]]:
    return [
        (f"{r}{c}", cells.iat[top + ri, left + ci])
        for ri, r in enumerate(ROWS)
        for ci, c in enumerate(COLUMNS)
    ]


def _find_long_header(
    cells: pd.DataFrame, well_col: str | None, value_col: str | None
) -> tuple[int, int, int] | None:
    """Find ``(header_row, well_col_idx, value_col_idx)`` for a long-format table."""
    strs = cells.map(_cell_str)
    for i in range(len(strs)):
        headers = [h.lower().replace(" ", "_") for h in strs.iloc[i].tolist()]
        if well_col is not None:
            w_idx = headers.index(well_col.lower()) if well_col.lower() in headers else None
        else:
            w_idx = next((k for k, h in enumerate(headers) if h in _WELL_COL_ALIASES), None)
        if w_idx is None:
            continue
        if value_col is not None:
            v_idx = headers.index(value_col.lower()) if value_col.lower() in headers else None
        else:
            v_idx = next(
                (k for k, h in enumerate(headers) if k != w_idx and _OD_COL_RE.match(h)), None
            )
        if v_idx is not None:
            return i, w_idx, v_idx
    return None


def read_plate(
    path: str | Path,
    *,
    plate_id: str | None = None,
    sheet: str | int | None = None,
    block: int | None = None,
    well_col: str | None = None,
    value_col: str | None = None,
) -> pd.DataFrame:
    """Read a plate-reader export in grid or long format (CSV, TSV, or Excel).

    The format is detected automatically:

    * **Grid**: an 8 x 12 block with a column-number header ``1 ... 12`` and
      row labels ``A ... H``. Metadata lines above or around the block are
      ignored. Some readers write several blocks, e.g. one per wavelength.
      In that case pass ``block`` (0-based) to choose one.
    * **Long**: a header row with a well column (``Well``, ``Position``, ...)
      and an OD column (``OD``, ``OD450``, ``Absorbance``, ``Value``, ...).
      Pass ``well_col`` / ``value_col`` to override the detection.

    Args:
        path: File to read.
        plate_id: Identifier stored in the ``plate_id`` column. Defaults to
            the file stem.
        sheet: Excel sheet name or index (default: first sheet).
        block: Which grid block to use when the file contains several.
        well_col: Explicit name of the well column (long format).
        value_col: Explicit name of the OD column (long format).

    Returns:
        A DataFrame with exactly 96 rows in plate order and the columns
        ``plate_id, well, row, column, od, od_status``. ``od_status`` is
        ``"ok"``, ``"missing"`` (no reading for that well), or ``"overflow"``
        (reader reported saturation; ``od`` is NaN).

    Raises:
        PlateFormatError: No recognizable table, ambiguous blocks, invalid
            well names, duplicate wells, or unparseable OD values.
    """
    path = Path(path)
    cells = _read_raw_cells(path, sheet)
    pid = plate_id if plate_id is not None else path.stem

    blocks = _find_grid_blocks(cells) if (well_col is None and value_col is None) else []
    pairs: list[tuple[str, object]]
    if blocks:
        if len(blocks) > 1 and block is None:
            raise PlateFormatError(
                path,
                [
                    f"found {len(blocks)} 8x12 grids (header rows "
                    f"{', '.join(str(b[0] + 1) for b in blocks)}); pass block=<0-based index> "
                    "to choose one, e.g. the primary wavelength"
                ],
            )
        idx = block if block is not None else 0
        if not 0 <= idx < len(blocks):
            raise PlateFormatError(path, [f"block={idx} requested but only {len(blocks)} found"])
        header_row, left = blocks[idx]
        pairs = _grid_to_long(cells, header_row + 1, left)
    elif cells.shape == (len(ROWS), len(COLUMNS)):
        # Bare 8x12 block with no labels: assume it is the plate as printed.
        pairs = _grid_to_long(cells, 0, 0)
    else:
        found = _find_long_header(cells, well_col, value_col)
        if found is None:
            raise PlateFormatError(
                path,
                [
                    "no plate data found: expected an 8x12 grid with a 1..12 header and "
                    "A..H row labels, or a table with a well column and an OD column"
                ],
            )
        header_row, w_idx, v_idx = found
        body = cells.iloc[header_row + 1 :]
        pairs = [
            (_cell_str(w), v)
            for w, v in zip(body.iloc[:, w_idx], body.iloc[:, v_idx], strict=True)
            if _cell_str(w) != ""
        ]

    return _pairs_to_plate(path, pid, pairs)


def _pairs_to_plate(path: Path, plate_id: str, pairs: list[tuple[str, object]]) -> pd.DataFrame:
    problems: list[str] = []
    records: dict[str, tuple[float, str]] = {}
    for raw_well, raw_value in pairs:
        well = normalize_well(raw_well)
        if well is None:
            problems.append(f"invalid well name {raw_well!r} (expected A1-H12)")
            continue
        if well in records:
            problems.append(f"well {well} appears more than once")
            continue
        parsed = _parse_od(raw_value)
        if parsed is None:
            problems.append(f"well {well}: cannot parse OD value {_cell_str(raw_value)!r}")
            continue
        records[well] = parsed
    if problems:
        raise PlateFormatError(path, problems)

    od = [records.get(w, (np.nan, "missing"))[0] for w in ALL_WELLS]
    status = [records.get(w, (np.nan, "missing"))[1] for w in ALL_WELLS]
    plate = pd.DataFrame(
        {
            "plate_id": plate_id,
            "well": list(ALL_WELLS),
            "row": [well_row(w) for w in ALL_WELLS],
            "column": [well_column(w) for w in ALL_WELLS],
            "od": np.asarray(od, dtype=float),
            "od_status": status,
        }
    )

    if (plate["od_status"] != "missing").sum() == 0:
        raise PlateFormatError(path, ["plate contains no OD readings"])
    negative = plate.loc[plate["od"] < 0, "well"].tolist()
    if negative:
        warnings.warn(
            f"{path.name}: {len(negative)} well(s) have negative raw OD "
            f"({', '.join(negative[:8])}{'...' if len(negative) > 8 else ''}). This is "
            "usually reader-side blanking. Values are kept as-is.",
            PlateWarning,
            stacklevel=2,
        )
    return plate


# --------------------------------------------------------------------------- #
# Layouts
# --------------------------------------------------------------------------- #
def _canonical_layout_columns(columns: Iterable[str]) -> dict[str, str]:
    """Map the user's column names to canonical names; unknown columns pass through."""
    mapping: dict[str, str] = {}
    taken: set[str] = set()
    for original in columns:
        key = str(original).strip().lower().replace(" ", "_").replace("-", "_")
        canonical = next(
            (c for c, aliases in _LAYOUT_ALIASES.items() if key in aliases and c not in taken), key
        )
        taken.add(canonical)
        mapping[original] = canonical
    return mapping


def read_layout(path: str | Path, *, sheet: str | int | None = None) -> pd.DataFrame:
    """Read and validate a plate layout (long format, CSV/TSV/Excel).

    Required columns: ``well`` and ``type``. Optional columns: ``sample_id``,
    ``concentration`` (standards only), ``dilution`` (default 1), and
    ``group``. Common synonyms are accepted, e.g. ``Well Type``, ``Conc``,
    ``Dilution Factor``, ``Treatment``. Any extra columns (``analyte``,
    ``notes``, ...) are kept unchanged. Wells not listed are treated as empty.

    Validation rules and their reasons:

    * ``type`` must be standard, sample, blank, or control (abbreviations
      like ``std``, ``unk``, ``blk``, and ``qc`` are accepted).
    * Standards need a positive concentration. A zero standard must be
      typed ``blank``. Zero cannot sit on a log-concentration axis, and it
      is what blank subtraction uses anyway.
    * Only standards may carry a concentration. A number next to a sample
      is either a typo or an expected value that belongs in the control
      specification (Phase 4). Guessing which would be unsafe.
    * Samples and controls need a ``sample_id``, which defines their
      replicate groups.
    * ``dilution`` is a fold-dilution >= 1 (``4`` means 1:4). Values below
      1 usually mean the user entered a fraction (0.25), which would scale
      concentrations the wrong way.
    * One ``sample_id`` must map to one type and one group, and a standard
      ID to one concentration. A conflict usually means a copy-paste error.

    Returns:
        A DataFrame with the columns ``well, row, column, type, sample_id,
        concentration, dilution, group`` plus any extra columns, sorted in
        plate order.

    Raises:
        LayoutError: Every problem found, listed in one exception.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    if _is_excel(path):
        raw = pd.read_excel(path, sheet_name=sheet if sheet is not None else 0, dtype=object)
    else:
        raw = pd.read_csv(
            _stdio.StringIO(_read_text(path)),
            sep=None,
            engine="python",
            dtype=object,
            skip_blank_lines=True,
        )
    raw = raw.dropna(how="all")
    layout = raw.rename(columns=_canonical_layout_columns(raw.columns))

    missing_cols = [c for c in ("well", "type") if c not in layout.columns]
    if missing_cols:
        raise LayoutError(
            path,
            [
                f"missing required column(s): {', '.join(missing_cols)} "
                f"(found: {', '.join(map(str, raw.columns))})"
            ],
        )
    for col in ("sample_id", "concentration", "dilution", "group"):
        if col not in layout.columns:
            layout[col] = np.nan

    problems: list[str] = []
    clean_rows: list[dict[str, object]] = []
    seen_wells: set[str] = set()
    extra_cols = [
        c
        for c in layout.columns
        if c not in ("well", "type", "sample_id", "concentration", "dilution", "group")
    ]

    for line_no, rec in zip(range(2, len(layout) + 2), layout.to_dict("records"), strict=True):
        where = f"line {line_no}"
        well = normalize_well(_cell_str(rec["well"]))
        if well is None:
            problems.append(f"{where}: invalid well name {_cell_str(rec['well'])!r}")
            continue
        where = f"well {well}"
        if well in seen_wells:
            problems.append(f"{where}: listed more than once")
            continue
        seen_wells.add(well)

        type_text = _cell_str(rec["type"]).lower()
        well_type = _TYPE_ALIASES.get(type_text)
        if well_type is None:
            problems.append(
                f"{where}: unknown type {type_text!r} (expected one of {', '.join(WELL_TYPES)})"
            )
            continue

        sample_id = _cell_str(rec["sample_id"]) or None
        group = _cell_str(rec["group"]) or None

        conc_text = _cell_str(rec["concentration"])
        conc = np.nan
        if conc_text:
            try:
                conc = float(conc_text)
            except ValueError:
                problems.append(f"{where}: concentration {conc_text!r} is not a number")
        if well_type == "standard":
            if not conc_text:
                problems.append(f"{where}: standard is missing its nominal concentration")
            elif np.isfinite(conc) and conc <= 0:
                problems.append(
                    f"{where}: standard concentration must be > 0 "
                    f"(got {conc_text}); type zero standards as 'blank'"
                )
            if sample_id is None and np.isfinite(conc):
                sample_id = f"STD_{conc:g}"
        elif conc_text:
            problems.append(
                f"{where}: only standards may have a concentration "
                f"(type is {well_type!r}); put control targets in the QC spec"
            )

        dil_text = _cell_str(rec["dilution"])
        dilution = 1.0
        if dil_text:
            try:
                dilution = float(dil_text)
            except ValueError:
                problems.append(f"{where}: dilution {dil_text!r} is not a number")
            else:
                if not np.isfinite(dilution) or dilution < 1:
                    problems.append(
                        f"{where}: dilution must be a fold-dilution >= 1 "
                        f"(got {dil_text}; write 4 for a 1:4 dilution)"
                    )

        if well_type in ("sample", "control") and sample_id is None:
            problems.append(f"{where}: {well_type} wells need a sample_id")
        if well_type == "blank" and sample_id is None:
            sample_id = "BLANK"

        clean_rows.append(
            {
                "well": well,
                "row": well_row(well),
                "column": well_column(well),
                "type": well_type,
                "sample_id": sample_id,
                "concentration": conc,
                "dilution": dilution,
                "group": group,
                **{c: rec[c] for c in extra_cols},
            }
        )

    if not clean_rows and not problems:
        problems.append("layout contains no wells")
    if problems:
        raise LayoutError(path, problems)

    layout_df = pd.DataFrame(clean_rows)
    problems.extend(_check_id_consistency(layout_df))
    if problems:
        raise LayoutError(path, problems)

    order = {w: i for i, w in enumerate(ALL_WELLS)}
    return layout_df.sort_values("well", key=lambda s: s.map(order)).reset_index(drop=True)


def _check_id_consistency(layout: pd.DataFrame) -> list[str]:
    """Each sample_id maps to one type and one group, and each standard ID to one level."""
    problems: list[str] = []
    for sid, grp in layout.groupby("sample_id", dropna=True, sort=True):
        wells = ", ".join(grp["well"])
        types = sorted(grp["type"].unique())
        if len(types) > 1:
            problems.append(
                f"sample_id {sid!r} is used for several types ({', '.join(types)}) in wells {wells}"
            )
            continue
        groups = sorted(grp["group"].dropna().unique())
        if len(groups) > 1:
            problems.append(
                f"sample_id {sid!r} is assigned to several groups "
                f"({', '.join(groups)}) in wells {wells}"
            )
        if types[0] == "standard":
            levels = grp["concentration"].dropna().unique()
            if len(levels) > 1:
                problems.append(
                    f"standard {sid!r} has conflicting concentrations "
                    f"({', '.join(f'{v:g}' for v in sorted(levels))}) "
                    f"in wells {wells}"
                )
    return problems


# --------------------------------------------------------------------------- #
# Merge + blank subtraction
# --------------------------------------------------------------------------- #
def merge_plate_layout(
    plate: pd.DataFrame, layout: pd.DataFrame, *, source: str = "plate/layout"
) -> pd.DataFrame:
    """Join OD readings onto the layout, keeping only wells the layout uses.

    * A layout well with **no reading** is an error. It usually means the
      layout was written for a different plate map or the export is
      truncated.
    * A well with a reading but **no layout entry** is dropped with a
      warning. Unused wells are common, but a whole missing column in the
      layout is worth noticing.
    * **Overflow** wells are kept (``od`` NaN, ``od_status`` "overflow").
      For a standard or sample, saturation carries information ("above
      range") that QC should report rather than hide.

    Returns:
        The layout columns plus ``plate_id, od, od_status``, in plate order.
    """
    plate_cols = ["plate_id", "well", "od", "od_status"]
    merged = layout.merge(plate[plate_cols], on="well", how="left", validate="one_to_one")

    no_reading = merged.loc[
        merged["od_status"].isna() | (merged["od_status"] == "missing"), "well"
    ].tolist()
    if no_reading:
        raise MergeError(
            source,
            [
                f"layout assigns {len(no_reading)} well(s) that have no OD reading: "
                f"{', '.join(no_reading)}"
            ],
        )

    unassigned = plate.loc[
        (plate["od_status"] != "missing") & ~plate["well"].isin(layout["well"]), "well"
    ].tolist()
    if unassigned:
        warnings.warn(
            f"{source}: {len(unassigned)} well(s) have readings but no layout entry "
            f"and were ignored: {', '.join(unassigned[:12])}"
            f"{'...' if len(unassigned) > 12 else ''}",
            PlateWarning,
            stacklevel=2,
        )

    front = [
        "plate_id",
        "well",
        "row",
        "column",
        "type",
        "sample_id",
        "concentration",
        "dilution",
        "group",
        "od",
        "od_status",
    ]
    rest = [c for c in merged.columns if c not in front]
    return merged[front + rest]


def subtract_blank(df: pd.DataFrame) -> pd.DataFrame:
    """Subtract each plate's mean blank OD. Adds ``od_net`` and blank statistics.

    **Why:** blank wells contain every reagent except the analyte. Their OD is
    the background from substrate auto-oxidation, non-specific binding of
    the detection antibody, and the plastic itself. Subtracting it makes a
    standard curve's lower asymptote represent "no analyte" rather than
    "no analyte + background". It also makes curves comparable across plates.

    **Why the mean:** the arithmetic mean is the conventional estimate and
    matches the LOD definition used in QC (mean blank + 3 SD). Blank outliers
    are detected in the QC stage, not hidden here by a robust estimator.

    **Why negative net ODs are kept:** a sample reading slightly below the
    blank is measurement noise around zero. Clipping to zero would bias
    low-end replicate means upward and shrink their SD. Back-calculation
    (Phase 3) flags these wells as below the curve's lower asymptote.

    Adds columns ``od_net``, ``blank_mean``, ``blank_sd`` (sample SD, NaN
    with a single blank), and ``n_blanks``. The blank statistics repeat on
    every row of a plate so they survive concatenation of plates.

    Raises:
        MergeError: A plate has no readable blank wells.
    """
    if "plate_id" not in df.columns:
        raise KeyError("subtract_blank expects a 'plate_id' column; use merge_plate_layout")
    out = df.copy()
    problems: list[str] = []
    stats: dict[object, tuple[float, float, int]] = {}
    for pid, grp in out.groupby("plate_id", sort=False):
        blanks = grp.loc[(grp["type"] == "blank") & (grp["od_status"] == "ok"), "od"]
        if blanks.empty:
            problems.append(
                f"plate {pid!r} has no readable blank wells; blank subtraction needs at least one"
            )
            continue
        sd = float(blanks.std(ddof=1)) if len(blanks) > 1 else np.nan
        stats[pid] = (float(blanks.mean()), sd, len(blanks))
        if len(blanks) < 2:
            warnings.warn(
                f"plate {pid!r} has a single blank well; blank SD and LOD cannot be estimated",
                PlateWarning,
                stacklevel=2,
            )
    if problems:
        raise MergeError("blank subtraction", problems)

    out["blank_mean"] = out["plate_id"].map(lambda p: stats[p][0])
    out["blank_sd"] = out["plate_id"].map(lambda p: stats[p][1])
    out["n_blanks"] = out["plate_id"].map(lambda p: stats[p][2])
    out["od_net"] = out["od"] - out["blank_mean"]
    return out


def load_plate(
    plate_path: str | Path,
    layout_path: str | Path,
    *,
    plate_id: str | None = None,
    plate_kwargs: dict[str, object] | None = None,
    layout_kwargs: dict[str, object] | None = None,
) -> pd.DataFrame:
    """Read a plate export and its layout, merge them, and blank-subtract.

    This is the single call most users need. The extra keyword dictionaries
    are forwarded to :func:`read_plate` and :func:`read_layout`.

    Example:
        >>> df = load_plate("data/simulated/example_plate_grid.csv",
        ...                 "data/layouts/example_layout.csv")  # doctest: +SKIP
    """
    plate = read_plate(plate_path, plate_id=plate_id, **(plate_kwargs or {}))  # type: ignore[arg-type]
    layout = read_layout(layout_path, **(layout_kwargs or {}))  # type: ignore[arg-type]
    merged = merge_plate_layout(
        plate, layout, source=f"{Path(plate_path).name} + {Path(layout_path).name}"
    )
    return subtract_blank(merged)


# --------------------------------------------------------------------------- #
# SoftMax Pro "Plate + Groups" text exports
# --------------------------------------------------------------------------- #
class SoftMaxPlate:
    """One plate from a SoftMax Pro text export.

    Attributes:
        name: The plate's name in SoftMax (not unique; often "Plate1" repeated).
        index: 0-based position of the plate in the file.
        plate: Tidy 96-row plate frame, as returned by :func:`read_plate`.
        raw: Per-wavelength raw ODs (one column per wavelength, indexed by well).
        groups: SoftMax group tables following the plate (e.g. ``"Standards"``,
            ``"Unk_Dilution"``, ``"Control"``), with the ``Sample`` column
            carried forward to every well row.
    """

    def __init__(
        self,
        name: str,
        index: int,
        plate: pd.DataFrame,
        raw: pd.DataFrame,
        groups: dict[str, pd.DataFrame],
    ):
        self.name, self.index, self.plate, self.raw, self.groups = (name, index, plate, raw, groups)

    def __repr__(self) -> str:
        return f"SoftMaxPlate({self.name!r}, index={self.index}, groups={list(self.groups)})"


def read_softmax_export(
    path: str | Path, *, plate_id_prefix: str | None = None, subtract_reference: bool = True
) -> list[SoftMaxPlate]:
    """Parse a SoftMax Pro "Plate + Groups" text export (tab-delimited, often Latin-1).

    Each ``Plate:`` block holds one row of 96 well values per wavelength. With a
    dual-wavelength read (e.g. ``450 620``), the plate's OD is
    ``OD(primary) − OD(reference)`` when ``subtract_reference`` is set. The
    reference wavelength corrects for scratches, fingerprints, and plastic
    imperfections, which absorb at every wavelength. This matches what
    SoftMax reports as the well "Value".

    ``Group:`` blocks after a plate are returned as DataFrames. They are where
    SoftMax stores the plate map (which wells are which standard or sample),
    along with the operator's masked wells.

    Returns:
        One :class:`SoftMaxPlate` per plate block, in file order. Plate IDs are
        ``<prefix>_p<n>`` (prefix defaults to the file stem).
    """
    path = Path(path)
    text = _read_text(path).replace("\r\n", "\n").replace("\r", "\n")
    lines = [ln.rstrip("\t ") for ln in text.split("\n")]
    prefix = plate_id_prefix or path.stem
    plates: list[SoftMaxPlate] = []
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("Plate:"):
            head = line.split("\t")
            name = head[1] if len(head) > 1 else f"plate{len(plates) + 1}"
            wl_field = next((f for f in head if re.fullmatch(r"\d{3}( \d{3})*\s*", f)), "")
            wavelengths = wl_field.split() or ["od"]
            header = lines[i + 1].split("\t")
            wells = [normalize_well(h) for h in header]
            well_idx = [k for k, w in enumerate(wells) if w]
            reads: list[list[str]] = []
            j = i + 2
            while j < len(lines) and not lines[j].startswith("~End"):
                if lines[j].strip():
                    reads.append(lines[j].split("\t"))
                j += 1
            if not reads:
                raise PlateFormatError(path, [f"plate block at line {i + 1} has no data"])
            raw = pd.DataFrame(
                {
                    wl: [_parse_od(r[k] if k < len(r) else "")[0] for k in well_idx]
                    for wl, r in zip(wavelengths, reads, strict=False)
                },
                index=[wells[k] for k in well_idx],
            )
            od = raw.iloc[:, 0]
            if subtract_reference and raw.shape[1] > 1:
                od = od - raw.iloc[:, 1]
            pid = f"{prefix}_p{len(plates) + 1}"
            pairs = [(w, f"{v}" if np.isfinite(v) else "") for w, v in od.items()]
            plate = _pairs_to_plate(path, pid, pairs)
            plates.append(SoftMaxPlate(name, len(plates), plate, raw, {}))
            i = j + 1
            continue
        if line.startswith("Group:") and plates:
            gname = line.split(":", 1)[1].strip()
            header = lines[i + 1].split("\t")
            rows, j = [], i + 2
            while j < len(lines) and lines[j].strip() and not lines[j].startswith("~End"):
                rows.append(lines[j].split("\t"))
                j += 1
            width = len(header)
            df = pd.DataFrame(
                [r[:width] + [""] * (width - len(r)) for r in rows], columns=header, dtype=object
            )
            df = df.apply(lambda s: s.str.strip())
            if "Sample" in df:
                df["Sample"] = df["Sample"].replace("", np.nan).ffill()
            plates[-1].groups[gname] = df
            while j < len(lines) and not lines[j].startswith("~End"):
                j += 1
            i = j + 1
            continue
        i += 1
    if not plates:
        raise PlateFormatError(path, ["no 'Plate:' blocks found; is this a SoftMax Pro export?"])
    return plates
