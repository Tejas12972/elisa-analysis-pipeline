"""Tests for plate and layout parsing (immunoassay.io)."""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from immunoassay.io import (
    ALL_WELLS,
    COLUMNS,
    ROWS,
    LayoutError,
    MergeError,
    PlateFormatError,
    PlateWarning,
    load_plate,
    merge_plate_layout,
    normalize_well,
    read_layout,
    read_plate,
    read_softmax_export,
    subtract_blank,
)

REPO = Path(__file__).resolve().parents[1]


def _od(well: str) -> float:
    """Deterministic, distinct OD per well so parsed positions can be checked exactly."""
    r = ROWS.index(well[0])
    c = int(well[1:])
    return round(0.1 + r * 0.25 + c * 0.01, 4)


def grid_lines(delim: str = ",", *, fill=_od, labels: bool = True) -> list[str]:
    lines = []
    if labels:
        lines.append(delim.join([""] + [str(c) for c in COLUMNS]))
    for r in ROWS:
        cells = [str(fill(f"{r}{c}")) for c in COLUMNS]
        lines.append(delim.join(([r] if labels else []) + cells))
    return lines


def write(path: Path, lines: list[str]) -> Path:
    path.write_text("\n".join(lines) + "\n")
    return path


def assert_plate_values(plate: pd.DataFrame) -> None:
    assert len(plate) == 96
    assert list(plate["well"]) == list(ALL_WELLS)
    np.testing.assert_allclose(plate["od"], [_od(w) for w in ALL_WELLS])
    assert (plate["od_status"] == "ok").all()


# --------------------------------------------------------------------------- #
# normalize_well
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        ("A1", "A1"),
        ("a01", "A1"),
        (" H 12 ", "H12"),
        ("B007", "B7"),
        ("c10", "C10"),
    ],
)
def test_normalize_well_accepts_common_spellings(raw, expected):
    assert normalize_well(raw) == expected


@pytest.mark.parametrize("raw", ["I1", "A0", "A13", "AA1", "", "1A", None, 3])
def test_normalize_well_rejects_non_wells(raw):
    assert normalize_well(raw) is None


# --------------------------------------------------------------------------- #
# read_plate: grid format
# --------------------------------------------------------------------------- #
def test_grid_csv(tmp_path):
    plate = read_plate(write(tmp_path / "p.csv", grid_lines()))
    assert_plate_values(plate)
    assert set(plate.columns) == {"plate_id", "well", "row", "column", "od", "od_status"}
    assert (plate["plate_id"] == "p").all()


def test_grid_with_instrument_metadata_and_trailing_rows(tmp_path):
    lines = [
        "Instrument: DemoReader 3000",
        "Date,2026-09-29",
        "",
        *grid_lines(),
        "",
        "Temperature,25.1",
    ]
    assert_plate_values(read_plate(write(tmp_path / "p.csv", lines)))


def test_grid_tab_delimited_and_offset_columns(tmp_path):
    lines = ["\t".join(["Plate 1"]), *["\t" + ln for ln in grid_lines("\t")]]
    assert_plate_values(read_plate(write(tmp_path / "p.txt", lines)))


def test_grid_semicolon_with_decimal_commas(tmp_path):
    lines = grid_lines(";", fill=lambda w: str(_od(w)).replace(".", ","))
    assert_plate_values(read_plate(write(tmp_path / "p.csv", lines)))


def test_bare_unlabelled_grid(tmp_path):
    assert_plate_values(read_plate(write(tmp_path / "p.csv", grid_lines(labels=False))))


def test_grid_excel(tmp_path):
    grid = pd.DataFrame(
        [[_od(f"{r}{c}") for c in COLUMNS] for r in ROWS], index=list(ROWS), columns=list(COLUMNS)
    )
    path = tmp_path / "p.xlsx"
    with pd.ExcelWriter(path) as xw:
        pd.DataFrame([["Reader export"]]).to_excel(
            xw, sheet_name="Results", header=False, index=False
        )
        grid.to_excel(xw, sheet_name="Results", startrow=3)
    assert_plate_values(read_plate(path, sheet="Results"))


def test_multiple_blocks_require_explicit_choice(tmp_path):
    lines = ["450 nm", *grid_lines(), "", "570 nm", *grid_lines(fill=lambda w: 0.04)]
    path = write(tmp_path / "p.csv", lines)
    with pytest.raises(PlateFormatError, match="found 2 8x12 grids"):
        read_plate(path)
    assert_plate_values(read_plate(path, block=0))
    assert (read_plate(path, block=1)["od"] == 0.04).all()
    with pytest.raises(PlateFormatError, match="only 2 found"):
        read_plate(path, block=2)


def test_overflow_and_empty_cells_are_statused_not_errors(tmp_path):
    special = {"A1": "OVRFLW", "B2": ">4.000", "C3": ""}
    lines = grid_lines(fill=lambda w: special.get(w, _od(w)))
    plate = read_plate(write(tmp_path / "p.csv", lines)).set_index("well")
    assert plate.loc["A1", "od_status"] == "overflow"
    assert plate.loc["B2", "od_status"] == "overflow"
    assert plate.loc["C3", "od_status"] == "missing"
    assert plate.loc[["A1", "B2", "C3"], "od"].isna().all()
    assert (plate.drop(["A1", "B2", "C3"])["od_status"] == "ok").all()


def test_unparseable_values_are_all_reported(tmp_path):
    bad = {"A1": "abc", "D7": "1.2.3", "H12": "n/a?"}
    lines = grid_lines(fill=lambda w: bad.get(w, _od(w)))
    with pytest.raises(PlateFormatError) as exc:
        read_plate(write(tmp_path / "p.csv", lines))
    assert len(exc.value.problems) == 3
    for well in bad:
        assert any(well in p for p in exc.value.problems)


def test_negative_raw_od_warns_but_keeps_value(tmp_path):
    lines = grid_lines(fill=lambda w: -0.01 if w == "E5" else _od(w))
    with pytest.warns(PlateWarning, match="negative raw OD"):
        plate = read_plate(write(tmp_path / "p.csv", lines))
    assert plate.set_index("well").loc["E5", "od"] == -0.01


def test_unrecognised_file_is_rejected(tmp_path):
    with pytest.raises(PlateFormatError, match="no plate data found"):
        read_plate(write(tmp_path / "p.csv", ["foo,bar", "1,2"]))


def test_all_empty_plate_is_rejected(tmp_path):
    with pytest.raises(PlateFormatError, match="no OD readings"):
        read_plate(write(tmp_path / "p.csv", grid_lines(fill=lambda w: "")))


def test_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_plate(tmp_path / "nope.csv")


# --------------------------------------------------------------------------- #
# read_plate: long format
# --------------------------------------------------------------------------- #
def test_long_csv_with_zero_padded_wells_and_extra_columns(tmp_path):
    lines = ["Plate,Well,Sample,OD450"] + [
        f"P1,{w[0]}{int(w[1:]):02d},x,{_od(w)}" for w in reversed(ALL_WELLS)
    ]
    assert_plate_values(read_plate(write(tmp_path / "p.csv", lines)))


def test_long_excel(tmp_path):
    path = tmp_path / "p.xlsx"
    pd.DataFrame(
        {"Well Position": list(ALL_WELLS), "Absorbance": [_od(w) for w in ALL_WELLS]}
    ).to_excel(path, index=False)
    assert_plate_values(read_plate(path))


def test_long_partial_plate_fills_missing(tmp_path):
    lines = ["well,od", "A1,0.5", "A2,0.6"]
    plate = read_plate(write(tmp_path / "p.csv", lines))
    assert len(plate) == 96
    assert (plate["od_status"] == "ok").sum() == 2


def test_long_explicit_column_names(tmp_path):
    lines = ["loc,reading_a,reading_b"] + [f"{w},9,{_od(w)}" for w in ALL_WELLS]
    plate = read_plate(write(tmp_path / "p.csv", lines), well_col="loc", value_col="reading_b")
    assert_plate_values(plate)


def test_long_duplicate_and_invalid_wells(tmp_path):
    lines = ["well,od", "A1,0.5", "a01,0.6", "Z9,0.1", "B1,0.2"]
    with pytest.raises(PlateFormatError) as exc:
        read_plate(write(tmp_path / "p.csv", lines))
    msgs = "\n".join(exc.value.problems)
    assert "A1 appears more than once" in msgs
    assert "'Z9'" in msgs


# --------------------------------------------------------------------------- #
# read_layout
# --------------------------------------------------------------------------- #
LAYOUT_HEADER = "well,type,sample_id,concentration,dilution,group"


def layout_file(tmp_path: Path, rows: list[str], header: str = LAYOUT_HEADER) -> Path:
    return write(tmp_path / "layout.csv", [header, *rows])


def test_layout_basic_and_defaults(tmp_path):
    layout = read_layout(
        layout_file(
            tmp_path,
            [
                "B1,standard,,100,,",
                "A1,std,S1,200,,",
                "H1,blank,,,,",
                "C3,unknown,P01,,4,vehicle",
                "C4,QC,CTRL,,,",
            ],
        )
    )
    assert list(layout["well"]) == ["A1", "B1", "C3", "C4", "H1"]  # plate order
    rec = layout.set_index("well")
    assert rec.loc["A1", "type"] == "standard"
    assert rec.loc["B1", "sample_id"] == "STD_100"  # generated from the level
    assert rec.loc["H1", "sample_id"] == "BLANK"
    assert rec.loc["C3", "type"] == "sample"
    assert rec.loc["C3", "dilution"] == 4.0
    assert rec.loc["C4", "dilution"] == 1.0
    assert rec.loc["C4", "type"] == "control"
    assert rec.loc["C3", "row"] == "C"
    assert rec.loc["C3", "column"] == 3


def test_layout_column_synonyms_and_extra_columns(tmp_path):
    layout = read_layout(
        layout_file(
            tmp_path,
            ["A1,Standard,S1,50,1,,IL-6", "A2,Sample,X,,2,Treated,IL-6"],
            header="Well,Well Type,Sample,Conc,Dilution Factor,Treatment,Analyte",
        )
    )
    assert {"type", "sample_id", "concentration", "dilution", "group", "analyte"} <= set(
        layout.columns
    )
    assert layout.set_index("well").loc["A2", "group"] == "Treated"


def test_layout_excel(tmp_path):
    path = tmp_path / "layout.xlsx"
    pd.DataFrame(
        {"well": ["A1", "H1"], "type": ["standard", "blank"], "concentration": [10, None]}
    ).to_excel(path, index=False)
    assert len(read_layout(path)) == 2


def test_layout_missing_required_column(tmp_path):
    with pytest.raises(LayoutError, match="missing required column"):
        read_layout(layout_file(tmp_path, ["A1,S1"], header="well,sample_id"))


@pytest.mark.parametrize(
    ("row", "fragment"),
    [
        ("Q1,sample,S,,,", "invalid well name"),
        ("A1,reagent,S,,,", "unknown type"),
        ("A1,standard,S,,,", "missing its nominal concentration"),
        ("A1,standard,S,0,,", "type zero standards as 'blank'"),
        ("A1,standard,S,abc,,", "not a number"),
        ("A1,sample,S,12,,", "only standards may have a concentration"),
        ("A1,sample,,,,", "need a sample_id"),
        ("A1,control,,,,", "need a sample_id"),
        ("A1,sample,S,,0.25,", "fold-dilution >= 1"),
        ("A1,sample,S,,x,", "dilution 'x' is not a number"),
    ],
)
def test_layout_row_validation(tmp_path, row, fragment):
    with pytest.raises(LayoutError, match=fragment):
        read_layout(layout_file(tmp_path, [row]))


def test_layout_reports_every_problem_at_once(tmp_path):
    with pytest.raises(LayoutError) as exc:
        read_layout(
            layout_file(
                tmp_path,
                [
                    "A1,sample,,,,",
                    "A2,foo,S,,,",
                    "A3,standard,S,,,",
                    "A1,blank,,,,",
                ],
            )
        )
    assert len(exc.value.problems) == 4


@pytest.mark.parametrize(
    ("rows", "fragment"),
    [
        (["A1,sample,S1,,,a", "A2,control,S1,,,"], "several types"),
        (["A1,sample,S1,,,a", "A2,sample,S1,,,b"], "several groups"),
        (["A1,standard,STD1,100,,", "A2,standard,STD1,50,,"], "conflicting concentrations"),
    ],
)
def test_layout_id_consistency(tmp_path, rows, fragment):
    with pytest.raises(LayoutError, match=fragment):
        read_layout(layout_file(tmp_path, rows))


def test_layout_replicates_and_multi_dilution_are_allowed(tmp_path):
    layout = read_layout(
        layout_file(
            tmp_path,
            [
                "A1,sample,S1,,2,a",
                "A2,sample,S1,,2,a",
                "B1,sample,S1,,8,a",
                "B2,sample,S1,,8,",
            ],
        )
    )
    assert len(layout) == 4


def test_layout_empty(tmp_path):
    with pytest.raises(LayoutError, match="no wells"):
        read_layout(layout_file(tmp_path, []))


# --------------------------------------------------------------------------- #
# merge + blank subtraction
# --------------------------------------------------------------------------- #
@pytest.fixture
def plate(tmp_path) -> pd.DataFrame:
    return read_plate(write(tmp_path / "p.csv", grid_lines()), plate_id="P1")


@pytest.fixture
def small_layout(tmp_path) -> pd.DataFrame:
    return read_layout(
        layout_file(
            tmp_path,
            [
                "A1,standard,,100,,",
                "A2,standard,,100,,",
                "H1,blank,,,,",
                "H2,blank,,,,",
                "C5,sample,S1,,2,ctrl",
            ],
        )
    )


def test_merge_keeps_layout_wells_and_warns_on_unassigned(plate, small_layout):
    with pytest.warns(PlateWarning, match="91 well"):
        merged = merge_plate_layout(plate, small_layout)
    assert list(merged["well"]) == ["A1", "A2", "C5", "H1", "H2"]
    assert merged.set_index("well").loc["C5", "od"] == _od("C5")
    assert (merged["plate_id"] == "P1").all()


def test_merge_errors_on_layout_wells_without_reading(tmp_path, small_layout):
    lines = grid_lines(fill=lambda w: "" if w in ("A2", "C5") else _od(w))
    plate = read_plate(write(tmp_path / "p.csv", lines))
    with pytest.raises(MergeError, match="A2, C5"):
        merge_plate_layout(plate, small_layout)


def test_merge_keeps_overflow_wells(tmp_path, small_layout):
    lines = grid_lines(fill=lambda w: "OVRFLW" if w == "A1" else _od(w))
    plate = read_plate(write(tmp_path / "p.csv", lines))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PlateWarning)
        merged = merge_plate_layout(plate, small_layout).set_index("well")
    assert merged.loc["A1", "od_status"] == "overflow"


def test_subtract_blank_math(plate, small_layout):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PlateWarning)
        out = subtract_blank(merge_plate_layout(plate, small_layout)).set_index("well")
    blanks = np.array([_od("H1"), _od("H2")])
    assert out["blank_mean"].iloc[0] == pytest.approx(blanks.mean())
    assert out["blank_sd"].iloc[0] == pytest.approx(blanks.std(ddof=1))
    assert (out["n_blanks"] == 2).all()
    assert out.loc["A1", "od_net"] == pytest.approx(_od("A1") - blanks.mean())
    # Negative net ODs are expected (A1 reads below the H-row blanks here) and kept.
    assert out.loc["A1", "od_net"] < 0


def test_subtract_blank_is_per_plate(plate, small_layout):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PlateWarning)
        m1 = merge_plate_layout(plate, small_layout)
    m2 = m1.assign(plate_id="P2", od=m1["od"] + 1.0)
    out = subtract_blank(pd.concat([m1, m2], ignore_index=True))
    np.testing.assert_allclose(
        out.loc[out.plate_id == "P1", "od_net"], out.loc[out.plate_id == "P2", "od_net"]
    )


def test_subtract_blank_requires_readable_blank(plate, tmp_path):
    layout = read_layout(layout_file(tmp_path, ["A1,standard,,100,,"]))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PlateWarning)
        merged = merge_plate_layout(plate, layout)
    with pytest.raises(MergeError, match="no readable blank"):
        subtract_blank(merged)


def test_single_blank_warns(plate, tmp_path):
    layout = read_layout(layout_file(tmp_path, ["A1,standard,,100,,", "H1,blank,,,,"]))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PlateWarning)
        merged = merge_plate_layout(plate, layout)
    with pytest.warns(PlateWarning, match="single blank"):
        out = subtract_blank(merged)
    assert out["blank_sd"].isna().all()


# --------------------------------------------------------------------------- #
# load_plate on the shipped demo files
# --------------------------------------------------------------------------- #
def test_load_plate_demo_files():
    demo = REPO / "data/simulated/demo"
    df = load_plate(demo / "plates/IL6_P1.csv", demo / "layouts/IL6_P1_layout.csv")
    assert df["plate_id"].iloc[0] == "IL6_P1"
    assert set(df["type"]) == {"standard", "blank", "control", "sample"}
    # Standards must increase in signal with concentration.
    std = df[df.type == "standard"].groupby("concentration")["od_net"].mean()
    assert std.is_monotonic_increasing


# --------------------------------------------------------------------------- #
# SoftMax Pro exports (real files from the ELISAtools package)
# --------------------------------------------------------------------------- #
SOFTMAX = REPO / "data/raw/elisatools_feng2019"


def test_softmax_export_plates_and_groups():
    plates = read_softmax_export(SOFTMAX / "Assay_3_and_4.txt")
    assert len(plates) == 2
    p = plates[0]
    assert list(p.raw.columns) == ["450", "620"]
    a1 = p.plate.set_index("well").loc["A1", "od"]
    assert a1 == pytest.approx(p.raw.loc["A1", "450"] - p.raw.loc["A1", "620"])
    assert {"Standards", "Unk_Dilution", "Control"} <= set(p.groups)
    std = p.groups["Standards"]
    assert std["Sample"].notna().all()  # carried forward
    assert p.plate["plate_id"].tolist()[0] == "Assay_3_and_4_p1"


def test_softmax_masked_wells_visible_in_groups():
    p = read_softmax_export(SOFTMAX / "Assay_2.txt")[0]
    std = p.groups["Standards"]
    assert set(std.loc[std["Value"] == "Masked", "Wells"]) == {"C1", "E1", "F1", "G2"}


def test_softmax_without_reference_subtraction():
    p = read_softmax_export(SOFTMAX / "Assay_2.txt", subtract_reference=False)[0]
    assert p.plate.set_index("well").loc["A1", "od"] == pytest.approx(2.2416)


def test_softmax_rejects_other_files(tmp_path):
    with pytest.raises(PlateFormatError, match="no 'Plate:' blocks"):
        read_softmax_export(write(tmp_path / "x.txt", ["hello"]))
