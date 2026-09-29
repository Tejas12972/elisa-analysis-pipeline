"""Synthetic ELISA plate generator with known ground truth.

A pipeline that is only tested on real data can only be checked for
*plausibility*, because nobody knows the true concentrations. Simulated
plates fix that. Every sample's true concentration is written to disk, so
the end-to-end tests can assert that the pipeline *recovers* it.

The simulator deliberately reproduces the failure modes the pipeline must
handle.

* **Heteroscedastic noise.** Each well's OD is drawn with
  ``SD = sqrt(sd_floor² + (cv · OD)²)``. There is a constant-CV component
  (pipetting, reagent volume) plus a small additive floor (reader
  electronics). This is why the curve fit is weighted.
* **Edge effects.** Outer wells evaporate faster and warm up sooner during
  incubation, so they read systematically high. ``edge_effect=0.1`` makes
  the 36 perimeter wells read 10% high.
* **Two kinds of plate-to-plate variation**, which is the key subtlety:

  - ``od_scale_cv``: the whole plate reads brighter or darker, e.g. from
    substrate development time or temperature. It multiplies *every* well,
    standards included, so the plate's own standard curve **absorbs** it.
    Back-calculated concentrations are unaffected.
  - ``calibrator_bias_cv``: the standard stock itself is off on this
    plate, e.g. from a reconstitution or serial-dilution error. The
    standards then contain ``k × nominal`` analyte, so every sample on the
    plate back-calculates at ``true / k``. The standard curve **cannot**
    detect this, because it is the reference that is wrong. Only samples
    shared between plates (bridge controls) can. That is the job of
    :mod:`immunoassay.normalize`.

* **Injected outliers.** Random wells get a gross error (bubble, splash,
  missed wash), which the QC outlier detector should catch.
* **High-dose hook effect (optional).** At extreme concentrations free
  analyte saturates *both* antibodies separately, so fewer sandwiches form
  and the signal paradoxically falls. Samples above ``hook_onset`` have
  their signal reduced. Diluting them further recovers the true value,
  which is what the dilution-linearity check detects.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from immunoassay.curves import four_pl
from immunoassay.io import ALL_WELLS, COLUMNS, ROWS

_EDGE_WELLS = frozenset(w for w in ALL_WELLS if w[0] in "AH" or int(w[1:]) in (1, 12))


@dataclass(frozen=True)
class AnalyteSpec:
    """One analyte's assay characteristics and biology.

    Attributes:
        name: Analyte name, e.g. ``"IL-6"``.
        a, b, c, d: True 4PL parameters in net OD (``a`` floor, ``d``
            saturation, ``c`` EC50 in ``units``, ``b`` Hill slope).
        top_standard: Highest standard concentration. Lower levels are made
            by serial dilution by ``standard_dilution``.
        n_standards: Number of non-zero standard levels.
        group_geo_means: Geometric-mean *neat* sample concentration per
            treatment group. Cytokine levels are roughly log-normal, so
            groups differ by a *fold change*.
        geo_sd: Geometric SD of the between-subject variation.
        controls: Neat concentration of each bridge/QC control sample.
        units: Concentration units (used in reports).
    """

    name: str = "IL-6"
    a: float = 0.02
    b: float = 1.1
    c: float = 350.0
    d: float = 3.2
    top_standard: float = 1000.0
    standard_dilution: float = 2.0
    n_standards: int = 7
    group_geo_means: dict[str, float] = field(
        default_factory=lambda: {"vehicle": 80.0, "LPS": 400.0, "LPS+drug": 180.0}
    )
    geo_sd: float = 1.6
    controls: dict[str, float] = field(
        default_factory=lambda: {"CTRL_LO": 40.0, "CTRL_MID": 160.0, "CTRL_HI": 600.0}
    )
    units: str = "pg/mL"

    @property
    def standard_levels(self) -> np.ndarray:
        return self.top_standard / self.standard_dilution ** np.arange(self.n_standards)


@dataclass(frozen=True)
class NoiseSpec:
    """Measurement-error model. See the module docstring for the meaning of each term."""

    cv: float = 0.05
    sd_floor: float = 0.005
    background: float = 0.06
    edge_effect: float = 0.0
    od_scale_cv: float = 0.10
    calibrator_bias_cv: float = 0.0
    outlier_rate: float = 0.0
    outlier_magnitude: float = 0.5
    hook_onset: float | None = None
    reader_max_od: float = 4.0


@dataclass
class SimulatedStudy:
    """Everything the simulator produces.

    Attributes:
        plates: One ``(plate_id, analyte, plate_grid, layout)`` tuple per
            plate. ``plate_grid`` is an 8 x 12 DataFrame of raw OD strings,
            as a reader would export them.
        truth: True neat concentration per ``(analyte, sample_id)``, with
            groups.
        plate_factors: True calibrator bias ``k`` and OD scale per plate.
        injected_outliers: ``(plate_id, well)`` of every gross-error well.
        units: Units per analyte.
    """

    plates: list[tuple[str, str, pd.DataFrame, pd.DataFrame]]
    truth: pd.DataFrame
    plate_factors: pd.DataFrame
    injected_outliers: pd.DataFrame
    units: dict[str, str]

    def control_ranges(self, tolerance: float = 0.25) -> pd.DataFrame:
        """Acceptance range per control: nominal ± ``tolerance``.

        This mimics the ranges printed on a kit's certificate of analysis.
        """
        ctrl = self.truth[self.truth["group"].isna()]
        return pd.DataFrame(
            {
                "analyte": ctrl["analyte"],
                "sample_id": ctrl["sample_id"],
                "nominal": ctrl["true_conc"],
                "low": ctrl["true_conc"] * (1 - tolerance),
                "high": ctrl["true_conc"] * (1 + tolerance),
            }
        ).reset_index(drop=True)

    def write(self, outdir: str | Path) -> Path:
        """Write grid-format plate CSVs, layouts, truth, and a manifest to ``outdir``.

        Returns:
            The path of ``manifest.csv``, which :func:`immunoassay.pipeline.run_pipeline`
            accepts directly. ``control_ranges.csv`` is written next to it.
        """
        outdir = Path(outdir)
        (outdir / "plates").mkdir(parents=True, exist_ok=True)
        (outdir / "layouts").mkdir(parents=True, exist_ok=True)
        manifest = []
        for plate_id, analyte, grid, layout in self.plates:
            plate_path = outdir / "plates" / f"{plate_id}.csv"
            layout_path = outdir / "layouts" / f"{plate_id}_layout.csv"
            lines = [
                f"Simulated plate export,{plate_id}",
                f"Analyte,{analyte}",
                "Wavelength,450 nm",
                "",
                "," + ",".join(str(c) for c in COLUMNS),
            ]
            lines += [f"{r}," + ",".join(grid.loc[r].tolist()) for r in ROWS]
            plate_path.write_text("\n".join(lines) + "\n")
            layout.to_csv(layout_path, index=False)
            manifest.append(
                {
                    "plate_id": plate_id,
                    "analyte": analyte,
                    "plate_file": f"plates/{plate_path.name}",
                    "layout_file": f"layouts/{layout_path.name}",
                    "units": self.units[analyte],
                }
            )
        pd.DataFrame(manifest).to_csv(outdir / "manifest.csv", index=False)
        self.control_ranges().to_csv(outdir / "control_ranges.csv", index=False)
        self.truth.to_csv(outdir / "truth.csv", index=False)
        self.plate_factors.to_csv(outdir / "plate_factors.csv", index=False)
        self.injected_outliers.to_csv(outdir / "injected_outliers.csv", index=False)
        return outdir / "manifest.csv"


def _plate_design(
    spec: AnalyteSpec,
    sample_slots: list[tuple[str, str, float]],
) -> list[dict[str, object]]:
    """Assign wells column-pair by column-pair (duplicates side by side).

    Columns 1-2: standards top-down in rows A-G, blanks in row H.
    Then bridge controls, then samples, each in horizontally adjacent
    duplicate wells.
    """
    rows: list[dict[str, object]] = []
    for i, level in enumerate(spec.standard_levels):
        for col in (1, 2):
            rows.append(
                {
                    "well": f"{ROWS[i]}{col}",
                    "type": "standard",
                    "sample_id": f"STD{i + 1}",
                    "concentration": float(level),
                    "dilution": 1.0,
                    "group": None,
                }
            )
    for col in (1, 2):
        rows.append(
            {
                "well": f"H{col}",
                "type": "blank",
                "sample_id": "BLANK",
                "concentration": np.nan,
                "dilution": 1.0,
                "group": None,
            }
        )
    pair_wells = [(f"{r}{c}", f"{r}{c + 1}") for c in range(3, 12, 2) for r in ROWS]
    entries: list[tuple[str, str, float, str | None]] = [
        ("control", cid, 1.0, None) for cid in spec.controls
    ] + [("sample", sid, dil, grp) for sid, grp, dil in sample_slots]
    if len(entries) > len(pair_wells):
        raise ValueError(
            f"{len(entries)} duplicate pairs do not fit on one plate ({len(pair_wells)} free)"
        )
    for (w1, w2), (typ, sid, dil, grp) in zip(pair_wells, entries, strict=False):
        for w in (w1, w2):
            rows.append(
                {
                    "well": w,
                    "type": typ,
                    "sample_id": sid,
                    "concentration": np.nan,
                    "dilution": dil,
                    "group": grp,
                }
            )
    return rows


def simulate_study(
    analytes: list[AnalyteSpec] | None = None,
    *,
    n_per_group: int = 12,
    sample_dilution: float = 2.0,
    linearity_dilution: float = 8.0,
    n_linearity_per_plate: int = 4,
    noise: NoiseSpec | None = None,
    seed: int = 0,
) -> SimulatedStudy:
    """Simulate a multi-plate, multi-analyte study with known truth.

    The same subjects (``n_per_group`` per treatment group) are measured for
    every analyte. Each analyte gets as many plates as it needs. Subjects
    are shuffled across plates so treatment group is not confounded with
    plate (randomised plate assignment is standard practice). The first
    ``n_linearity_per_plate`` samples on each plate are also run at
    ``linearity_dilution`` to exercise dilution-linearity checks.

    Returns:
        A :class:`SimulatedStudy`. Call ``.write(outdir)`` to get files the
        pipeline can read.
    """
    analytes = analytes or [AnalyteSpec()]
    noise = noise or NoiseSpec()
    rng = np.random.default_rng(seed)

    groups = list(analytes[0].group_geo_means)
    subjects = [
        (f"S{k + 1:03d}", grp)
        for k, grp in enumerate(g for g in groups for _ in range(n_per_group))
    ]

    truth_rows, plates, factor_rows, outlier_rows = [], [], [], []
    for spec in analytes:
        if list(spec.group_geo_means) != groups:
            raise ValueError("all analytes must define the same treatment groups")
        true_conc = {
            sid: float(np.exp(rng.normal(np.log(spec.group_geo_means[grp]), np.log(spec.geo_sd))))
            for sid, grp in subjects
        }
        truth_rows += [
            {"analyte": spec.name, "sample_id": sid, "group": grp, "true_conc": true_conc[sid]}
            for sid, grp in subjects
        ]
        truth_rows += [
            {"analyte": spec.name, "sample_id": cid, "group": None, "true_conc": conc}
            for cid, conc in spec.controls.items()
        ]

        capacity = 40 - len(spec.controls)
        order = [subjects[i] for i in rng.permutation(len(subjects))]
        per_plate = capacity - n_linearity_per_plate
        n_plates = -(-len(order) // per_plate)  # ceil: balance subjects across plates
        chunks = [list(c) for c in np.array_split(np.array(order, dtype=object), n_plates)]
        chunks = [[tuple(x) for x in c] for c in chunks]
        for p, chunk in enumerate(chunks):
            ascii_name = re.sub(r"[^A-Za-z0-9]", "", spec.name.replace("α", "a"))
            plate_id = f"{ascii_name}_P{p + 1}"  # file-name safe, e.g. TNFa_P1
            slots = [(sid, grp, sample_dilution) for sid, grp in chunk]
            slots += [(sid, grp, linearity_dilution) for sid, grp in chunk[:n_linearity_per_plate]]
            design = _plate_design(spec, slots)

            k = float(np.exp(rng.normal(0, noise.calibrator_bias_cv)))
            scale = float(np.exp(rng.normal(0, noise.od_scale_cv)))
            factor_rows.append(
                {
                    "plate_id": plate_id,
                    "analyte": spec.name,
                    "calibrator_bias": k,
                    "od_scale": scale,
                }
            )

            grid = pd.DataFrame("", index=list(ROWS), columns=list(COLUMNS), dtype=object)
            for well_row in design:
                well = str(well_row["well"])
                typ = well_row["type"]
                if typ == "standard":
                    in_well = k * float(well_row["concentration"])  # type: ignore[arg-type]
                elif typ == "blank":
                    in_well = 0.0
                else:
                    neat = (
                        spec.controls[str(well_row["sample_id"])]
                        if typ == "control"
                        else true_conc[str(well_row["sample_id"])]
                    )
                    in_well = neat / float(well_row["dilution"])  # type: ignore[arg-type]
                signal = float(four_pl(in_well, spec.a, spec.b, spec.c, spec.d))
                # Kits place their top standard below the hook, so only samples hook.
                hooked = typ in ("sample", "control") and noise.hook_onset is not None
                if hooked and in_well > noise.hook_onset:
                    signal *= (noise.hook_onset / in_well) ** 0.5
                od = (signal + noise.background) * scale
                if noise.edge_effect and well in _EDGE_WELLS:
                    od *= 1 + noise.edge_effect
                od += rng.normal(0, np.hypot(noise.sd_floor, noise.cv * od))
                if noise.outlier_rate and typ != "blank" and rng.random() < noise.outlier_rate:
                    od *= 1 + rng.choice([-1, 1]) * noise.outlier_magnitude
                    outlier_rows.append({"plate_id": plate_id, "well": well})
                grid.loc[well[0], int(well[1:])] = (
                    "OVRFLW" if od >= noise.reader_max_od else f"{od:.4f}"
                )

            layout = pd.DataFrame(design)
            layout["concentration"] = layout["concentration"].map(
                lambda v: "" if pd.isna(v) else f"{v:g}"
            )
            layout["dilution"] = layout["dilution"].map(lambda v: f"{v:g}")
            layout["group"] = layout["group"].fillna("")
            layout.loc[layout["type"] == "blank", "sample_id"] = ""
            plates.append((plate_id, spec.name, grid, layout))

    return SimulatedStudy(
        plates=plates,
        truth=pd.DataFrame(truth_rows),
        plate_factors=pd.DataFrame(factor_rows),
        injected_outliers=pd.DataFrame(outlier_rows, columns=["plate_id", "well"]),
        units={s.name: s.units for s in analytes},
    )


def default_demo_study(seed: int = 17) -> SimulatedStudy:
    """The two-cytokine demo study used in the README and example report.

    IL-6 and TNF-α are measured in the same 36 subjects (vehicle, LPS, and
    LPS + anti-inflammatory drug; 12 each). The effects are chosen to be
    biologically plausible for an LPS challenge. IL-6 rises strongly and
    the drug blunts it partially. TNF-α rises less and the drug has little
    effect. The noise includes an 8% plate-level calibrator bias, so
    inter-plate normalization has something to fix, plus edge effects and
    occasional outlier wells.
    """
    il6 = AnalyteSpec()
    tnf = AnalyteSpec(
        name="TNF-α",
        a=0.03,
        b=1.2,
        c=180.0,
        d=2.8,
        top_standard=500.0,
        group_geo_means={"vehicle": 30.0, "LPS": 110.0, "LPS+drug": 95.0},
        geo_sd=1.5,
        controls={"CTRL_LO": 20.0, "CTRL_MID": 70.0, "CTRL_HI": 250.0},
    )
    noise = NoiseSpec(
        cv=0.05,
        edge_effect=0.03,
        od_scale_cv=0.12,
        calibrator_bias_cv=0.08,
        outlier_rate=0.01,
        outlier_magnitude=0.6,
    )
    return simulate_study([il6, tnf], n_per_group=12, noise=noise, seed=seed)
