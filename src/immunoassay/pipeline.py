"""End-to-end orchestration: files in, QC'd concentrations and statistics out.

Stage order matters, and each step exists for a reason:

1. **Parse + blank-subtract** each plate (:mod:`immunoassay.io`).
2. **Replicate outliers** (n ≥ 3) are flagged *before* fitting, so a gross
   error cannot distort the curve.
3. **Fit** a 4PL/5PL per plate. Then **standard residual outliers**
   (duplicate standards that disagree, where one well is far off the curve)
   are flagged, and the curve is refitted without them.
4. **Back-calculate** every well, compute **% recovery** of the standards,
   then **LOD/LLOQ/ULOQ**.
5. **Replicate summaries** with CV and range flags. **Dilution correction**
   to neat concentrations. **Dilution linearity**. One **reportable result**
   per sample.
6. **Controls** are checked against expected ranges *before*
   normalization, because normalization uses the controls and would
   otherwise hide their failures. Then the **plate pass/fail** verdict.
7. **Inter-plate normalization** from bridge controls.
8. **Statistics** on the normalized values, excluding failed plates and
   handling censored samples explicitly.
"""

from __future__ import annotations

import json
import warnings
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from immunoassay import curves, dilution, normalize, qc
from immunoassay import io as pio
from immunoassay import stats as st


@dataclass
class PipelineConfig:
    """All tunable settings in one place (serialisable into the report)."""

    model: str = "auto"
    weighting: str = "1/y^2"
    include_zero_standard: bool = True
    qc: qc.QCConfig = field(default_factory=qc.QCConfig)
    stats: st.StatsConfig = field(default_factory=st.StatsConfig)
    normalize: bool = True
    linearity_max_deviation_pct: float = 20.0
    control_ranges: pd.DataFrame | None = None
    group_order: list[str] | None = None
    units: dict[str, str] = field(default_factory=dict)

    def as_dict(self) -> dict:
        return {
            "model": self.model,
            "weighting": self.weighting,
            "include_zero_standard": self.include_zero_standard,
            "qc": asdict(self.qc),
            "stats": asdict(self.stats),
            "normalize": self.normalize,
            "linearity_max_deviation_pct": self.linearity_max_deviation_pct,
            "group_order": self.group_order,
            "control_ranges": (
                None
                if self.control_ranges is None
                else self.control_ranges.to_dict(orient="records")
            ),
        }


@dataclass
class PlateInput:
    """One plate export + its layout."""

    plate_path: Path
    layout_path: Path
    plate_id: str | None = None
    analyte: str | None = None
    plate_kwargs: dict = field(default_factory=dict)


@dataclass
class PipelineResult:
    """Every intermediate table, so each step can be audited."""

    wells: pd.DataFrame
    fits: dict[tuple[str, str], curves.CurveFit]
    recovery: pd.DataFrame
    limits: pd.DataFrame
    samples: pd.DataFrame
    linearity: pd.DataFrame
    results: pd.DataFrame
    plate_factors: pd.DataFrame
    controls: pd.DataFrame
    plate_qc: pd.DataFrame
    stats: st.StatsResult
    bridge_cv: pd.DataFrame
    config: PipelineConfig
    warnings: list[str] = field(default_factory=list)

    def curve_table(self) -> pd.DataFrame:
        rows = [
            {
                "plate_id": p,
                "analyte": a,
                **f.summary(),
                "aicc_4pl": f.comparison.get("4pl", np.nan),
                "aicc_5pl": f.comparison.get("5pl", np.nan),
            }
            for (p, a), f in self.fits.items()
        ]
        return pd.DataFrame(rows)

    def write_tables(self, outdir: str | Path) -> dict[str, Path]:
        """Write every result table as CSV and return their paths."""
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        tables = {
            "wells": self.wells,
            "standard_curves": self.curve_table(),
            "standard_recovery": self.recovery,
            "limits": self.limits,
            "replicates": self.samples,
            "dilution_linearity": self.linearity,
            "results": self.results,
            "plate_factors": self.plate_factors,
            "controls": self.controls,
            "plate_qc": self.plate_qc,
            "stats_omnibus": self.stats.omnibus,
            "stats_pairwise": self.stats.pairwise,
            "stats_input": self.stats.data,
        }
        paths = {}
        for name, df in tables.items():
            paths[name] = outdir / f"{name}.csv"
            df.to_csv(paths[name], index=False)
        (outdir / "config.json").write_text(
            json.dumps(self.config.as_dict(), indent=2, default=str)
        )
        return paths


def read_manifest(path: str | Path) -> list[PlateInput]:
    """Read a manifest CSV: ``plate_file, layout_file`` [, ``plate_id``, ``analyte``].

    Relative paths are resolved against the manifest's own directory.
    """
    path = Path(path)
    m = pd.read_csv(path, dtype=str).fillna("")
    missing = {"plate_file", "layout_file"} - set(m.columns)
    if missing:
        raise pio.InputFileError(path, [f"manifest missing column(s): {', '.join(missing)}"])
    base = path.parent
    return [
        PlateInput(
            plate_path=(base / r["plate_file"]).resolve(),
            layout_path=(base / r["layout_file"]).resolve(),
            plate_id=r.get("plate_id") or None,
            analyte=r.get("analyte") or None,
        )
        for r in m.to_dict("records")
    ]


def load_wells(inputs: list[PlateInput]) -> pd.DataFrame:
    """Load and blank-subtract all plates into one tidy frame with an ``analyte`` column."""
    frames = []
    for inp in inputs:
        df = pio.load_plate(
            inp.plate_path, inp.layout_path, plate_id=inp.plate_id, plate_kwargs=inp.plate_kwargs
        )
        if inp.analyte:
            df["analyte"] = inp.analyte
        elif "analyte" not in df.columns:
            df["analyte"] = "analyte"
        frames.append(df)
    wells = pd.concat(frames, ignore_index=True)
    dup = wells.groupby("plate_id")["well"].apply(lambda w: w.duplicated().any())
    if dup.any():
        raise pio.MergeError(
            "manifest", [f"plate_id {p!r} is used by more than one file" for p in dup[dup].index]
        )
    return wells


def run_pipeline(
    inputs: list[PlateInput] | str | Path, config: PipelineConfig | None = None
) -> PipelineResult:
    """Run the full pipeline on a manifest path or a list of :class:`PlateInput`."""
    if isinstance(inputs, (str, Path)):
        inputs = read_manifest(inputs)
    return analyze_wells(load_wells(inputs), config)


def analyze_wells(wells: pd.DataFrame, config: PipelineConfig | None = None) -> PipelineResult:
    """Run every stage after parsing, on an already-loaded tidy well frame.

    Useful when real data needed custom parsing (see ``notebooks/``).
    """
    cfg = config or PipelineConfig()
    caught: list[str] = []
    with warnings.catch_warnings(record=True) as wlist:
        warnings.simplefilter("always")
        if "analyte" not in wells.columns:
            wells = wells.assign(analyte="analyte")

        wells = qc.flag_replicate_outliers(wells, cfg.qc)
        fits = curves.fit_standard_curves(
            wells, model=cfg.model, weighting=cfg.weighting, include_zero=cfg.include_zero_standard
        )
        flagged = qc.flag_standard_residual_outliers(wells, fits, cfg.qc)
        if flagged["outlier"].sum() > wells["outlier"].sum():
            wells = flagged
            fits = curves.fit_standard_curves(
                wells,
                model=cfg.model,
                weighting=cfg.weighting,
                include_zero=cfg.include_zero_standard,
            )

        wells = curves.back_calculate(wells, fits)
        recovery = curves.standard_recovery(wells)
        limits = qc.plate_limits(wells, fits, recovery, cfg.qc)

        samples = qc.summarize_replicates(wells, limits, cfg.qc)
        samples = dilution.apply_dilution(samples)
        linearity = dilution.dilution_linearity(samples, cfg.linearity_max_deviation_pct)
        results = dilution.reportable_results(samples, linearity)

        controls = qc.check_controls(samples, cfg.control_ranges)
        plate_qc = qc.plate_summary(fits, limits, controls, samples, cfg.qc)

        if cfg.normalize:
            factors = normalize.estimate_plate_factors(results)
        else:
            factors = (
                results[["plate_id", "analyte"]]
                .drop_duplicates()
                .assign(factor=1.0, log_factor_se=np.nan, n_bridges=0, bridge_ids="")
            )
        results = normalize.apply_plate_factors(results, factors)
        bcv = normalize.bridge_cv(results, "conc").merge(
            normalize.bridge_cv(results, "conc_norm"),
            on=["analyte", "sample_id"],
            suffixes=("_before", "_after"),
            how="outer",
        )

        data = st.prepare_values(results, plate_qc, cfg.stats)
        stats_result = st.compare_groups(data, cfg.stats, cfg.group_order)
        caught = [str(w.message) for w in wlist]

    return PipelineResult(
        wells=wells,
        fits=fits,
        recovery=recovery,
        limits=limits,
        samples=samples,
        linearity=linearity,
        results=results,
        plate_factors=factors,
        controls=controls,
        plate_qc=plate_qc,
        stats=stats_result,
        bridge_cv=bcv,
        config=cfg,
        warnings=caught,
    )
