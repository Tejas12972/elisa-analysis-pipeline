"""Streamlit UI: pick or upload plates, tune QC thresholds, review results, download.

Runs two ways:

* locally: ``streamlit run app/streamlit_app.py``
* in the browser with no server, via stlite on GitHub Pages
  (``scripts/build_site.py`` bundles this file with the package and data).
"""

from __future__ import annotations

import io
import json
import logging
import platform
import tempfile
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import streamlit as st

import immunoassay
from immunoassay import plotting as P
from immunoassay import qc
from immunoassay import stats as stx
from immunoassay.io import InputFileError
from immunoassay.pipeline import PipelineConfig, PlateInput, read_manifest, run_pipeline
from immunoassay.report import build_report
from immunoassay.simulate import default_demo_study

# In the browser (stlite) the root logger runs at DEBUG; keep library chatter out of the console.
for _name in ("matplotlib", "PIL"):
    logging.getLogger(_name).setLevel(logging.WARNING)

HERE = Path(__file__).resolve().parent
# Repo checkout: app/ sits one level below data/. Browser bundle: everything at the root.
ROOT = next((p for p in (HERE, HERE.parent) if (p / "data").is_dir()), HERE.parent)
REPO_URL = "https://github.com/Tejas12972/elisa-analysis-pipeline"

DEMO = "Demo study: LPS challenge (simulated, known truth)"
REAL_IL6 = "Real: mouse IL-6 after sleep fragmentation (Dryad, CC0)"
REAL_SOFTMAX = "Real: SoftMax Pro plates, two kit lots (ELISAtools, MIT)"
UPLOAD = "Upload my own plates"

DATASETS = {
    REAL_IL6: (ROOT / "data/processed/sleep_fragmentation_il6/manifest.csv", False, None),
    REAL_SOFTMAX: (ROOT / "data/processed/elisatools_feng2019/manifest.csv", True, None),
}

st.set_page_config(page_title="ELISA Analysis Pipeline", page_icon="🧪", layout="wide")
st.title("🧪 ELISA Analysis Pipeline")
st.markdown(
    "Raw 96-well plate-reader data → **weighted 4PL/5PL standard curves** → **QC** "
    "(LOD/LLOQ/ULOQ, %CV, recovery, outliers, controls) → **dilution & inter-plate "
    "correction** → **group statistics** → a downloadable report. "
    f"Source, tests and methods: [GitHub]({REPO_URL})."
)
with st.expander("New here? How to use this app", expanded=False):
    st.markdown(
        """
1. **Pick a dataset** in the sidebar. The *demo study* is simulated, so the true answers are
   known. The two *real* datasets are public ELISA data from a published study and from an
   R package.
2. **Change a QC threshold** (e.g. max %CV) and watch plates pass or fail and samples get
   flagged. Everything reruns instantly.
3. **Explore the tabs**: standard curves with residuals, flagged samples, dilution linearity,
   plate-to-plate normalization, and the statistics with effect sizes and confidence intervals.
4. **Download** the full HTML report or every result table as CSV.

Want to try your own data? Choose *Upload my own plates*. Sample files are provided there.
The whole app runs in *your browser* (Python compiled to WebAssembly), and nothing is uploaded
to a server.
"""
    )

# --------------------------------------------------------------------------- #
# Sidebar
# --------------------------------------------------------------------------- #
with st.sidebar:
    st.header("1 · Data")
    source = st.radio(
        "Dataset", [DEMO, REAL_IL6, REAL_SOFTMAX, UPLOAD], label_visibility="collapsed"
    )
    uploads, layouts, ranges_file = [], [], None
    analyte, units = "analyte", "pg/mL"
    if source == UPLOAD:
        uploads = st.file_uploader(
            "Plate exports (CSV/TSV/Excel, 8×12 grid or long format)", accept_multiple_files=True
        )
        layouts = st.file_uploader(
            "Layouts, one per plate (name them <plate>_layout.csv)", accept_multiple_files=True
        )
        ranges_file = st.file_uploader("Control ranges (optional: analyte,sample_id,low,high)")
        analyte = st.text_input("Analyte", "IL-6")
        units = st.text_input("Units", "pg/mL")

    st.header("2 · Curve fitting")
    model = st.selectbox(
        "Model", ["auto", "4pl", "5pl"], help="auto: 5PL only if it lowers AICc by more than 2"
    )
    weighting = st.selectbox(
        "Weighting", ["1/y^2", "1/y", "none"], help="1/y² assumes constant-CV noise (typical)"
    )

    st.header("3 · QC thresholds")
    max_cv = st.slider("Max replicate %CV", 5.0, 40.0, 20.0, 1.0)
    rec_lo, rec_hi = st.slider("Standard recovery window (%)", 50, 150, (80, 120))
    min_r2 = st.slider("Min R²", 0.90, 1.0, 0.98, 0.005)
    min_std = st.slider("Min fraction of standard levels passing", 0.5, 1.0, 0.75, 0.05)

    st.header("4 · Statistics")
    transform = st.selectbox("Scale", ["log10", "none"])
    censored = st.selectbox(
        "Below-LLOQ values",
        ["substitute", "exclude", "extrapolate"],
        help="substitute: LLOQ/√2 (default) · extrapolate: sensitivity analysis only",
    )
    harmonize = st.checkbox("Common LLOQ across plates", True)
    exclude_failed = st.checkbox("Exclude plates that fail QC", True)
    normalize = st.checkbox("Inter-plate normalization (bridge controls)", True)
    group_order = st.text_input("Group order (first = reference)", "")
    st.caption(
        f"immunoassay {immunoassay.__version__} · Python {platform.python_version()} · "
        f"pandas {pd.__version__} · numpy {np.__version__}"
    )


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def _save_uploads(tmp: Path) -> tuple[list[PlateInput], pd.DataFrame | None]:
    if not uploads or not layouts:
        raise InputFileError("upload", ["upload at least one plate file and its layout"])
    lay_by_stem = {Path(f.name).stem.lower().replace("_layout", ""): f for f in layouts}
    inputs = []
    for i, f in enumerate(uploads):
        plate_path = tmp / f.name
        plate_path.write_bytes(f.getvalue())
        stem = Path(f.name).stem.lower()
        lay = lay_by_stem.get(stem) or (layouts[i] if len(layouts) == len(uploads) else None)
        if lay is None:
            raise InputFileError(f.name, [f"no layout named {stem}_layout.* was uploaded"])
        lay_path = tmp / f"layout_{i}_{lay.name}"
        lay_path.write_bytes(lay.getvalue())
        inputs.append(PlateInput(plate_path, lay_path, analyte=analyte))
    ranges = pd.read_csv(ranges_file) if ranges_file is not None else None
    return inputs, ranges


def _config(ranges: pd.DataFrame | None, unit_map: dict[str, str]) -> PipelineConfig:
    return PipelineConfig(
        model=model,
        weighting=weighting,
        qc=qc.QCConfig(
            max_cv_pct=max_cv,
            recovery_low_pct=rec_lo,
            recovery_high_pct=rec_hi,
            min_r_squared=min_r2,
            min_std_pass_frac=min_std,
        ),
        stats=stx.StatsConfig(
            transform=transform,
            censored=censored,
            harmonize_lloq=harmonize,
            exclude_failed_plates=exclude_failed,
        ),
        normalize=normalize,
        control_ranges=ranges,
        group_order=[g.strip() for g in group_order.split(",") if g.strip()] or None,
        units=unit_map,
    )


@st.cache_resource(show_spinner="Running pipeline…", max_entries=16)
def _cached_run(source: str, settings: str, _inputs, _cfg):
    """Built-in datasets: rerun only when the dataset or a setting changes."""
    return run_pipeline(_inputs, _cfg)


@st.cache_resource(show_spinner="Simulating demo study…")
def _demo_study():
    """Simulate the demo study once and write it to a temp dir (also used for sample files)."""
    study = default_demo_study()
    manifest = study.write(Path(tempfile.mkdtemp()))
    return manifest, study.control_ranges(), study.units


def _units_from_manifest(path: Path) -> dict[str, str]:
    m = pd.read_csv(path, dtype=str)
    if {"analyte", "units"} <= set(m.columns):
        return dict(zip(m["analyte"], m["units"], strict=True))
    return {}


# --------------------------------------------------------------------------- #
# Run
# --------------------------------------------------------------------------- #
try:
    tmpdir = Path(tempfile.mkdtemp())
    if source == DEMO:
        manifest, ranges, unit_map = _demo_study()
        inputs = read_manifest(manifest)
        if not group_order:
            group_order = "vehicle,LPS,LPS+drug"
    elif source in DATASETS:
        manifest, use_norm, ranges = DATASETS[source]
        if not manifest.exists():
            st.error(f"Missing {manifest.relative_to(ROOT)}. Run the scripts in `scripts/` first.")
            st.stop()
        inputs, unit_map = read_manifest(manifest), _units_from_manifest(manifest)
        normalize = normalize and use_norm
    else:
        manifest_dir = _demo_study()[0].parent
        st.info(
            "Upload plate exports and matching layouts in the sidebar. No files handy? "
            "Download a sample plate and layout from the simulated demo study, then upload them."
        )
        c1, c2, _ = st.columns([1, 1, 2])
        c1.download_button(
            "Sample plate (IL6_P1.csv)",
            (manifest_dir / "plates/IL6_P1.csv").read_bytes(),
            "IL6_P1.csv",
            "text/csv",
        )
        c2.download_button(
            "Sample layout (IL6_P1_layout.csv)",
            (manifest_dir / "layouts/IL6_P1_layout.csv").read_bytes(),
            "IL6_P1_layout.csv",
            "text/csv",
        )
        if not uploads:
            st.stop()
        inputs, ranges = _save_uploads(tmpdir)
        unit_map = {analyte: units}
    cfg = _config(ranges, unit_map)
    if source == UPLOAD:
        result = run_pipeline(inputs, cfg)
    else:
        result = _cached_run(source, json.dumps(cfg.as_dict(), default=str), inputs, cfg)
except InputFileError as exc:
    st.error(f"**{exc.source}**: the file could not be used as-is.")
    for p in exc.problems:
        st.write(f"- {p}")
    st.stop()
except RuntimeError as exc:
    st.error(str(exc))
    st.stop()

# --------------------------------------------------------------------------- #
# Results
# --------------------------------------------------------------------------- #
q = result.plate_qc
c1, c2, c3, c4 = st.columns(4)
c1.metric("Plates passing QC", f"{int(q.passed.sum())}/{len(q)}")
c2.metric("Samples", result.results.query("type == 'sample'").sample_id.nunique())
c3.metric("Flagged replicate groups", int((result.samples["flags"] != "").sum()))
c4.metric("Wells excluded as outliers", int(result.wells.outlier.sum()))

tabs = st.tabs(
    [
        "Plate QC",
        "Standard curves",
        "Samples & flags",
        "Dilution",
        "Normalization",
        "Statistics",
        "Download",
    ]
)
with tabs[0]:
    st.caption(
        "A plate passes when R² ≥ threshold, ≥ 75% of standard levels recover within the window "
        "with acceptable CV, a quantifiable range (LLOQ–ULOQ) exists, and controls fall inside "
        "their expected ranges. Limits are in-well concentrations."
    )
    st.dataframe(q, width="stretch", hide_index=True)
    st.dataframe(result.limits, width="stretch", hide_index=True)
with tabs[1]:
    st.caption(
        "Weighted fit through the standards (log x-axis). The shaded band is the validated "
        "quantification range. Weighted residuals should scatter evenly around zero."
    )
    keys = list(result.fits)
    key = st.selectbox("Plate", keys, format_func=lambda k: f"{k[0]} · {k[1]}")
    w = result.wells[(result.wells.plate_id == key[0]) & (result.wells.analyte == key[1])]
    lim = result.limits.set_index(["plate_id", "analyte"]).loc[key]
    left, right = st.columns(2)
    left.pyplot(
        P.plot_standard_curve(
            result.fits[key], w, lim, title=f"{key[0]} · {key[1]}", units=unit_map.get(key[1], "")
        )
    )
    right.pyplot(P.plot_plate_heatmap(w, title=f"{key[0]} raw OD"))
    ct = result.curve_table()
    st.dataframe(ct[(ct.plate_id == key[0]) & (ct.analyte == key[1])], hide_index=True)
    rec = result.recovery
    st.markdown("**Standard recovery** (back-calculated ÷ nominal)")
    st.dataframe(rec[(rec.plate_id == key[0]) & (rec.analyte == key[1])], hide_index=True)
with tabs[2]:
    st.caption(
        "below_lloq / above_uloq: outside the validated range (censored, not reported as a "
        "number) · below_lod: not distinguishable from blank · cv_high: replicates disagree."
    )
    only_flagged = st.checkbox("Flagged only", True)
    s = result.samples[result.samples["flags"] != ""] if only_flagged else result.samples
    st.dataframe(s, width="stretch", hide_index=True)
    st.markdown("**Wells excluded as outliers**")
    st.dataframe(result.wells[result.wells.outlier], width="stretch", hide_index=True)
with tabs[3]:
    st.caption(
        "Samples run at several dilutions should agree after correction (±20%). A concentration "
        "that rises with dilution suggests matrix interference or the high-dose hook effect."
    )
    st.dataframe(result.linearity, width="stretch", hide_index=True)
with tabs[4]:
    st.caption(
        "Bridge controls shared by every plate estimate a multiplicative plate factor; "
        "corrected concentration = raw ÷ factor. Lines should be flatter after correction."
    )
    st.dataframe(result.plate_factors, width="stretch", hide_index=True)
    st.dataframe(result.bridge_cv, width="stretch", hide_index=True)
    multi = result.results.groupby("analyte").plate_id.nunique()
    for a in multi[multi > 1].index:
        st.pyplot(P.plot_normalization(result.results, a, units=unit_map.get(a, "")))
with tabs[5]:
    if result.stats.pairwise.empty:
        st.info("This dataset has no treatment groups to compare (or too few quantifiable values).")
    st.caption(
        "Shapiro-Wilk picks Welch t / Welch ANOVA vs Mann-Whitney / Kruskal-Wallis. Omnibus "
        "p-values are Benjamini-Hochberg adjusted across analytes, pairwise p-values Holm "
        "adjusted. Hollow points are censored values that were substituted."
    )
    st.dataframe(result.stats.omnibus, width="stretch", hide_index=True)
    st.dataframe(result.stats.pairwise, width="stretch", hide_index=True)
    analytes = list(result.stats.data.analyte.unique()) if len(result.stats.data) else []
    if analytes:
        cols = st.columns(len(analytes))
        for col, a in zip(cols, analytes, strict=True):
            col.pyplot(
                P.plot_group_comparison(
                    result.stats.data,
                    a,
                    group_order=cfg.group_order,
                    pairwise=result.stats.pairwise,
                    units=unit_map.get(a, ""),
                    log=cfg.stats.transform == "log10",
                )
            )
        st.dataframe(
            result.stats.data[["plate_id", "analyte", "sample_id", "group", "value", "handling"]],
            hide_index=True,
        )
with tabs[6]:
    report_key = f"report::{source}::{json.dumps(cfg.as_dict(), default=str)}"
    if st.button("Build HTML report"):
        st.session_state[report_key] = build_report(
            result, title="ELISA analysis report", subtitle=source, units=unit_map
        )
    if report_key in st.session_state:
        st.download_button(
            "Download HTML report",
            st.session_state[report_key],
            "elisa_report.html",
            "text/html",
        )
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        out = tmpdir / "tables"
        for name, path in result.write_tables(out).items():
            zf.write(path, f"{name}.csv")
        zf.write(out / "config.json", "config.json")
    st.download_button(
        "Download all result tables (zip of CSVs)",
        buf.getvalue(),
        "elisa_results.zip",
        "application/zip",
    )
    st.download_button(
        "Download results.csv", result.results.to_csv(index=False), "results.csv", "text/csv"
    )
