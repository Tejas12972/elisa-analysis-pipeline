"""Build the static website deployed to GitHub Pages (``_site/``).

The site has three parts:

* ``index.html`` and ``img/``: the landing page, copied from ``docs/``.
* ``reports/``: example HTML reports, regenerated from the current code on every build so
  they can never drift from the pipeline.
* ``app/``: the Streamlit app running entirely in the visitor's browser via
  `stlite <https://github.com/whitphx/stlite>`_ (Streamlit on Pyodide/WebAssembly).
  ``bundle.zip`` holds the app script, the ``immunoassay`` package, and the processed real
  datasets. stlite unpacks it into the in-browser file system, and Streamlit puts the
  script's directory on ``sys.path``, so ``import immunoassay`` just works.

Usage::

    python scripts/build_site.py            # -> _site/
    python -m http.server -d _site 8000     # preview at http://localhost:8000
"""

from __future__ import annotations

import argparse
import json
import shutil
import warnings
import zipfile
from pathlib import Path

import pandas as pd

from immunoassay import stats as st
from immunoassay.pipeline import PipelineConfig, run_pipeline
from immunoassay.report import write_report

REPO = Path(__file__).resolve().parents[1]
STLITE_VERSION = "1.9.2"
# Packages micropip installs from the Pyodide distribution; Streamlit ships with stlite.
REQUIREMENTS = ["numpy", "pandas", "scipy", "statsmodels", "matplotlib", "openpyxl"]

# The app page is a str.format template ({version}, {requirements}; CSS braces doubled).
APP_TEMPLATE = Path(__file__).with_name("site_app_template.html")


def build_bundle(dest: Path) -> list[str]:
    """Zip the app, the package source and the processed datasets for stlite."""
    members: list[tuple[Path, str]] = [(REPO / "app/streamlit_app.py", "streamlit_app.py")]
    pkg = REPO / "src/immunoassay"
    members += [(p, f"immunoassay/{p.relative_to(pkg)}") for p in sorted(pkg.rglob("*.py"))]
    for dataset in ("sleep_fragmentation_il6", "elisatools_feng2019"):
        root = REPO / "data/processed" / dataset
        members += [
            (p, f"data/processed/{dataset}/{p.relative_to(root)}")
            for p in sorted(root.rglob("*"))
            if p.is_file()
        ]
    with zipfile.ZipFile(dest, "w", zipfile.ZIP_DEFLATED) as zf:
        for src, arc in members:
            zf.write(src, arc)
    return [arc for _, arc in members]


def build_reports(dest: Path) -> list[Path]:
    """Regenerate the three example reports from the current code."""
    dest.mkdir(parents=True, exist_ok=True)
    demo = REPO / "data/simulated/demo"
    il6 = REPO / "data/processed/sleep_fragmentation_il6"
    softmax = REPO / "data/processed/elisatools_feng2019"
    jobs = [
        (
            demo / "manifest.csv",
            PipelineConfig(
                control_ranges=pd.read_csv(demo / "control_ranges.csv"),
                group_order=["vehicle", "LPS", "LPS+drug"],
                units={"IL-6": "pg/mL", "TNF-α": "pg/mL"},
            ),
            "demo_report.html",
            "Demo study: LPS challenge, IL-6 and TNF-α",
            "Simulated data with known ground truth (immunoassay.simulate.default_demo_study)",
        ),
        (
            il6 / "manifest.csv",
            PipelineConfig(
                normalize=False,
                units={"IL-6": "pg/mL"},
                stats=st.StatsConfig(value_col="conc"),
            ),
            "sleep_fragmentation_il6_report.html",
            "Mouse serum IL-6 after acute sleep fragmentation (real data)",
            "Nguyen, Fields & Ashley, Dryad doi:10.5061/dryad.tdz08kq3p (CC0). Groups = "
            "treatment × time; the time-matched analysis is in notebooks/01.",
        ),
        (
            softmax / "manifest.csv",
            PipelineConfig(units={"analyte (undisclosed)": "pg/mL"}),
            "elisatools_softmax_report.html",
            "SoftMax Pro plates across two kit lots (real data)",
            "ELISAtools 0.1.8 example data, Feng Feng, MIT license; 5 plates, 2 kit lots",
        ),
    ]
    out = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for manifest, cfg, name, title, subtitle in jobs:
            out.append(
                write_report(
                    run_pipeline(manifest, cfg), dest / name, title=title, subtitle=subtitle
                )
            )
    return out


def main(out: Path) -> None:
    if out.exists():
        shutil.rmtree(out)
    shutil.copytree(REPO / "docs", out)
    (out / ".nojekyll").touch()
    app = out / "app"
    app.mkdir()
    files = build_bundle(app / "bundle.zip")
    (app / "index.html").write_text(
        APP_TEMPLATE.read_text(encoding="utf-8").format(
            version=STLITE_VERSION, requirements=json.dumps(REQUIREMENTS)
        ),
        encoding="utf-8",
    )
    reports = build_reports(out / "reports")
    size_kb = (app / "bundle.zip").stat().st_size / 1024
    print(f"site -> {out.relative_to(REPO)}")
    print(f"  app bundle: {len(files)} files, {size_kb:.0f} KB (stlite {STLITE_VERSION})")
    print(f"  reports: {', '.join(p.name for p in reports)}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=REPO / "_site")
    main(ap.parse_args().out)
