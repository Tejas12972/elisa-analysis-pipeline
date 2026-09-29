# ELISA Analysis Pipeline

Python pipeline that turns raw 96-well sandwich-ELISA plate-reader exports into QC'd
concentrations, group statistics, and an HTML report. Portfolio project for biotech
data-analysis co-op applications, so scientific correctness, tests, and explanatory docstrings
matter as much as working code. Public repo: github.com/Tejas12972/elisa-analysis-pipeline
(live site + in-browser app on GitHub Pages).

## Status
Complete. The package (import name `immunoassay`, distribution `elisa-analysis-pipeline`),
two real-data notebooks, CI, and a GitHub Pages site built by `scripts/build_site.py`:
the landing page (`docs/`), freshly generated reports, and the Streamlit app running in-browser
via stlite (`scripts/site_app_template.html`). `.github/workflows/pages.yml` gates each deploy
on `scripts/check_site.py`, which boots the WebAssembly app in headless Chromium.

## Run / test
```bash
python3.11 -m venv .venv && .venv/bin/pip install -e ".[dev,app,notebooks]"  # plain pip, NOT uv (x86-wheel issue on this Mac)
.venv/bin/pytest --cov=immunoassay                  # 135 tests
.venv/bin/ruff check . && .venv/bin/ruff format --check .
.venv/bin/immunoassay simulate --out data/simulated/demo
.venv/bin/immunoassay run data/simulated/demo/manifest.csv --out reports/demo
.venv/bin/streamlit run app/streamlit_app.py
.venv/bin/python scripts/prepare_sleep_fragmentation.py   # data/raw -> data/processed
.venv/bin/python scripts/prepare_elisatools.py
.venv/bin/python scripts/build_notebooks.py              # regenerates + executes notebooks/
.venv/bin/python scripts/build_site.py                   # -> _site/ (Pages)
python scripts/check_site.py --channel chrome            # needs `pip install playwright`
```

## Conventions / constraints
- src layout; type hints; ruff (line length 100) is the formatter. It also lints the notebooks,
  so edit notebook *source* in `scripts/build_notebooks.py`, never the .ipynb directly.
- The tidy frame is one row per well (`od`, `od_status` ok|missing|overflow, `od_net`, `blank_*`).
  Replicate key = `(plate_id, analyte, type, sample_id, dilution)`.
- Input validation collects *all* problems into one `InputFileError` subclass. Legitimate
  oddities raise a `PlateWarning` instead.
- Stage order in `pipeline.analyze_wells` is deliberate (see its docstring). Controls are judged
  BEFORE normalization.
- Scientific choices are documented in module docstrings AND the README table. Keep them in
  sync, and re-verify notebook prose against executed outputs whenever numbers change.
- Curve fits include the zero standard (blank wells at x = 0). Weights are 1/y² on the *gross* OD.
- Stats default: log10, Welch/MW auto, BH across analytes, LLOQ/√2 substitution with a
  study-wide common LLOQ. `censored="extrapolate"` is for labelled sensitivity analyses only.
- `data/raw/` is byte-for-byte upstream (see each PROVENANCE.md). Never edit it; derive into
  `data/processed/` via `scripts/`.
- Local pandas is 3.x, but the browser app runs Pyodide (Python 3.13, pandas 2.3, SciPy 1.14),
  pinned in `requirements/browser-stack.txt` and tested in CI. Code must work on both.
  `DataFrame.flags` is a pandas attribute, so always use `df["flags"]`.
- git here is 2.15: no `git branch --show-current`, no `git init -b`.
