"""Command-line interface: ``immunoassay simulate`` and ``immunoassay run``."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from immunoassay import qc, stats
from immunoassay.pipeline import PipelineConfig, run_pipeline
from immunoassay.report import write_report
from immunoassay.simulate import default_demo_study


def _simulate(args: argparse.Namespace) -> int:
    study = default_demo_study(seed=args.seed)
    manifest = study.write(args.out)
    print(f"wrote {len(study.plates)} plates to {args.out} (manifest: {manifest})")
    return 0


def _run(args: argparse.Namespace) -> int:
    manifest = Path(args.manifest)
    units = {}
    m = pd.read_csv(manifest, dtype=str)
    if "units" in m and "analyte" in m:
        units = dict(zip(m["analyte"], m["units"], strict=True))
    ranges = None
    ctrl_path = Path(args.controls) if args.controls else manifest.parent / "control_ranges.csv"
    if ctrl_path.exists():
        ranges = pd.read_csv(ctrl_path)
    cfg = PipelineConfig(
        model=args.model,
        weighting=args.weighting,
        qc=qc.QCConfig(max_cv_pct=args.max_cv),
        stats=stats.StatsConfig(transform=args.transform, censored=args.censored),
        normalize=not args.no_normalize,
        control_ranges=ranges,
        group_order=args.group_order.split(",") if args.group_order else None,
        units=units,
    )
    result = run_pipeline(manifest, cfg)
    out = Path(args.out)
    result.write_tables(out)
    report = write_report(result, out / "report.html", title=args.title, subtitle=str(manifest))
    print(result.plate_qc[["plate_id", "analyte", "passed", "reasons"]].to_string(index=False))
    print(f"\ntables + report written to {out} (open {report})")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="immunoassay", description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    s = sub.add_parser("simulate", help="write the two-cytokine demo study")
    s.add_argument("--out", default="data/simulated/demo")
    s.add_argument("--seed", type=int, default=17)
    s.set_defaults(func=_simulate)

    r = sub.add_parser("run", help="run the pipeline on a manifest CSV")
    r.add_argument("manifest")
    r.add_argument("--out", default="reports/latest")
    r.add_argument("--controls", help="control ranges CSV (analyte,sample_id,low,high)")
    r.add_argument("--model", default="auto", choices=["auto", "4pl", "5pl"])
    r.add_argument("--weighting", default="1/y^2", choices=["1/y^2", "1/y", "none"])
    r.add_argument("--max-cv", type=float, default=20.0)
    r.add_argument("--transform", default="log10", choices=["log10", "none"])
    r.add_argument("--censored", default="substitute", choices=["substitute", "exclude"])
    r.add_argument("--no-normalize", action="store_true")
    r.add_argument("--group-order", help="comma-separated; first group is the reference")
    r.add_argument("--title", default="ELISA analysis report")
    r.set_defaults(func=_run)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
