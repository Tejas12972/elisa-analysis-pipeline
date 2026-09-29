"""End-to-end: simulate plates with known truth, run everything, check recovery."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from immunoassay.cli import main as cli_main
from immunoassay.pipeline import PipelineConfig, run_pipeline
from immunoassay.report import build_report
from immunoassay.simulate import AnalyteSpec, NoiseSpec, default_demo_study, simulate_study

GROUPS = ["vehicle", "LPS", "LPS+drug"]


@pytest.fixture(scope="module")
def demo(tmp_path_factory):
    study = default_demo_study()
    manifest = study.write(tmp_path_factory.mktemp("demo"))
    cfg = PipelineConfig(
        control_ranges=study.control_ranges(), group_order=GROUPS, units=study.units
    )
    return study, run_pipeline(manifest, cfg)


def _accuracy(result, truth, col):
    res = result.results[(result.results["type"] == "sample") & result.results[col].notna()]
    m = res.merge(truth, on=["analyte", "sample_id"])
    return m.assign(ratio=m[col] / m["true_conc"])


def test_all_plates_pass(demo):
    _, r = demo
    assert r.plate_qc["passed"].all(), r.plate_qc[["plate_id", "reasons"]]


def test_recovered_concentrations_within_tolerance(demo):
    """Normalized results equal truth × (geometric-mean calibrator bias)⁻¹ within noise.

    Bridge normalization aligns every plate to the *average* plate, so any bias
    shared by all plates remains. That shared part is removed here before
    judging accuracy.
    """
    study, r = demo
    acc = _accuracy(r, study.truth, "conc_norm")
    shared = study.plate_factors.groupby("analyte")["calibrator_bias"].apply(
        lambda k: np.exp(np.log(k).mean())
    )
    acc["ratio_adj"] = acc["ratio"] * acc["analyte"].map(shared)
    err = np.abs(np.log(acc["ratio_adj"]))
    assert len(acc) >= 60
    assert np.exp(err.median()) - 1 < 0.08
    assert (np.exp(err) - 1 < 0.25).mean() > 0.95


def test_normalization_reduces_bias_and_bridge_cv(demo):
    study, r = demo
    raw = _accuracy(r, study.truth, "conc")
    norm = _accuracy(r, study.truth, "conc_norm")
    assert np.log(norm["ratio"]).std() < np.log(raw["ratio"]).std()
    cv = r.bridge_cv.dropna()
    assert (cv["inter_plate_cv_pct_after"] < cv["inter_plate_cv_pct_before"]).mean() >= 0.75


def test_plate_factors_track_calibrator_bias(demo):
    study, r = demo
    f = r.plate_factors.merge(study.plate_factors, on=["plate_id", "analyte"])
    for _, g in f.groupby("analyte"):
        # factor ∝ 1 / k (samples read at true / k)
        rel_est = g["factor"] / g["factor"].iloc[0]
        rel_true = g["calibrator_bias"].iloc[0] / g["calibrator_bias"]
        np.testing.assert_allclose(rel_est, rel_true, rtol=0.12)


def test_injected_standard_outliers_detected(demo):
    study, r = demo
    w = r.wells.set_index(["plate_id", "well"])
    injected = [tuple(x) for x in study.injected_outliers.to_numpy()]
    std_injected = [x for x in injected if w.loc[x, "type"] == "standard"]
    for key in std_injected:
        assert bool(w.loc[key, "outlier"]), key


def test_biology_detected(demo):
    _, r = demo
    om = r.stats.omnibus.set_index("analyte")
    assert bool(om.loc["IL-6", "significant"])
    pw = r.stats.pairwise.set_index(["analyte", "reference", "group"])
    il6_lps = pw.loc[("IL-6", "vehicle", "LPS")]
    assert il6_lps["significant"] and 3 < il6_lps["estimate"] < 8  # true fold change = 5
    assert il6_lps["ci_low"] < 5 < il6_lps["ci_high"]


def test_no_calibrator_bias_means_raw_is_accurate(tmp_path):
    noise = NoiseSpec(calibrator_bias_cv=0.0, od_scale_cv=0.15)
    study = simulate_study(noise=noise, seed=11)
    r = run_pipeline(study.write(tmp_path), PipelineConfig(control_ranges=study.control_ranges()))
    acc = _accuracy(r, study.truth, "conc")
    assert np.exp(np.abs(np.log(acc["ratio"])).median()) - 1 < 0.06


def test_hook_effect_flagged_by_linearity(tmp_path):
    spec = AnalyteSpec(
        group_geo_means={"low": 30.0, "high": 4000.0}, geo_sd=1.2, controls={"CTRL": 100.0}
    )
    noise = NoiseSpec(hook_onset=600.0, od_scale_cv=0.0)
    study = simulate_study(
        [spec],
        n_per_group=8,
        n_linearity_per_plate=16,
        linearity_dilution=16.0,
        noise=noise,
        seed=4,
    )
    r = run_pipeline(study.write(tmp_path))
    lin = r.linearity.merge(study.truth, on=["analyte", "sample_id", "group"])
    high = lin[lin["group"] == "high"]
    assert len(high) > 0
    # Hooked high samples either rise with dilution or only quantify when diluted.
    assert (
        high["trend"].eq("rises_with_dilution")
        | high["interpretation"].str.contains("higher dilution")
    ).mean() > 0.5


def test_report_and_tables(demo, tmp_path):
    _, r = demo
    html = build_report(r, title="t", units={"IL-6": "pg/mL"})
    assert html.startswith("<!doctype html>")
    for section in (
        "Plate acceptance",
        "Standard curves",
        "Dilution linearity",
        "Inter-plate normalization",
        "Group comparisons",
    ):
        assert section in html
    assert html.count("data:image/png;base64") >= 4 + 4 + 2 + 2
    paths = r.write_tables(tmp_path)
    assert pd.read_csv(paths["results"]).shape[0] == len(r.results)
    assert (tmp_path / "config.json").exists()


def test_cli_simulate_and_run(tmp_path, capsys):
    assert cli_main(["simulate", "--out", str(tmp_path / "sim")]) == 0
    assert (
        cli_main(
            [
                "run",
                str(tmp_path / "sim" / "manifest.csv"),
                "--out",
                str(tmp_path / "out"),
                "--group-order",
                ",".join(GROUPS),
            ]
        )
        == 0
    )
    assert (tmp_path / "out" / "report.html").exists()
    assert "True" in capsys.readouterr().out
