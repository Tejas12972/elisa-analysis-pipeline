"""Tests for QC: outliers, limits, replicate flags, controls, plate acceptance."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from immunoassay import qc
from immunoassay.curves import back_calculate, fit_standard_curves, four_pl, standard_recovery

P = (0.02, 1.1, 150.0, 3.0)
LEVELS = 1000 / 2 ** np.arange(7)


def test_modified_z_flags_gross_error_and_handles_zero_spread():
    z = qc.modified_z(np.array([1.00, 1.01, 1.02, 1.60]))
    assert abs(z[-1]) > 3.5 and np.all(np.abs(z[:3]) < 3.5)
    assert np.all(qc.modified_z(np.array([1.0, 1.0, 1.0])) == 0)


def _wells(ods: dict[str, list[float]], typ="sample") -> pd.DataFrame:
    rows = []
    for sid, vals in ods.items():
        for i, v in enumerate(vals):
            rows.append(
                {
                    "plate_id": "P",
                    "analyte": "X",
                    "type": typ,
                    "sample_id": sid,
                    "dilution": 1.0,
                    "od": v,
                    "od_status": "ok",
                    "well": f"{sid}{i}",
                }
            )
    return pd.DataFrame(rows)


def test_replicate_outlier_requires_failing_cv():
    w = _wells({"bad": [1.0, 1.02, 1.9], "fine": [1.0, 1.01, 1.08], "dup": [1.0, 2.0]})
    out = qc.flag_replicate_outliers(w)
    flagged = out.loc[out.outlier, "well"].tolist()
    assert flagged == ["bad2"]  # 'fine' passes CV, 'dup' has n = 2
    assert "modified z" in out.loc[out.outlier, "outlier_reason"].iloc[0]


def test_replicate_outlier_never_removes_majority():
    w = _wells({"s": [1.0, 3.0, 5.0]})
    out = qc.flag_replicate_outliers(w)
    assert out["outlier"].sum() <= 1


def _plate(noise_seed=None, spike: dict | None = None) -> pd.DataFrame:
    rng = np.random.default_rng(noise_seed)
    rows = []
    for i, c in enumerate(LEVELS):
        for rep in range(2):
            y = float(four_pl(c, *P))
            if noise_seed is not None:
                y *= 1 + rng.normal(0, 0.02)
            rows.append(
                {
                    "plate_id": "P",
                    "analyte": "X",
                    "type": "standard",
                    "sample_id": f"STD{i}",
                    "concentration": c,
                    "dilution": 1.0,
                    "od_net": y,
                    "od": y + 0.05,
                    "od_status": "ok",
                    "well": f"{'ABCDEFG'[i]}{rep + 1}",
                    "blank_mean": 0.05,
                    "blank_sd": 0.004,
                    "group": None,
                }
            )
    df = pd.DataFrame(rows)
    for well, factor in (spike or {}).items():
        df.loc[df.well == well, "od_net"] *= factor
    df["outlier"] = False
    df["outlier_reason"] = ""
    return df


def test_standard_residual_outlier_found_and_curve_improves():
    wells = _plate(noise_seed=1, spike={"D1": 0.5})
    fits = fit_standard_curves(wells, model="4pl")
    flagged = qc.flag_standard_residual_outliers(wells, fits)
    assert flagged.loc[flagged.outlier, "well"].tolist() == ["D1"]
    refit = fit_standard_curves(flagged, model="4pl")[("P", "X")]
    assert refit.params["c"] == pytest.approx(150, rel=0.05)


def test_two_standard_outliers_both_found():
    wells = _plate(noise_seed=2, spike={"B1": 1.4, "E2": 0.55})
    fits = fit_standard_curves(wells, model="4pl")
    flagged = qc.flag_standard_residual_outliers(wells, fits)
    assert set(flagged.loc[flagged.outlier, "well"]) == {"B1", "E2"}


def test_clean_standards_not_flagged():
    wells = _plate(noise_seed=3)
    fits = fit_standard_curves(wells, model="4pl")
    assert not qc.flag_standard_residual_outliers(wells, fits)["outlier"].any()


def _limits_for(wells):
    fits = fit_standard_curves(wells, model="4pl")
    wells = back_calculate(wells, fits)
    rec = standard_recovery(wells)
    return wells, fits, rec, qc.plate_limits(wells, fits, rec)


def test_limits_on_clean_plate():
    _, fits, _, lim = _limits_for(_plate(noise_seed=4))
    row = lim.iloc[0]
    assert row["lloq"] == pytest.approx(LEVELS.min())
    assert row["uloq"] == pytest.approx(LEVELS.max())
    assert row["n_std_levels_pass"] == 7
    # True floor a = 0.02 > 3·SD = 0.012, so the curve-floor rule applies.
    fit = fits[("P", "X")]
    assert row["lod_basis"] == "curve floor + k·SD"
    assert row["lod"] == pytest.approx(fit.back_calculate(fit.params["a"] + 0.012)[0][0])
    assert 0 < row["lod"] < row["lloq"]


def test_lod_uses_blank_rule_when_above_floor():
    wells = _plate(noise_seed=4).assign(blank_sd=0.05)
    _, fits, _, lim = _limits_for(wells)
    assert lim.iloc[0]["lod_basis"] == "mean blank + k·SD"
    assert lim.iloc[0]["lod"] == pytest.approx(fits[("P", "X")].back_calculate(0.15)[0][0])


def test_lloq_uses_longest_contiguous_passing_run():
    wells, fits, rec, _ = _limits_for(_plate(noise_seed=4))
    rec = rec.copy()
    lowest_two = rec["concentration"].nsmallest(2).index
    rec.loc[lowest_two[1], "recovery_pct"] = 150  # second-lowest level fails
    lim = qc.plate_limits(wells, fits, rec).iloc[0]
    assert lim["lloq"] == pytest.approx(sorted(LEVELS)[2])
    assert "15.62" in lim["failing_levels"] or "31.25" in lim["failing_levels"]


def _with_samples(conc_by_sample: dict[str, tuple[float, float]]):
    wells = _plate(noise_seed=5)
    rows = []
    for sid, (c1, c2) in conc_by_sample.items():
        for rep, c in enumerate((c1, c2)):
            y = float(four_pl(c, *P))
            rows.append(
                {
                    "plate_id": "P",
                    "analyte": "X",
                    "type": "sample",
                    "sample_id": sid,
                    "concentration": np.nan,
                    "dilution": 2.0,
                    "od_net": y,
                    "od": y + 0.05,
                    "od_status": "ok",
                    "well": f"{sid}_{rep}",
                    "blank_mean": 0.05,
                    "blank_sd": 0.004,
                    "group": "g",
                    "outlier": False,
                    "outlier_reason": "",
                }
            )
    wells = pd.concat([wells, pd.DataFrame(rows)], ignore_index=True)
    wells, _, _, lim = _limits_for(wells)
    return qc.summarize_replicates(wells, lim).set_index("sample_id"), lim


def test_replicate_summary_flags():
    s, _ = _with_samples(
        {
            "good": (100, 104),
            "noisy": (60, 140),
            "low": (3, 3.2),
            "high": (5000, 5200),
            "floor": (0, 0),
        }
    )
    assert s.loc["good", "flags"] == ""
    assert s.loc["good", "conc_well"] == pytest.approx(102, rel=0.03)
    assert s.loc["noisy", "cv_high"]
    assert s.loc["low", "range_flag"] == "below_lloq"
    assert s.loc["high", "range_flag"] == "above_uloq"
    assert s.loc["floor", "below_lod"]
    assert s.loc["good", "dilution"] == 2.0


def test_controls_and_plate_summary():
    s, lim = _with_samples({"good": (100, 104)})
    samples = s.reset_index()
    samples.loc[0, "type"] = "control"
    samples["conc_neat"] = samples["conc_well"] * samples["dilution"]
    ranges = pd.DataFrame({"analyte": ["X"], "sample_id": ["good"], "low": [150], "high": [250]})
    ctrl = qc.check_controls(samples, ranges)
    assert bool(ctrl["passed"].iloc[0])
    wells = _plate(noise_seed=5)
    fits = fit_standard_curves(wells, model="4pl")
    summary = qc.plate_summary(fits, lim, ctrl, samples)
    assert bool(summary["passed"].iloc[0])

    bad_ranges = ranges.assign(low=500, high=600)
    ctrl_bad = qc.check_controls(samples, bad_ranges)
    summary_bad = qc.plate_summary(fits, lim, ctrl_bad, samples)
    assert not bool(summary_bad["passed"].iloc[0])
    assert "controls out of range" in summary_bad["reasons"].iloc[0]


def test_controls_without_ranges_are_not_judged():
    s, _ = _with_samples({"c": (100, 104)})
    samples = s.reset_index().assign(type="control")
    samples["conc_neat"] = samples["conc_well"] * samples["dilution"]
    ctrl = qc.check_controls(samples, None)
    assert ctrl["passed"].isna().all()


def test_edge_tolerance_applies_only_at_the_limits():
    wells, fits, rec, _ = _limits_for(_plate(noise_seed=4))
    rec = rec.sort_values("concentration").reset_index(drop=True)
    edge_case = rec.copy()
    edge_case.loc[0, "recovery_pct"] = 123  # lowest level: within ±25% → still the LLOQ
    lim = qc.plate_limits(wells, fits, edge_case).iloc[0]
    assert lim["lloq"] == pytest.approx(rec.loc[0, "concentration"])
    assert lim["n_std_levels_pass"] == 7
    interior = rec.copy()
    interior.loc[2, "recovery_pct"] = 123  # interior level: must meet ±20% → breaks the run
    lim2 = qc.plate_limits(wells, fits, interior).iloc[0]
    assert lim2["lloq"] == pytest.approx(rec.loc[2, "concentration"])  # 3rd level becomes an edge
    assert lim2["n_std_levels_pass"] == 6


def test_quantification_range_helper():
    t, f = True, False
    arr = np.array
    assert qc._quantification_range(arr([t, t, t]), arr([t, t, t])) == (0, 2)
    assert qc._quantification_range(arr([f, t, f, t]), arr([t, t, t, t])) == (0, 2)
    assert qc._quantification_range(arr([f, f]), arr([f, f])) == (None, None)
