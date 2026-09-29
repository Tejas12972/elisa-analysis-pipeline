"""Tests for 4PL/5PL fitting and back-calculation (immunoassay.curves)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from immunoassay.curves import (
    back_calculate,
    compute_weights,
    fit_curve,
    fit_standard_curves,
    five_pl,
    four_pl,
    inverse_five_pl,
    inverse_four_pl,
    standard_recovery,
)

LEVELS = 1000 / 2 ** np.arange(8)
TRUE = {"a": 0.02, "b": 1.1, "c": 150.0, "d": 3.0}


def noisy(y: np.ndarray, cv: float, seed: int = 0) -> np.ndarray:
    return y * (1 + np.random.default_rng(seed).normal(0, cv, y.shape))


def test_four_pl_limits_and_midpoint():
    a, b, c, d = TRUE.values()
    assert four_pl(0.0, a, b, c, d) == pytest.approx(a)
    assert four_pl(1e12, a, b, c, d) == pytest.approx(d, rel=1e-6)
    assert four_pl(c, a, b, c, d) == pytest.approx((a + d) / 2)


def test_five_pl_reduces_to_four_pl():
    x = np.geomspace(1, 1e4, 20)
    np.testing.assert_allclose(five_pl(x, *TRUE.values(), 1.0), four_pl(x, *TRUE.values()))


@pytest.mark.parametrize("g", [1.0, 0.6, 1.8])
def test_inverse_round_trips(g):
    x = np.geomspace(1, 5000, 25)
    y = five_pl(x, *TRUE.values(), g)
    np.testing.assert_allclose(inverse_five_pl(y, *TRUE.values(), g), x, rtol=1e-8)


def test_inverse_is_nan_outside_asymptotes():
    out = inverse_four_pl(np.array([0.01, 0.02, 3.0, 3.5]), *TRUE.values())
    assert np.isnan(out).all()


def test_weights_normalised_and_floored():
    w = compute_weights(np.array([0.0, 0.5, 1.0, 2.0]), "1/y^2")
    assert w.mean() == pytest.approx(1.0)
    assert np.isfinite(w).all()
    assert w[0] == w.max()
    np.testing.assert_allclose(compute_weights(np.array([1.0, 2.0]), "none"), [1.0, 1.0])
    with pytest.raises(ValueError, match="unknown weighting"):
        compute_weights(np.array([1.0]), "1/x")


def test_weight_offset_uses_gross_signal():
    y = np.array([0.02, 1.0])
    ratio_net = compute_weights(y, "1/y^2")[0] / compute_weights(y, "1/y^2")[1]
    ratio_gross = compute_weights(y, "1/y^2", 0.08)[0] / compute_weights(y, "1/y^2", 0.08)[1]
    assert ratio_gross < ratio_net


def test_fit_recovers_known_parameters():
    x = np.repeat(LEVELS, 2)
    y = noisy(four_pl(x, *TRUE.values()), 0.02)
    fit = fit_curve(x, y, model="4pl")
    assert fit.model == "4pl"
    assert fit.params["c"] == pytest.approx(TRUE["c"], rel=0.1)
    assert fit.params["b"] == pytest.approx(TRUE["b"], rel=0.1)
    assert fit.r_squared > 0.99
    assert all(np.isfinite(v) for v in fit.stderr.values())
    assert fit.increasing


def test_auto_prefers_4pl_for_symmetric_data():
    x = np.repeat(LEVELS, 2)
    y = noisy(four_pl(x, *TRUE.values()), 0.02, seed=3)
    fit = fit_curve(x, y, model="auto")
    assert fit.model == "4pl"
    assert set(fit.comparison) == {"4pl", "5pl"}


def test_auto_picks_5pl_for_strongly_asymmetric_data():
    x = np.repeat(np.geomspace(2, 20000, 10), 3)
    y = noisy(five_pl(x, 0.02, 1.8, 400.0, 3.0, 0.25), 0.01, seed=1)
    fit = fit_curve(x, y, model="auto")
    assert fit.model == "5pl"
    assert fit.comparison["5pl"] < fit.comparison["4pl"] - 2


def test_decreasing_competitive_curve():
    x = np.repeat(LEVELS, 2)
    y = noisy(four_pl(x, 2.5, 1.0, 100.0, 0.1), 0.02)
    fit = fit_curve(x, y, model="4pl")
    assert not fit.increasing
    conc, flag = fit.back_calculate(four_pl(np.array([50.0]), 2.5, 1.0, 100.0, 0.1))
    assert conc[0] == pytest.approx(50, rel=0.1)
    assert flag[0] == "ok"


def test_weighting_improves_low_end_accuracy():
    """Constant-CV noise: 1/y² weighting should back-calculate the lowest standards better."""
    errs = {"1/y^2": [], "none": []}
    x = np.repeat(LEVELS, 2)
    for seed in range(30):
        y = noisy(four_pl(x, *TRUE.values()), 0.08, seed)
        for w in errs:
            fit = fit_curve(x, y, model="4pl", weighting=w)
            est, _ = fit.back_calculate(four_pl(LEVELS[-2:], *TRUE.values()))
            errs[w].append(np.nanmean(np.abs(np.log(est / LEVELS[-2:]))))
    assert np.mean(errs["1/y^2"]) < np.mean(errs["none"])


def test_back_calculate_flags():
    x = np.repeat(LEVELS, 2)
    fit = fit_curve(x, four_pl(x, *TRUE.values()), model="4pl")
    conc, flag = fit.back_calculate(np.array([-0.5, np.nan, 5.0, 1.0]))
    assert list(flag) == ["below_curve", "no_signal", "above_curve", "ok"]
    assert np.isnan(conc[:3]).all() and np.isfinite(conc[3])


def test_fit_requires_enough_points():
    with pytest.raises(ValueError, match="need >= 5"):
        fit_curve(np.array([1, 2, 3, 4.0]), np.array([0.1, 0.2, 0.3, 0.4]))
    with pytest.raises(ValueError, match="4 non-zero levels"):
        fit_curve(np.array([1, 1, 2, 2, 3, 3.0]), np.array([0.1, 0.1, 0.2, 0.2, 0.3, 0.3]))


def _std_frame(plate="P1", analyte="X", noise=0.0):
    x = np.repeat(LEVELS[:7], 2)
    y = noisy(four_pl(x, *TRUE.values()), noise) if noise else four_pl(x, *TRUE.values())
    return pd.DataFrame(
        {
            "plate_id": plate,
            "analyte": analyte,
            "type": "standard",
            "concentration": x,
            "od_net": y,
            "od_status": "ok",
            "blank_mean": 0.05,
            "outlier": False,
            "well": [f"W{i}" for i in range(len(x))],
        }
    )


def test_plate_level_fit_backcalc_and_recovery():
    std = _std_frame()
    samples = pd.DataFrame(
        {
            "plate_id": "P1",
            "analyte": "X",
            "type": "sample",
            "concentration": np.nan,
            "od_net": [four_pl(100.0, *TRUE.values()), np.nan],
            "od_status": ["ok", "overflow"],
            "blank_mean": 0.05,
            "outlier": False,
            "well": ["S1", "S2"],
        }
    )
    wells = pd.concat([std, samples], ignore_index=True)
    fits = fit_standard_curves(wells, model="4pl")
    assert fits[("P1", "X")].weight_offset == 0.05
    out = back_calculate(wells, fits)
    s = out[out.type == "sample"].set_index("well")
    assert s.loc["S1", "conc_well"] == pytest.approx(100, rel=1e-3)
    assert s.loc["S2", "curve_flag"] == "above_curve"
    rec = standard_recovery(out)
    assert len(rec) == 7
    np.testing.assert_allclose(rec["recovery_pct"], 100, rtol=1e-3)


def test_outlier_standards_are_not_fitted():
    std = _std_frame()
    std.loc[0, "od_net"] = 0.01
    std.loc[0, "outlier"] = True
    fit = fit_standard_curves(std, model="4pl")[("P1", "X")]
    assert fit.n == 13


def test_zero_standard_anchors_lower_asymptote():
    """A biased floor extrapolation is pulled back to the measured zero standard."""
    x = np.repeat(LEVELS[:6], 2)  # lowest standard ~31 pg/mL: the floor is extrapolated
    y = noisy(four_pl(x, *TRUE.values()), 0.03, seed=7)
    blank = pd.DataFrame(
        {
            "plate_id": "P",
            "analyte": "X",
            "type": "blank",
            "concentration": np.nan,
            "od_net": [0.021, 0.019],
            "od_status": "ok",
            "blank_mean": 0.05,
            "outlier": False,
            "well": ["H1", "H2"],
        }
    )
    std = pd.DataFrame(
        {
            "plate_id": "P",
            "analyte": "X",
            "type": "standard",
            "concentration": x,
            "od_net": y,
            "od_status": "ok",
            "blank_mean": 0.05,
            "outlier": False,
            "well": [f"W{i}" for i in range(len(x))],
        }
    )
    wells = pd.concat([std, blank], ignore_index=True)
    with_zero = fit_standard_curves(wells, model="4pl")[("P", "X")]
    without = fit_standard_curves(wells, model="4pl", include_zero=False)[("P", "X")]
    assert with_zero.n == without.n + 2
    assert abs(with_zero.params["a"] - 0.02) < abs(without.params["a"] - 0.02) + 1e-9
    assert with_zero.params["a"] == pytest.approx(0.02, abs=0.01)
