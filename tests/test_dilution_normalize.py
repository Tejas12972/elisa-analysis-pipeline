"""Tests for dilution handling and inter-plate normalization."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from immunoassay import dilution, normalize


def _rep(sid, dil, conc_well, flag="", cv_high=False, plate="P1", typ="sample", group="g"):
    return {
        "plate_id": plate,
        "analyte": "X",
        "type": typ,
        "sample_id": sid,
        "group": group,
        "dilution": float(dil),
        "conc_well": conc_well,
        "lloq_well": 10.0,
        "uloq_well": 1000.0,
        "range_flag": flag,
        "cv_high": cv_high,
        "flags": ";".join(f for f in (flag, "cv_high" if cv_high else "") if f),
    }


def test_apply_dilution_scales_values_and_limits():
    s = dilution.apply_dilution(pd.DataFrame([_rep("A", 4, 50.0)]))
    assert s.loc[0, "conc_neat"] == 200
    assert s.loc[0, "lloq_neat"] == 40
    assert s.loc[0, "uloq_neat"] == 4000


def test_linearity_flat_rising_and_single():
    s = dilution.apply_dilution(
        pd.DataFrame(
            [
                _rep("lin", 2, 100.0),
                _rep("lin", 8, 26.0),  # 200 vs 208
                _rep("hook", 2, 50.0),
                _rep("hook", 8, 40.0),  # 100 vs 320: rises with dilution
                _rep("one", 2, 800.0, "above_uloq"),
                _rep("one", 8, 150.0),
                _rep("solo", 2, 100.0),
            ]
        )
    )
    lin = dilution.dilution_linearity(s).set_index("sample_id")
    assert "solo" not in lin.index
    assert bool(lin.loc["lin", "linear"]) and lin.loc["lin", "trend"] == "flat"
    assert lin.loc["hook", "trend"] == "rises_with_dilution"
    assert not bool(lin.loc["hook", "linear"])
    assert "hook" in lin.loc["hook", "interpretation"]
    assert np.isnan(lin.loc["one", "linear"])
    assert "higher dilution" in lin.loc["one", "interpretation"]


def test_reportable_results_rules():
    s = dilution.apply_dilution(
        pd.DataFrame(
            [
                _rep("lin", 2, 100.0),
                _rep("lin", 8, 26.0),
                _rep("hook", 2, 50.0),
                _rep("hook", 8, 40.0),
                _rep("low", 2, np.nan, "below_lloq"),
                _rep("low", 8, np.nan, "below_lloq"),
                _rep("high", 2, np.nan, "above_uloq"),
                _rep("noisy", 2, 100.0, cv_high=True),
                _rep("mixed", 2, np.nan, "above_uloq"),
                _rep("mixed", 8, np.nan, "below_lloq"),
            ]
        )
    )
    res = dilution.reportable_results(s, dilution.dilution_linearity(s)).set_index("sample_id")
    assert res.loc["lin", "conc"] == pytest.approx(204)
    assert res.loc["hook", "conc"] == pytest.approx(320)  # most dilute when non-linear
    assert "most dilute" in res.loc["hook", "basis"]
    assert res.loc["low", "censor"] == "<LLOQ" and res.loc["low", "censor_limit"] == 20
    assert res.loc["high", "censor"] == ">ULOQ" and res.loc["high", "censor_limit"] == 2000
    assert res.loc["noisy", "censor"] == "cv_fail"
    assert res.loc["mixed", "censor"] == "inconsistent"
    assert res[["conc"]].loc[["low", "high", "noisy", "mixed"]].isna().all().all()


def _bridge_results(bias: dict[str, float], levels=(40.0, 160.0, 600.0), noise=0.0, seed=0):
    rng = np.random.default_rng(seed)
    rows = []
    for plate, k in bias.items():
        for i, lvl in enumerate(levels):
            rows.append(
                {
                    "plate_id": plate,
                    "analyte": "X",
                    "type": "control",
                    "sample_id": f"C{i}",
                    "conc": lvl / k * np.exp(rng.normal(0, noise)),
                    "censor": "",
                    "censor_limit": np.nan,
                }
            )
        rows.append(
            {
                "plate_id": plate,
                "analyte": "X",
                "type": "sample",
                "sample_id": f"S_{plate}",
                "conc": 100.0 / k,
                "censor": "",
                "censor_limit": np.nan,
            }
        )
    return pd.DataFrame(rows)


def test_plate_factors_recover_known_bias():
    bias = {"P1": 1.2, "P2": 0.8, "P3": 1.0}
    res = _bridge_results(bias)
    f = normalize.estimate_plate_factors(res).set_index("plate_id")
    gm = np.exp(np.mean(np.log(list(bias.values()))))
    for p, k in bias.items():
        assert f.loc[p, "factor"] == pytest.approx(gm / k, rel=1e-9)
        assert f.loc[p, "n_bridges"] == 3
    out = normalize.apply_plate_factors(res, f.reset_index())
    ctrl = out[out.type == "control"]
    assert ctrl.groupby("sample_id")["conc_norm"].std().max() == pytest.approx(0, abs=1e-9)
    s = out[out.type == "sample"]
    assert s["conc_norm"].std() == pytest.approx(0, abs=1e-9)


def test_normalization_reduces_bridge_cv_with_noise():
    res = _bridge_results({"P1": 1.25, "P2": 0.8, "P3": 1.1, "P4": 0.9}, noise=0.03, seed=2)
    f = normalize.estimate_plate_factors(res)
    out = normalize.apply_plate_factors(res, f)
    before = normalize.bridge_cv(out, "conc")["inter_plate_cv_pct"].mean()
    after = normalize.bridge_cv(out, "conc_norm")["inter_plate_cv_pct"].mean()
    assert after < before / 3
    assert f["log_factor_se"].notna().all()


def test_single_plate_and_missing_bridges():
    res = _bridge_results({"P1": 1.0})
    f = normalize.estimate_plate_factors(res)
    assert (f["factor"] == 1).all()
    res2 = _bridge_results({"P1": 1.2, "P2": 0.8})
    res2 = pd.concat(
        [
            res2,
            pd.DataFrame(
                [
                    {
                        "plate_id": "P3",
                        "analyte": "X",
                        "type": "sample",
                        "sample_id": "lonely",
                        "conc": 5.0,
                        "censor": "",
                        "censor_limit": np.nan,
                    }
                ]
            ),
        ],
        ignore_index=True,
    )
    with pytest.warns(UserWarning, match="no usable bridge"):
        f2 = normalize.estimate_plate_factors(res2).set_index("plate_id")
    assert f2.loc["P3", "factor"] == 1.0


def test_explicit_bridge_ids():
    res = _bridge_results({"P1": 1.2, "P2": 0.8})
    f = normalize.estimate_plate_factors(res, bridge_ids=["C1"]).set_index("plate_id")
    assert (f["n_bridges"] == 1).all()
    assert f.loc["P1", "factor"] / f.loc["P2", "factor"] == pytest.approx(0.8 / 1.2)
