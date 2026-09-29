"""Tests for group comparisons (immunoassay.stats)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import stats as sps

from immunoassay import stats as st


def _results(groups: dict[str, np.ndarray], analyte="X", censor=None) -> pd.DataFrame:
    rows = []
    for g, vals in groups.items():
        for i, v in enumerate(vals):
            rows.append(
                {
                    "plate_id": "P1",
                    "analyte": analyte,
                    "type": "sample",
                    "sample_id": f"{g}{i}",
                    "group": g,
                    "conc": v,
                    "conc_norm": v,
                    "censor": "",
                    "censor_limit": np.nan,
                    "censor_limit_norm": np.nan,
                }
            )
    df = pd.DataFrame(rows)
    for sid, (kind, limit) in (censor or {}).items():
        m = df.sample_id == sid
        df.loc[m, ["conc", "conc_norm"]] = np.nan
        df.loc[m, "censor"] = kind
        df.loc[m, ["censor_limit", "censor_limit_norm"]] = limit
    return df


def _lognormal(gm, n, seed, gsd=1.4):
    """Exactly log-normal-shaped sample (normal quantiles), shuffled by ``seed``."""
    q = sps.norm.ppf((np.arange(n) + 0.5) / n)
    return np.random.default_rng(seed).permutation(np.exp(np.log(gm) + np.log(gsd) * q))


def test_two_groups_welch_fold_change_matches_scipy():
    a, b = _lognormal(100, 12, 1), _lognormal(300, 12, 2)
    data = st.prepare_values(_results({"ctrl": a, "trt": b}))
    res = st.compare_groups(data)
    pw = res.pairwise.iloc[0]
    assert pw["test"] == "Welch t-test"
    ref = sps.ttest_ind(np.log10(b), np.log10(a), equal_var=False)
    assert pw["p"] == pytest.approx(ref.pvalue)
    gm_ratio = np.exp(np.mean(np.log(b)) - np.mean(np.log(a)))
    assert pw["estimate"] == pytest.approx(gm_ratio)
    assert pw["ci_low"] < gm_ratio < pw["ci_high"]
    assert pw["std_effect_name"] == "Hedges' g" and pw["std_effect"] > 0
    assert res.omnibus.iloc[0]["p_adj_bh"] == pytest.approx(pw["p"])


def test_non_normal_uses_mann_whitney_with_hodges_lehmann():
    a = np.array([1, 1.1, 1.2, 1.3, 1.1, 1.2, 1.0, 50, 60, 70.0])  # log scale is bimodal
    b = a * 3
    res = st.compare_groups(st.prepare_values(_results({"a": a, "b": b})))
    pw = res.pairwise.iloc[0]
    assert pw["test"] == "Mann-Whitney U"
    assert pw["effect_name"] == "fold change (HL)"
    assert pw["estimate"] == pytest.approx(3.0)
    assert -1 <= pw["std_effect"] <= 1


def test_three_groups_welch_anova_and_holm():
    g = {
        "veh": _lognormal(80, 12, 3),
        "lps": _lognormal(400, 12, 4),
        "drug": _lognormal(180, 12, 5),
    }
    res = st.compare_groups(st.prepare_values(_results(g)), group_order=["veh", "lps", "drug"])
    om = res.omnibus.iloc[0]
    assert om["test"] == "Welch ANOVA" and om["p"] < 1e-4
    assert 0 < om["omnibus_effect"] < 1
    assert len(res.pairwise) == 3
    assert list(res.pairwise["reference"]) == ["veh", "veh", "lps"]
    assert (res.pairwise["p_holm"] >= res.pairwise["p"]).all()


def test_bh_across_analytes():
    df = pd.concat(
        [
            _results({"c": _lognormal(100, 10, 6), "t": _lognormal(400, 10, 7)}, analyte="A"),
            _results({"c": _lognormal(100, 10, 8), "t": _lognormal(105, 10, 9)}, analyte="B"),
        ]
    )
    res = st.compare_groups(st.prepare_values(df))
    om = res.omnibus.set_index("analyte")
    assert (om["p_adj_bh"] >= om["p"]).all()
    assert bool(om.loc["A", "significant"]) and not bool(om.loc["B", "significant"])


def test_censoring_policies_and_failed_plates():
    df = _results(
        {"c": np.array([10, 12, 14.0]), "t": np.array([20, 25, 30.0])},
        censor={"c0": ("<LLOQ", 8.0), "t0": (">ULOQ", 40.0), "t1": ("cv_fail", np.nan)},
    )
    sub = st.prepare_values(df).set_index("sample_id")
    assert sub.loc["c0", "value"] == pytest.approx(8 / np.sqrt(2))
    assert sub.loc["t0", "value"] == 40
    assert np.isnan(sub.loc["t1", "value"]) and sub.loc["t1", "handling"] == "excluded: cv_fail"
    exc = st.prepare_values(df, cfg=st.StatsConfig(censored="exclude")).set_index("sample_id")
    assert np.isnan(exc.loc["c0", "value"])
    qc = pd.DataFrame({"plate_id": ["P1"], "analyte": ["X"], "passed": [False]})
    failed = st.prepare_values(df, qc)
    assert failed["value"].isna().all()
    assert (failed["handling"] == "excluded: plate failed QC").all()


def test_raw_scale_reports_mean_difference():
    a, b = np.array([1, 2, 3, 4, 5.0]), np.array([3, 4, 5, 6, 7.0])
    cfg = st.StatsConfig(transform="none")
    pw = st.compare_groups(st.prepare_values(_results({"a": a, "b": b}), cfg=cfg), cfg).pairwise
    assert pw.iloc[0]["effect_name"] == "mean difference"
    assert pw.iloc[0]["estimate"] == pytest.approx(2.0)


def test_hodges_lehmann_ci_contains_estimate():
    x, y = np.arange(10.0), np.arange(10.0) + 5
    est, lo, hi = st.hodges_lehmann(x, y)
    assert est == 5 and lo <= est <= hi


def test_too_few_samples_not_tested():
    res = st.compare_groups(
        st.prepare_values(_results({"a": np.array([1.0]), "b": np.array([2.0])}))
    )
    assert res.omnibus.iloc[0]["test"] == "not tested"
    assert res.pairwise.empty


def test_harmonized_lloq_prevents_cross_plate_misranking():
    df = _results(
        {"c": np.array([10.0, 12, 14]), "t": np.array([35.0, 40, 90])},
        censor={"c0": ("<LLOQ", 62.5)},
    )  # insensitive plate: c0 could be up to 62.5
    data = st.prepare_values(df).set_index("sample_id")
    common = 62.5 / np.sqrt(2)
    # Every value below 62.5 now shares the censored value, measured or not.
    for sid in ("c0", "c1", "c2", "t0", "t1"):
        assert data.loc[sid, "value"] == pytest.approx(common)
        assert "common LLOQ" in data.loc[sid, "handling"]
    assert data.loc["t2", "value"] == 90
    off = st.prepare_values(df, cfg=st.StatsConfig(harmonize_lloq=False)).set_index("sample_id")
    assert off.loc["c1", "value"] == 12


def test_extrapolate_policy_uses_sub_lloq_values():
    df = _results(
        {"c": np.array([10.0, 12]), "t": np.array([30.0, 40])},
        censor={"c0": ("<LLOQ", 15.0), "c1": ("<LLOQ", 15.0)},
    )
    df["conc_extrap"] = [4.0, np.nan, 30.0, 40.0]
    df["conc_extrap_norm"] = df["conc_extrap"]
    data = st.prepare_values(df, cfg=st.StatsConfig(censored="extrapolate")).set_index("sample_id")
    assert data.loc["c0", "value"] == 4.0
    assert data.loc["c0", "handling"].startswith("extrapolated")
    assert data.loc["c1", "value"] == 2.0  # below curve floor: half the smallest value
    assert "lowest rank" in data.loc["c1", "handling"]


def test_all_tied_values_are_not_tested_instead_of_crashing():
    """All values censored at one common LLOQ: SciPy 1.14's kruskal raises on this."""
    tied = {g: np.full(5, 7.0) for g in ("a", "b", "c")}
    res = st.compare_groups(st.prepare_values(_results(tied)))
    assert res.omnibus.iloc[0]["test"] == "not tested"
    assert "identical" in res.omnibus.iloc[0]["note"]
    assert res.pairwise.empty


def test_single_tied_pair_gets_p_one():
    g = {"a": np.full(5, 7.0), "b": np.full(5, 7.0), "c": _lognormal(40, 8, 1)}
    res = st.compare_groups(st.prepare_values(_results(g)))
    pw = res.pairwise.set_index(["reference", "group"])
    assert pw.loc[("a", "b"), "p"] == 1.0
    assert pw.loc[("a", "b"), "test"].startswith("not tested")
