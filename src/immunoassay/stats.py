"""Group comparisons of analyte concentrations.

Decisions, and why
------------------
**Log scale by default.** Cytokine concentrations are right-skewed and
roughly log-normal between subjects, and biological effects act
multiplicatively ("LPS raised IL-6 five-fold"). Testing log10 concentrations
makes the data closer to normal and variances closer to equal. It also turns
a difference in means into an interpretable **fold change** (ratio of
geometric means) with a confidence interval. Set ``transform="none"`` to
test raw concentrations.

**Test choice.** Normality is checked per group with Shapiro-Wilk. With
n ≈ 10 per group Shapiro-Wilk has little power, so a non-significant result
is weak reassurance. The fallback is therefore rank-based rather than
"prove normality first".

* 2 groups, all normal: **Welch's t-test**. It does not assume equal
  variances and loses almost no power when they *are* equal, so there is no
  reason to pre-test variances and switch to Student's t.
* 2 groups, otherwise: **Mann-Whitney U**, with the Hodges-Lehmann shift
  estimate (the median of all pairwise differences) and its
  distribution-free CI as the effect size.
* 3+ groups, all normal: **Welch's one-way ANOVA**. 3+ groups, otherwise:
  **Kruskal-Wallis**. Pairwise follow-ups use the same two-group tests, with
  Holm correction within each analyte.

**Multiple analytes.** Testing IL-6, TNF-α, IL-1β, ... at α = 0.05 each
inflates false discoveries. The omnibus p-values across analytes are
adjusted with **Benjamini-Hochberg**, which controls the expected
proportion of false discoveries (FDR). That is the usual standard for
biomarker panels, and it is less conservative than Bonferroni.

**Effect sizes.** A p-value only says an effect is unlikely to be zero.
Every comparison also reports how *big* it is, with a CI:

* fold change (log scale) or mean difference (raw scale);
* Hedges' g (standardised difference, small-sample corrected) for the
  parametric path;
* rank-biserial correlation for Mann-Whitney;
* η² (ANOVA) or ε² (Kruskal-Wallis) for the omnibus tests.

Handling flagged samples
------------------------
Every sample's fate is recorded in a ``handling`` column, so the report can
show exactly what entered the tests.

* **Plate failed QC**: excluded. A plate failing QC means its values are
  not trustworthy.
* **cv_fail / inconsistent / no_range**: excluded. The value is unknown.
* **<LLOQ**: substituted with ``LLOQ / √2`` by default (the conventional
  midpoint on the log scale for values known to lie in (0, LLOQ)). Dropping
  them would *remove the lowest values* and bias the group upwards, a
  worse error than substitution. Rank tests are unaffected by the exact
  substitute, because a censored value is simply the lowest rank.
  ``censored="exclude"`` is available and is reported as such.
* **>ULOQ**: substituted with ULOQ (a lower bound) and flagged. In a real
  study these samples must be re-assayed at a higher dilution.
* **Plates with different LLOQs** (``harmonize_lloq=True``, the default).
  A "<LLOQ" on a plate whose LLOQ is 62.5 pg/mL may hide a value *higher*
  than a measured 35 pg/mL on a more sensitive plate. Substituting each
  plate's own LLOQ/√2 would then rank samples wrongly across plates. So
  every value below the *study-wide* LLOQ (the largest censoring limit among
  the samples compared) is censored at that one common limit. This is
  conservative: it discards resolution the better plates had, so that no
  comparison depends on which plate a sample happened to be on.
* ``censored="extrapolate"`` uses the curve's back-calculated value below
  LLOQ (``conc_extrap``) instead. That is common practice, and what many
  papers do implicitly, but those values lie outside the validated range.
  Use it only as a labelled sensitivity analysis.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import Literal

import numpy as np
import pandas as pd
from scipy import stats as sps
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.oneway import anova_oneway


@dataclass(frozen=True)
class StatsConfig:
    """Statistics options (see the module docstring for the reasoning)."""

    transform: Literal["log10", "none"] = "log10"
    censored: Literal["substitute", "exclude", "extrapolate"] = "substitute"
    harmonize_lloq: bool = True
    lloq_divisor: float = float(np.sqrt(2))
    alpha: float = 0.05
    normality_alpha: float = 0.05
    exclude_failed_plates: bool = True
    value_col: str = "conc_norm"


@dataclass
class StatsResult:
    omnibus: pd.DataFrame
    pairwise: pd.DataFrame
    data: pd.DataFrame


def prepare_values(
    results: pd.DataFrame,
    plate_qc: pd.DataFrame | None = None,
    cfg: StatsConfig | None = None,
) -> pd.DataFrame:
    """Turn reportable results into analysis values plus a ``handling`` note per sample."""
    cfg = cfg or StatsConfig()
    df = results[(results["type"] == "sample") & results["group"].notna()].copy()
    df = df[df["group"].astype(str) != ""]
    value_col = cfg.value_col if cfg.value_col in df else "conc"
    limit_col = "censor_limit_norm" if value_col == "conc_norm" else "censor_limit"
    failed: set[tuple[str, str]] = set()
    if plate_qc is not None and cfg.exclude_failed_plates:
        failed = {(r.plate_id, r.analyte) for r in plate_qc.itertuples() if not bool(r.passed)}
    values, handling = [], []
    for r in df.itertuples():
        v, lim = getattr(r, value_col), getattr(r, limit_col)
        if (r.plate_id, r.analyte) in failed:
            values.append(np.nan)
            handling.append("excluded: plate failed QC")
        elif r.censor == "":
            values.append(v)
            handling.append("measured")
        elif r.censor == "<LLOQ" and cfg.censored == "extrapolate":
            ext = getattr(
                r, "conc_extrap_norm" if value_col == "conc_norm" else "conc_extrap", np.nan
            )
            values.append(ext)
            handling.append(
                "extrapolated: <LLOQ (outside validated range)"
                if np.isfinite(ext)
                else "below curve floor"
            )
        elif r.censor == "<LLOQ" and cfg.censored == "substitute":
            values.append(lim / cfg.lloq_divisor)
            handling.append(f"substituted: <LLOQ → LLOQ/{cfg.lloq_divisor:.3g}")
        elif r.censor == ">ULOQ" and cfg.censored == "substitute":
            values.append(lim)
            handling.append("substituted: >ULOQ → ULOQ (lower bound)")
        else:
            values.append(np.nan)
            handling.append(f"excluded: {r.censor}")
    df["value"] = values
    df["handling"] = handling
    if cfg.harmonize_lloq and cfg.censored == "substitute":
        df = _harmonize_lloq(df, limit_col, cfg)
    if cfg.censored == "extrapolate":
        # Below the curve floor no concentration is defined, but the sample is certainly
        # lower than every extrapolated one. Rank it lowest: half the smallest value.
        for _, grp in df.groupby("analyte"):
            floor = grp.index[grp["handling"] == "below curve floor"]
            pos = grp.loc[grp["value"] > 0, "value"]
            if len(floor) and len(pos):
                df.loc[floor, "value"] = float(pos.min()) / 2
                df.loc[floor, "handling"] = "below curve floor → lowest rank (min/2)"
    return df


def _harmonize_lloq(df: pd.DataFrame, limit_col: str, cfg: StatsConfig) -> pd.DataFrame:
    """Censor every value below the study-wide LLOQ at that common limit (per analyte)."""
    df = df.copy()
    for _, grp in df.groupby("analyte"):
        cens = grp[(grp["censor"] == "<LLOQ") & ~grp["handling"].str.startswith("excluded")]
        if cens.empty:
            continue
        common = float(cens[limit_col].max())
        own = grp[limit_col].fillna(common)
        below = grp.index[
            grp["value"].notna()
            & ((grp["value"] < common) | (grp["censor"] == "<LLOQ"))
            & (own.to_numpy() <= common)
        ]
        df.loc[below, "value"] = common / cfg.lloq_divisor
        df.loc[below, "handling"] = [
            f"substituted: below common LLOQ {common:.3g} → /{cfg.lloq_divisor:.3g}" for _ in below
        ]
    return df


# --------------------------------------------------------------------------- #
# Effect-size helpers
# --------------------------------------------------------------------------- #
def hedges_g(x: np.ndarray, y: np.ndarray, alpha: float = 0.05) -> tuple[float, float, float]:
    """Hedges' g for ``y - x`` with an approximate normal-theory CI."""
    n1, n2 = len(x), len(y)
    sp = np.sqrt(((n1 - 1) * np.var(x, ddof=1) + (n2 - 1) * np.var(y, ddof=1)) / (n1 + n2 - 2))
    d = (np.mean(y) - np.mean(x)) / sp if sp > 0 else np.nan
    j = 1 - 3 / (4 * (n1 + n2) - 9)
    g = d * j
    se = np.sqrt((n1 + n2) / (n1 * n2) + g**2 / (2 * (n1 + n2)))
    z = sps.norm.ppf(1 - alpha / 2)
    return float(g), float(g - z * se), float(g + z * se)


def hodges_lehmann(x: np.ndarray, y: np.ndarray, alpha: float = 0.05) -> tuple[float, float, float]:
    """Hodges-Lehmann shift estimate of ``y - x`` with its distribution-free CI.

    The CI comes from the order statistics of the n1·n2 pairwise differences.
    The rank ``k`` is taken from the normal approximation to the Mann-Whitney
    null distribution (Bauer 1972).
    """
    diffs = np.sort(np.subtract.outer(y, x).ravel())
    n1, n2 = len(x), len(y)
    est = float(np.median(diffs))
    z = sps.norm.ppf(1 - alpha / 2)
    k = int(np.floor(n1 * n2 / 2 - z * np.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)))
    k = max(k, 0)
    lo, hi = diffs[k], diffs[len(diffs) - k - 1]
    return est, float(lo), float(hi)


def _two_group(x: np.ndarray, y: np.ndarray, parametric: bool, cfg: StatsConfig) -> dict:
    """Test ``y`` vs ``x`` (reference) and return test statistics plus effect sizes."""
    log = cfg.transform == "log10"
    if np.ptp(np.concatenate([x, y])) == 0:
        # Every value tied (typically: all censored at the same common LLOQ). There is no
        # evidence of a difference, and SciPy versions disagree on whether to raise here.
        return {
            "test": "not tested (all values tied)",
            "statistic": np.nan,
            "p": 1.0,
            "estimate": 1.0 if log else 0.0,
            "ci_low": np.nan,
            "ci_high": np.nan,
            "effect_name": "fold change" if log else "mean difference",
            "std_effect_name": "",
            "std_effect": np.nan,
            "std_effect_ci": "",
        }
    if parametric:
        res = sps.ttest_ind(y, x, equal_var=False)
        ci = res.confidence_interval(1 - cfg.alpha)
        diff = float(np.mean(y) - np.mean(x))
        g, g_lo, g_hi = hedges_g(x, y, cfg.alpha)
        out = {
            "test": "Welch t-test",
            "statistic": float(res.statistic),
            "p": float(res.pvalue),
            "estimate": diff,
            "ci_low": float(ci.low),
            "ci_high": float(ci.high),
            "std_effect_name": "Hedges' g",
            "std_effect": g,
            "std_effect_ci": f"[{g_lo:.2f}, {g_hi:.2f}]",
        }
    else:
        res = sps.mannwhitneyu(y, x, alternative="two-sided")
        est, lo, hi = hodges_lehmann(x, y, cfg.alpha)
        rbc = 2 * float(res.statistic) / (len(x) * len(y)) - 1
        out = {
            "test": "Mann-Whitney U",
            "statistic": float(res.statistic),
            "p": float(res.pvalue),
            "estimate": est,
            "ci_low": lo,
            "ci_high": hi,
            "std_effect_name": "rank-biserial r",
            "std_effect": rbc,
            "std_effect_ci": "",
        }
    if log:
        out["effect_name"] = "fold change" if parametric else "fold change (HL)"
        for k in ("estimate", "ci_low", "ci_high"):
            out[k] = float(10 ** out[k])
    else:
        out["effect_name"] = "mean difference" if parametric else "median shift (HL)"
    return out


# --------------------------------------------------------------------------- #
# Main entry point
# --------------------------------------------------------------------------- #
def compare_groups(
    data: pd.DataFrame,
    cfg: StatsConfig | None = None,
    group_order: list[str] | None = None,
) -> StatsResult:
    """Compare treatment groups for every analyte in ``data`` (from :func:`prepare_values`).

    The first group in ``group_order`` (default: order of appearance) is the
    reference. Fold changes are "group / reference".
    """
    cfg = cfg or StatsConfig()
    omni_rows, pair_rows = [], []
    for analyte, df in data.groupby("analyte", sort=False):
        order = group_order or list(dict.fromkeys(df["group"]))
        vals: dict[str, np.ndarray] = {}
        for g in order:
            v = df.loc[(df["group"] == g) & np.isfinite(df["value"]), "value"].to_numpy(float)
            if cfg.transform == "log10":
                v = np.log10(v[v > 0])
            if len(v) >= 2:
                vals[g] = v
        groups = list(vals)
        if len(groups) < 2:
            omni_rows.append(
                {
                    "analyte": analyte,
                    "test": "not tested",
                    "note": "fewer than two groups with n >= 2",
                }
            )
            continue

        if np.ptp(np.concatenate(list(vals.values()))) == 0:
            omni_rows.append(
                {
                    "analyte": analyte,
                    "test": "not tested",
                    "note": "all values identical (e.g. all censored at the common LLOQ)",
                }
            )
            continue

        shapiro = {
            g: float(sps.shapiro(v).pvalue) if len(v) >= 3 and np.ptp(v) > 0 else np.nan
            for g, v in vals.items()
        }
        normal = all(np.isfinite(p) and p > cfg.normality_alpha for p in shapiro.values())
        norm_note = ", ".join(f"{g}: p={p:.2g}" for g, p in shapiro.items())
        ns = {g: len(v) for g, v in vals.items()}
        base = {
            "analyte": analyte,
            "groups": " vs ".join(groups),
            "n": ", ".join(f"{g}={n}" for g, n in ns.items()),
            "shapiro": norm_note,
            "normal": normal,
            "scale": cfg.transform,
        }

        if len(groups) == 2:
            r = _two_group(vals[groups[0]], vals[groups[1]], normal, cfg)
            omni_rows.append(
                {
                    **base,
                    "test": r["test"],
                    "statistic": r["statistic"],
                    "p": r["p"],
                    "omnibus_effect_name": r["std_effect_name"],
                    "omnibus_effect": r["std_effect"],
                }
            )
        else:
            samples = [vals[g] for g in groups]
            n_tot = sum(ns.values())
            grand = np.concatenate(samples)
            if normal:
                res = anova_oneway(samples, use_var="unequal")
                ss_b = sum(len(s) * (s.mean() - grand.mean()) ** 2 for s in samples)
                ss_t = float(((grand - grand.mean()) ** 2).sum())
                omni_rows.append(
                    {
                        **base,
                        "test": "Welch ANOVA",
                        "statistic": float(res.statistic),
                        "p": float(res.pvalue),
                        "omnibus_effect_name": "eta²",
                        "omnibus_effect": ss_b / ss_t if ss_t > 0 else np.nan,
                    }
                )
            else:
                h, p = sps.kruskal(*samples)
                omni_rows.append(
                    {
                        **base,
                        "test": "Kruskal-Wallis",
                        "statistic": float(h),
                        "p": float(p),
                        "omnibus_effect_name": "epsilon²",
                        "omnibus_effect": float(h / (n_tot - 1)),
                    }
                )

        pairs = list(combinations(groups, 2))
        results = [(_two_group(vals[a], vals[b], normal, cfg), a, b) for a, b in pairs]
        raw_p = [r["p"] for r, _, _ in results]
        holm = multipletests(raw_p, method="holm")[1] if len(raw_p) > 1 else raw_p
        for (r, a, b), ph in zip(results, holm, strict=True):
            pair_rows.append(
                {
                    "analyte": analyte,
                    "reference": a,
                    "group": b,
                    "n_ref": ns[a],
                    "n_group": ns[b],
                    "test": r["test"],
                    "effect_name": r["effect_name"],
                    "estimate": r["estimate"],
                    "ci_low": r["ci_low"],
                    "ci_high": r["ci_high"],
                    "std_effect_name": r["std_effect_name"],
                    "std_effect": r["std_effect"],
                    "std_effect_ci": r["std_effect_ci"],
                    "p": r["p"],
                    "p_holm": float(ph),
                    "significant": bool(ph < cfg.alpha),
                }
            )

    omnibus = pd.DataFrame(omni_rows)
    # Stable schema even when nothing could be tested (e.g. every value censored).
    for col in ("test", "note", "statistic", "p", "p_adj_bh"):
        if col not in omnibus:
            omnibus[col] = np.nan
    omnibus["significant"] = False
    if omnibus["p"].notna().any():
        tested = omnibus["p"].notna()
        omnibus.loc[tested, "p_adj_bh"] = multipletests(omnibus.loc[tested, "p"], method="fdr_bh")[
            1
        ]
        omnibus["significant"] = omnibus["p_adj_bh"] < cfg.alpha
    return StatsResult(omnibus=omnibus, pairwise=pd.DataFrame(pair_rows), data=data)
