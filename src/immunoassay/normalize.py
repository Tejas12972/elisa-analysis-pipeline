"""Inter-plate normalization using shared bridge samples.

Each plate carries its own standard curve, which corrects plate-wide changes
in *signal*, such as a brighter substrate lot or longer development. It
cannot correct errors in the *standard itself*. If the calibrator on plate 2
was reconstituted 15% too concentrated, every sample on plate 2 reads 15%
low, and the curve on that plate looks perfect. The same happens with
anything that affects samples but not standards (sample thaw handling,
plate-position drift).

The remedy is to run the same **bridge samples** (inter-plate controls or
pooled reference samples) on every plate. Any systematic difference in
their measured concentrations between plates is a plate effect.

Model
-----
Effects are multiplicative, so the model works on the log scale::

    log(conc[p, s]) = mu[s] + beta[p] + noise

where ``mu[s]`` is the bridge sample's true level and ``beta[p]`` the plate's
log bias, with ``sum(beta) = 0``, so the *average* plate is the reference.
Ordinary least squares on this two-way layout handles plates missing some
bridges. Each plate's correction factor is ``exp(beta[p])``, and corrected
concentrations are ``conc / factor``. With one bridge per plate this is
simple ratio scaling. With several, plate effects are averaged over
bridges, which reduces the noise of the correction by about ``1/sqrt(n)``.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd


def estimate_plate_factors(
    results: pd.DataFrame,
    *,
    bridge_ids: list[str] | None = None,
    bridge_type: str = "control",
) -> pd.DataFrame:
    """Estimate a multiplicative correction factor per plate, per analyte.

    Args:
        results: Reportable results (``plate_id, analyte, type, sample_id, conc``).
        bridge_ids: Sample IDs to use as bridges. Default: every sample of
            ``bridge_type`` that was quantified on at least two plates.
        bridge_type: Which ``type`` supplies bridges by default.

    Returns:
        ``plate_id, analyte, factor, log_factor_se, n_bridges, bridge_ids``.
        Plates with no usable bridge get ``factor = 1`` and a warning. If an
        analyte has only one plate, all its factors are 1 (nothing to align).
    """
    rows = []
    for analyte, df in results.groupby("analyte", sort=False):
        plates = list(dict.fromkeys(df["plate_id"]))
        usable = df[np.isfinite(df["conc"]) & (df["conc"] > 0)]
        if bridge_ids is not None:
            usable = usable[usable["sample_id"].isin(bridge_ids)]
        else:
            usable = usable[usable["type"] == bridge_type]
        counts = usable.groupby("sample_id")["plate_id"].nunique()
        usable = usable[usable["sample_id"].isin(counts[counts >= 2].index)]

        bridged = sorted(set(usable["plate_id"]), key=plates.index)
        for p in plates:
            if p not in bridged and len(plates) > 1:
                warnings.warn(
                    f"{analyte}: plate {p!r} has no usable bridge samples; left uncorrected",
                    UserWarning,
                    stacklevel=2,
                )
        if len(bridged) < 2:
            rows += [
                {
                    "plate_id": p,
                    "analyte": analyte,
                    "factor": 1.0,
                    "log_factor_se": np.nan,
                    "n_bridges": 0,
                    "bridge_ids": "",
                }
                for p in plates
            ]
            continue

        samples = sorted(usable["sample_id"].unique())
        n_s, n_p = len(samples), len(bridged)
        # Design: one column per bridge level (mu) + (n_p - 1) sum-to-zero plate columns.
        x = np.zeros((len(usable), n_s + n_p - 1))
        y = np.log(usable["conc"].to_numpy(dtype=float))
        for i, (sid, pid) in enumerate(zip(usable["sample_id"], usable["plate_id"], strict=True)):
            x[i, samples.index(sid)] = 1.0
            j = bridged.index(pid)
            if j < n_p - 1:
                x[i, n_s + j] = 1.0
            else:
                x[i, n_s:] = -1.0
        coef, *_ = np.linalg.lstsq(x, y, rcond=None)
        beta = np.append(coef[n_s:], -coef[n_s:].sum())
        dof = len(y) - x.shape[1]
        se = np.full(n_p, np.nan)
        if dof > 0:
            resid = y - x @ coef
            s2 = float(resid @ resid) / dof
            cov = s2 * np.linalg.pinv(x.T @ x)
            contrast = np.zeros((n_p, x.shape[1]))
            for j in range(n_p):
                if j < n_p - 1:
                    contrast[j, n_s + j] = 1.0
                else:
                    contrast[j, n_s:] = -1.0
            se = np.sqrt(np.einsum("ij,jk,ik->i", contrast, cov, contrast))
        for p in plates:
            if p in bridged:
                j = bridged.index(p)
                ids = usable.loc[usable["plate_id"] == p, "sample_id"]
                rows.append(
                    {
                        "plate_id": p,
                        "analyte": analyte,
                        "factor": float(np.exp(beta[j])),
                        "log_factor_se": float(se[j]),
                        "n_bridges": int(ids.nunique()),
                        "bridge_ids": ",".join(sorted(ids)),
                    }
                )
            else:
                rows.append(
                    {
                        "plate_id": p,
                        "analyte": analyte,
                        "factor": 1.0,
                        "log_factor_se": np.nan,
                        "n_bridges": 0,
                        "bridge_ids": "",
                    }
                )
    return pd.DataFrame(rows)


def apply_plate_factors(results: pd.DataFrame, factors: pd.DataFrame) -> pd.DataFrame:
    """Divide concentrations (and censoring limits) by each plate's factor.

    Adds ``plate_factor`` and ``conc_norm`` / ``censor_limit_norm``. The raw
    ``conc`` is kept, so the correction can always be audited.
    """
    f = factors.set_index(["plate_id", "analyte"])["factor"]
    out = results.copy()
    idx = pd.MultiIndex.from_frame(out[["plate_id", "analyte"]])
    out["plate_factor"] = f.reindex(idx).fillna(1.0).to_numpy()
    out["conc_norm"] = out["conc"] / out["plate_factor"]
    if "conc_extrap" in out:
        out["conc_extrap_norm"] = out["conc_extrap"] / out["plate_factor"]
    out["censor_limit_norm"] = out["censor_limit"] / out["plate_factor"]
    return out


def bridge_cv(results: pd.DataFrame, value: str) -> pd.DataFrame:
    """Between-plate %CV of each bridge sample for ``value`` (``conc`` or ``conc_norm``).

    This is the headline before/after metric. Normalization should shrink it.
    """
    ctrl = results[(results["type"] == "control") & results[value].notna()]
    return (
        ctrl.groupby(["analyte", "sample_id"])[value]
        .agg(lambda v: 100 * v.std(ddof=1) / v.mean() if len(v) > 1 else np.nan)
        .rename("inter_plate_cv_pct")
        .reset_index()
    )
