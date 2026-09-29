"""Dilution correction, dilution linearity, and choosing the reportable result.

Samples are diluted before assay for two reasons. The analyte must land in
the quantifiable part of the curve, and serum/plasma *matrix* (other
proteins, lipids, heterophilic antibodies) interferes less when diluted. The
in-well concentration is converted back to the **neat** (undiluted)
concentration by multiplying by the fold-dilution.

Dilution linearity
------------------
If an assay behaves ideally, a sample run at 1:2 and at 1:8 gives the same
*neat* concentration after correction. When the corrected values disagree
by more than 20%, something other than the analyte is shaping the signal.

* **Corrected concentration rises with dilution** is the classic signature
  of *matrix interference*: something in the sample suppresses signal, and
  diluting it away "unmasks" analyte. It is also the signature of the
  **high-dose hook effect**. At extreme concentrations, free analyte
  saturates the capture and detection antibodies *separately*, fewer
  sandwiches form, and the signal paradoxically drops. A hooked sample can
  even look in range at low dilution.
* **Corrected concentration falls with dilution** suggests the sample's
  analyte and the calibrator do not react identically (non-parallelism),
  for example different isoforms or binding proteins. It can also mean that
  the most-diluted result sits near the LLOQ, where accuracy is weakest.

In both cases the most-diluted in-range result is the most trustworthy,
because it carries the least matrix and is furthest from any hook. It is
reported, and the sample is flagged.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

SAMPLE_KEYS = ["plate_id", "analyte", "type", "sample_id"]


def apply_dilution(samples: pd.DataFrame) -> pd.DataFrame:
    """Add neat concentrations: ``conc_neat``, ``lloq_neat``, ``uloq_neat``.

    The quantification limits scale with dilution too. A 1:4 sample's LLOQ
    is 4 × the in-well LLOQ, and that is the censoring value used when the
    sample reads below range.
    """
    out = samples.copy()
    out["conc_neat"] = out["conc_well"] * out["dilution"]
    out["lloq_neat"] = out["lloq_well"] * out["dilution"]
    out["uloq_neat"] = out["uloq_well"] * out["dilution"]
    return out


def _in_range(df: pd.DataFrame) -> pd.Series:
    return df["range_flag"].eq("") & ~df["cv_high"].astype(bool) & df["conc_neat"].notna()


def dilution_linearity(samples: pd.DataFrame, max_deviation_pct: float = 20.0) -> pd.DataFrame:
    """Compare neat concentrations of samples measured at two or more dilutions.

    Only dilutions that are in range with acceptable CV are compared.
    Out-of-range dilutions are expected, since that is why samples get run
    at several.

    Returns one row per multiply-diluted sample with the columns:
    ``dilutions`` and ``neat_concs`` (``;``-joined), ``n_in_range``,
    ``max_deviation_pct`` (largest deviation from the median of in-range
    results), ``trend`` (``rises_with_dilution`` / ``falls_with_dilution``
    / ``flat``), ``linear`` (bool, NaN if < 2 in range), and
    ``interpretation``.
    """
    df = samples[samples["type"] == "sample"]
    rows = []
    for key, grp in df.groupby(SAMPLE_KEYS, sort=False):
        if grp["dilution"].nunique() < 2:
            continue
        grp = grp.sort_values("dilution")
        ok = grp[_in_range(grp)]
        record = dict(zip(SAMPLE_KEYS, key, strict=True))
        record.update(
            group=grp["group"].iloc[0],
            dilutions=";".join(f"{d:g}" for d in grp["dilution"]),
            neat_concs=";".join(
                "NA" if not np.isfinite(c) else f"{c:.4g}" for c in grp["conc_neat"]
            ),
            n_in_range=len(ok),
        )
        if len(ok) < 2:
            hooked = grp["range_flag"].iloc[-1] == "" and grp["range_flag"].iloc[0] == "above_uloq"
            record.update(
                max_deviation_pct=np.nan,
                trend="",
                linear=np.nan,
                interpretation=(
                    "only the higher dilution quantifies (expected for high samples)"
                    if hooked
                    else "fewer than two dilutions in range; linearity not assessable"
                ),
            )
            rows.append(record)
            continue
        conc = ok["conc_neat"].to_numpy()
        med = float(np.median(conc))
        dev = float(np.max(np.abs(conc - med)) / med * 100)
        change = (conc[-1] - conc[0]) / conc[0] * 100
        trend = (
            "rises_with_dilution"
            if change > max_deviation_pct
            else "falls_with_dilution"
            if change < -max_deviation_pct
            else "flat"
        )
        linear = dev <= max_deviation_pct and trend == "flat"
        interpretation = {
            "rises_with_dilution": "possible matrix interference or high-dose hook effect",
            "falls_with_dilution": "possible non-parallelism (analyte vs calibrator)",
            "flat": "dilution-linear" if linear else "scatter between dilutions",
        }[trend]
        record.update(
            max_deviation_pct=dev, trend=trend, linear=linear, interpretation=interpretation
        )
        rows.append(record)
    cols = [
        *SAMPLE_KEYS,
        "group",
        "dilutions",
        "neat_concs",
        "n_in_range",
        "max_deviation_pct",
        "trend",
        "linear",
        "interpretation",
    ]
    return pd.DataFrame(rows, columns=cols)


def reportable_results(
    samples: pd.DataFrame, linearity: pd.DataFrame | None = None
) -> pd.DataFrame:
    """Collapse dilutions into one reportable result per (plate, analyte, sample).

    Rules:

    1. If ≥ 1 dilution is in range (quantifiable, CV OK), the result is the
       mean neat concentration of the in-range dilutions. If linearity
       failed, the most-diluted in-range result is used instead (see module
       docstring).
    2. If none is in range, the result is **censored**: ``censor = "<LLOQ"``
       with ``censor_limit`` = the smallest neat LLOQ among the dilutions
       (the most sensitive measurement), or ``">ULOQ"`` with the largest
       neat ULOQ. The numeric ``conc`` stays NaN. How a censored value
       enters statistics is a separate, explicit choice made in
       :mod:`immunoassay.stats`.
    3. Groups that failed CV at every dilution get ``censor = "cv_fail"``.

    ``conc_extrap`` is the reportable value when there is one. Otherwise it is
    the mean back-calculated neat concentration *outside* the validated range
    (still on the curve, so below LLOQ but above the curve floor). It is
    never used by default. It exists for sensitivity analyses that reproduce
    how many labs, and many published papers, report sub-LLOQ values.
    """
    lin = (
        linearity.set_index(SAMPLE_KEYS)["linear"]
        if linearity is not None and not linearity.empty
        else pd.Series(dtype=object)
    )
    rows = []
    for key, grp in samples.groupby(SAMPLE_KEYS, sort=False):
        grp = grp.sort_values("dilution")
        ok = grp[_in_range(grp)]
        linear = lin.get(key, np.nan) if len(lin) else np.nan
        censor, limit = "", np.nan
        if len(ok):
            if linear is False:
                conc = float(ok["conc_neat"].iloc[-1])
                basis = f"1:{ok['dilution'].iloc[-1]:g} (non-linear; most dilute used)"
            else:
                conc = float(ok["conc_neat"].mean())
                basis = " & ".join(f"1:{d:g}" for d in ok["dilution"])
        else:
            conc, basis = np.nan, ""
            flags = set(grp["range_flag"])
            if "" in flags:
                censor = "cv_fail"  # quantifiable at some dilution, but imprecise
            elif flags == {"below_lloq"}:
                censor, limit = "<LLOQ", float(grp["lloq_neat"].min())
            elif flags == {"above_uloq"}:
                censor, limit = ">ULOQ", float(grp["uloq_neat"].max())
            elif flags == {"no_range"}:
                censor = "no_range"
            else:
                censor = "inconsistent"  # e.g. >ULOQ at 1:2 but <LLOQ at 1:8
        finite = grp.loc[grp["conc_neat"].notna(), "conc_neat"]
        conc_extrap = (
            conc if np.isfinite(conc) else (float(finite.mean()) if len(finite) else np.nan)
        )
        rows.append(
            {
                **dict(zip(SAMPLE_KEYS, key, strict=True)),
                "group": grp["group"].iloc[0],
                "conc": conc,
                "conc_extrap": conc_extrap,
                "censor": censor,
                "censor_limit": limit,
                "basis": basis,
                "n_dilutions": len(grp),
                "flags": ";".join(sorted({f for fl in grp["flags"] for f in fl.split(";") if f})),
                "linear": linear,
            }
        )
    return pd.DataFrame(rows)
