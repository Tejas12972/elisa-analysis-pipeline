"""Quality control: outliers, replicate precision, assay limits, controls, plate acceptance.

Definitions used throughout
---------------------------
**%CV** (coefficient of variation) is ``100 × SD / mean`` of replicate wells.
It measures precision in a way that does not depend on scale. Ligand-binding
assay guidance (FDA *Bioanalytical Method Validation*, 2018; ICH M10, 2022)
treats ≤ 20% as acceptable. Where every replicate quantifies, CV is computed
on back-calculated concentration, because that is the reported quantity.
Otherwise it falls back to OD.

**% recovery** is ``100 × back-calculated / nominal`` for a standard. It is
how well the fitted curve reproduces the concentrations it was built from.
The acceptance window is 80–120%.

**LOD** (limit of detection) is the smallest concentration distinguishable
from zero: mean blank + 3 SD of the blanks, converted to concentration
through the curve. After blank subtraction the mean blank is 0, so the OD
threshold is simply ``3 × blank_sd``. If the fitted curve's floor lies above
that threshold, the floor + 3 SD is used, and ``lod_basis`` records which
rule applied. With only 2-3 blank wells the SD is very uncertain, so treat
the LOD as an order of magnitude. A signal above the LOD says "analyte
is present", but not *how much* with acceptable accuracy.

**LLOQ / ULOQ** (lower / upper limit of quantification) bound the range where
concentrations are both accurate and precise. Here they are the lowest and
highest standards of the widest *contiguous* run of levels in which every
interior level has recovery within 80–120% and CV ≤ 20%. The two end levels
(the LLOQ and ULOQ themselves) may use the relaxed 75–125% / CV ≤ 25%
allowed at the limits by FDA (2018) and ICH M10. Requiring contiguity stops
one lucky low-end standard from extending the range across a failing level.
Samples outside [LLOQ, ULOQ] are flagged, not reported as numbers. Below
LLOQ they are censored (see :mod:`immunoassay.stats`). Above ULOQ they
should be re-assayed at a higher dilution.

Outlier detection: why a median/MAD rule and not Grubbs
-------------------------------------------------------
ELISA replicates come in twos or threes. Grubbs' test assumes normality,
and at n = 3 its critical value (1.153) is barely below the largest value
the statistic can reach (1.155), so it almost never fires. The modified
z-score of Iglewicz & Hoaglin (1993),
``M = 0.6745 · (x − median) / MAD`` with ``|M| > 3.5`` flagged, uses the
median and median absolute deviation. One bad well cannot inflate them. The
rule works from n = 3.

At n = 3 the MAD is simply the *smaller* of the two deviations from the
median, so on its own the rule fires on any slightly lopsided triplicate.
On the real ELISAtools plates the ungated rule flagged 23 wells, only 6 of
which the analyst had masked. The z-score test is therefore **gated on precision**: it is
only applied when the group's raw-OD CV already fails ``max_cv_pct``. An
"outlier" in a group that meets the precision criterion changes nothing
that matters, and removing it would only make the CV look better than it is.

With **duplicates**, no statistic can tell which of two disagreeing wells is
wrong. For samples the pair is therefore flagged ``cv_high`` rather than
thinned. For **standards** there is extra information, because the fitted
curve says where the well *should* be. A standard well is excluded when
three things hold: its level's CV fails, it is the member of the pair
farther from the curve, and its weighted residual is > 3.5 robust SDs from
zero. At most one well per level and ``max_std_exclusions`` per plate are
excluded, so QC cannot quietly rebuild a curve to fit itself.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from immunoassay.curves import CurveFit, compute_weights, fit_curve

REPLICATE_KEYS = ["plate_id", "analyte", "type", "sample_id", "dilution"]


@dataclass(frozen=True)
class QCConfig:
    """QC thresholds. Defaults follow FDA (2018) / ICH M10 ligand-binding-assay guidance."""

    max_cv_pct: float = 20.0
    recovery_low_pct: float = 80.0
    recovery_high_pct: float = 120.0
    edge_tolerance_pct: float = 5.0
    lod_sd_multiplier: float = 3.0
    outlier_z: float = 3.5
    max_std_exclusions: int = 2
    min_r_squared: float = 0.98
    min_std_pass_frac: float = 0.75
    min_control_pass_frac: float = 2 / 3


def _cv(values: np.ndarray) -> float:
    values = values[np.isfinite(values)]
    if values.size < 2 or values.mean() == 0:
        return np.nan
    return float(100 * values.std(ddof=1) / abs(values.mean()))


def modified_z(values: np.ndarray) -> np.ndarray:
    """Iglewicz-Hoaglin modified z-scores. All-zero when the data have no spread.

    MAD is floored at 1% of |median| (and 1e-3) so that three nearly identical
    replicates do not turn a trivial difference into a huge z-score.
    """
    values = np.asarray(values, dtype=float)
    med = np.nanmedian(values)
    mad = np.nanmedian(np.abs(values - med))
    mad = max(mad, 0.01 * abs(med), 1e-3)
    return 0.6745 * (values - med) / mad


# --------------------------------------------------------------------------- #
# Outliers
# --------------------------------------------------------------------------- #
def flag_replicate_outliers(wells: pd.DataFrame, cfg: QCConfig | None = None) -> pd.DataFrame:
    """Flag wells whose OD is a modified-z outlier within their replicate group (n >= 3).

    A group is only tested when its raw-OD CV exceeds ``max_cv_pct`` (see the
    module docstring). Adds ``outlier`` (bool) and ``outlier_reason``. Blank
    wells are included, because one bad blank shifts every net OD on the plate.
    """
    cfg = cfg or QCConfig()
    out = wells.copy()
    out["outlier"] = False
    out["outlier_reason"] = ""
    ok = out["od_status"] == "ok"
    for _, grp in out[ok].groupby(REPLICATE_KEYS, sort=False, dropna=False):
        if len(grp) < 3 or not _cv(grp["od"].to_numpy()) > cfg.max_cv_pct:
            continue
        z = modified_z(grp["od"].to_numpy())
        bad = grp.index[np.abs(z) > cfg.outlier_z]
        # Never discard a majority: with n=3 at most one well can go.
        if len(bad) and len(grp) - len(bad) >= 2:
            out.loc[bad, "outlier"] = True
            out.loc[bad, "outlier_reason"] = (
                f"replicate CV > {cfg.max_cv_pct:g}% and |modified z| > {cfg.outlier_z:g}"
            )
    return out


def flag_standard_residual_outliers(
    wells: pd.DataFrame,
    fits: dict[tuple[str, str], CurveFit],
    cfg: QCConfig | None = None,
) -> pd.DataFrame:
    """Exclude gross-error standard wells using leave-one-out residuals.

    Only wells in a level whose replicate CV fails are candidates. Each
    candidate is judged by its **deletion residual**: the curve is refitted
    *without* it, and its weighted residual from that curve is divided by the
    refit's residual SD. A plain residual from the full fit would understate
    the error, because a gross error drags the curve towards itself (and two
    of them can drag it far enough to hide both). The worst candidate above
    ``outlier_z`` is excluded, and the search repeats up to
    ``max_std_exclusions`` times per plate.

    Returns a copy of ``wells`` with ``outlier`` updated. Callers should refit
    the curves if any new standard wells were flagged.
    """
    cfg = cfg or QCConfig()
    out = wells.copy()
    if "outlier" not in out.columns:
        out["outlier"] = False
        out["outlier_reason"] = ""
    for (pid, analyte), fit in fits.items():
        for _ in range(cfg.max_std_exclusions):
            mask = (
                (out["plate_id"] == pid)
                & (out["analyte"] == analyte)
                & (out["type"] == "standard")
                & (out["od_status"] == "ok")
                & ~out["outlier"]
            )
            std = out[mask]
            conc = std["concentration"].to_numpy(dtype=float)
            y = std["od_net"].to_numpy(dtype=float)
            w = compute_weights(y, fit.weighting, fit.weight_offset)
            candidates = [
                i
                for _, lvl in std.groupby("concentration")
                if len(lvl) >= 2 and _cv(lvl["od_net"].to_numpy()) > cfg.max_cv_pct
                for i in lvl.index
            ]
            best_t, best_idx = 0.0, None
            for idx in candidates:
                pos = std.index.get_loc(idx)
                keep = np.arange(len(std)) != pos
                try:
                    refit = fit_curve(
                        conc[keep],
                        y[keep],
                        model="4pl",
                        weighting=fit.weighting,
                        weight_offset=fit.weight_offset,
                    )
                except (ValueError, RuntimeError):
                    continue
                r_keep = y[keep] - refit.predict(conc[keep])
                dof = keep.sum() - 4
                if dof <= 0:
                    continue
                s = np.sqrt(np.sum(w[keep] * r_keep**2) / dof)
                r_i = y[pos] - float(refit.predict(conc[pos]))
                t = abs(r_i) * np.sqrt(w[pos]) / max(s, 1e-12)
                if t > best_t:
                    best_t, best_idx = t, idx
            if best_idx is None or best_t <= cfg.outlier_z:
                break
            out.loc[best_idx, "outlier"] = True
            out.loc[best_idx, "outlier_reason"] = (
                f"standard pair CV > {cfg.max_cv_pct:g}% and deletion residual "
                f"{best_t:.1f} SD > {cfg.outlier_z:g}"
            )
    return out


def _level_passes(rec: pd.DataFrame, cfg: QCConfig) -> tuple[np.ndarray, np.ndarray]:
    """Per standard level: (strict pass, relaxed "edge" pass).

    Strict is recovery within ``recovery_low/high_pct`` and CV ≤ ``max_cv_pct``.
    Edge widens both by ``edge_tolerance_pct``. That is the FDA (2018) / ICH M10
    allowance for the LLOQ and ULOQ calibrators (±25% instead of ±20%). A level
    with a single quantified well has no CV and is judged on recovery only.
    """
    single = rec["n_quantified"].to_numpy() == 1
    cv = rec["cv_pct"].to_numpy(dtype=float)
    recv = rec["recovery_pct"].to_numpy(dtype=float)
    t = cfg.edge_tolerance_pct
    strict = (
        (recv >= cfg.recovery_low_pct)
        & (recv <= cfg.recovery_high_pct)
        & ((cv <= cfg.max_cv_pct) | single)
    )
    edge = (
        (recv >= cfg.recovery_low_pct - t)
        & (recv <= cfg.recovery_high_pct + t)
        & ((cv <= cfg.max_cv_pct + t) | single)
    )
    return strict, edge


def _quantification_range(strict: np.ndarray, edge: np.ndarray) -> tuple[int | None, int | None]:
    """Widest run of levels whose interior passes strictly and whose ends pass relaxed.

    Returns index positions (LLOQ, ULOQ) in ascending-concentration order, or
    ``(None, None)`` if no run of at least two levels qualifies. Ties go to the
    run with the lower LLOQ, because sensitivity is usually what is scarce.
    """
    best: tuple[int, int] | None = None
    n = len(strict)
    for i in range(n):
        if not edge[i]:
            continue
        for j in range(i + 1, n):
            if not all(strict[i + 1 : j]):
                break
            if edge[j] and (best is None or j - i > best[1] - best[0]):
                best = (i, j)
    return best if best is not None else (None, None)


# --------------------------------------------------------------------------- #
# Limits
# --------------------------------------------------------------------------- #
def plate_limits(
    wells: pd.DataFrame,
    fits: dict[tuple[str, str], CurveFit],
    recovery: pd.DataFrame,
    cfg: QCConfig | None = None,
) -> pd.DataFrame:
    """LOD, LLOQ, and ULOQ per (plate_id, analyte), in in-well concentration units.

    Also reports which standard levels passed and the fraction passing. The
    plate-acceptance rule uses that fraction.
    """
    cfg = cfg or QCConfig()
    rows = []
    for (pid, analyte), fit in fits.items():
        plate = wells[(wells["plate_id"] == pid) & (wells["analyte"] == analyte)]
        blank_sd = float(plate["blank_sd"].iloc[0])
        # After blank subtraction the mean blank is 0, so the threshold is k·SD. When the
        # fitted zero-concentration asymptote sits above that (common: the curve floor is
        # estimated from standards, not blanks), no concentration maps to the threshold, so
        # the curve floor is used as the zero-analyte response instead: floor + k·SD.
        lod_od, lod_basis, lod_conc = np.nan, "no blank SD", np.nan
        if np.isfinite(blank_sd):
            lod_od, lod_basis = cfg.lod_sd_multiplier * blank_sd, "mean blank + k·SD"
            floor = fit.params["a"]
            if fit.increasing and lod_od <= floor:
                lod_od, lod_basis = floor + cfg.lod_sd_multiplier * blank_sd, "curve floor + k·SD"
            conc, _ = fit.back_calculate(lod_od)
            lod_conc = float(conc[0])

        rec = recovery[(recovery["plate_id"] == pid) & (recovery["analyte"] == analyte)]
        rec = rec.sort_values("concentration")
        strict, edge = _level_passes(rec, cfg)
        levels = rec["concentration"].to_numpy()
        lloq_i, uloq_i = _quantification_range(strict, edge)
        lloq = float(levels[lloq_i]) if lloq_i is not None else np.nan
        uloq = float(levels[uloq_i]) if uloq_i is not None else np.nan
        # Calibrator acceptance (the 75% rule): interior levels at ±20%, the lowest and
        # highest calibrator at the relaxed ±25%.
        passing = strict.copy()
        if len(levels):
            passing[0] |= edge[0]
            passing[-1] |= edge[-1]
        rows.append(
            {
                "plate_id": pid,
                "analyte": analyte,
                "lod_od": lod_od,
                "lod": lod_conc,
                "lod_basis": lod_basis,
                "lloq": lloq,
                "uloq": uloq,
                "n_std_levels": len(levels),
                "n_std_levels_pass": int(passing.sum()),
                "std_pass_frac": float(passing.mean()) if len(levels) else np.nan,
                "failing_levels": ", ".join(f"{v:g}" for v in levels[~passing]),
            }
        )
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Replicate summaries and sample flags
# --------------------------------------------------------------------------- #
def summarize_replicates(
    wells: pd.DataFrame, limits: pd.DataFrame, cfg: QCConfig | None = None
) -> pd.DataFrame:
    """One row per replicate group of samples and controls, with QC flags.

    Columns include ``n_wells``, ``n_used`` (after outlier removal),
    ``od_mean``, ``od_cv_pct``, ``conc_well`` (mean in-well concentration),
    ``cv_pct`` (the CV that is judged), ``range_flag`` (``""``, ``below_lloq``,
    or ``above_uloq``), ``below_lod``, ``cv_high``, and ``flags`` (a
    human-readable ``;``-joined list).

    Range logic, in order of precedence:

    1. Any used well above the curve or saturated gives ``above_uloq``. The
       group cannot be quantified at this dilution.
    2. Any used well below the curve's floor gives ``below_lloq``.
    3. Otherwise the mean in-well concentration is compared with the plate's
       LLOQ/ULOQ. A plate without a valid quantification range marks
       everything ``no_range``.
    """
    cfg = cfg or QCConfig()
    lim = limits.set_index(["plate_id", "analyte"])
    subset = wells[wells["type"].isin(["sample", "control"])]
    rows = []
    for key, grp in subset.groupby(REPLICATE_KEYS, sort=False, dropna=False):
        pid, analyte, typ, sid, dilution = key
        used = grp[~grp["outlier"].astype(bool)] if "outlier" in grp else grp
        conc = used["conc_well"].to_numpy(dtype=float)
        flags_curve = set(used["curve_flag"])
        od_cv = _cv(used["od_net"].to_numpy(dtype=float))
        finite = conc[np.isfinite(conc)]
        conc_mean = float(finite.mean()) if finite.size == len(used) and finite.size else np.nan
        cv_pct = _cv(conc) if finite.size == len(used) and finite.size >= 2 else od_cv

        lloq, uloq, lod = (
            lim.loc[(pid, analyte), ["lloq", "uloq", "lod"]]
            if (pid, analyte) in lim.index
            else (np.nan, np.nan, np.nan)
        )
        if "above_curve" in flags_curve:
            range_flag = "above_uloq"
        elif "below_curve" in flags_curve:
            range_flag = "below_lloq"
        elif not (np.isfinite(lloq) and np.isfinite(uloq)):
            range_flag = "no_range"
        elif conc_mean < lloq:
            range_flag = "below_lloq"
        elif conc_mean > uloq:
            range_flag = "above_uloq"
        else:
            range_flag = ""
        below_lod = "below_curve" in flags_curve or (
            np.isfinite(lod) and np.isfinite(conc_mean) and conc_mean < lod
        )
        cv_high = bool(np.isfinite(cv_pct) and cv_pct > cfg.max_cv_pct)
        n_outliers = len(grp) - len(used)

        flags = [f for f in (range_flag,) if f]
        if below_lod:
            flags.append("below_lod")
        if cv_high:
            flags.append("cv_high")
        if n_outliers:
            flags.append(f"{n_outliers}_outlier_well_removed")
        rows.append(
            {
                "plate_id": pid,
                "analyte": analyte,
                "type": typ,
                "sample_id": sid,
                "group": grp["group"].iloc[0],
                "dilution": float(dilution),
                "wells": ",".join(grp["well"]),
                "n_wells": len(grp),
                "n_used": len(used),
                "od_mean": float(used["od_net"].mean()),
                "od_cv_pct": od_cv,
                "conc_well": conc_mean,
                "cv_pct": cv_pct,
                "lloq_well": lloq,
                "uloq_well": uloq,
                "range_flag": range_flag,
                "below_lod": bool(below_lod),
                "cv_high": cv_high,
                "flags": ";".join(flags),
            }
        )
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Controls + plate acceptance
# --------------------------------------------------------------------------- #
def check_controls(samples: pd.DataFrame, control_ranges: pd.DataFrame | None) -> pd.DataFrame:
    """Compare each control's neat concentration with its expected ``[low, high]``.

    ``control_ranges`` needs columns ``analyte, sample_id, low, high`` in neat
    concentration units. Controls without a specified range are reported with
    ``pass = NaN``. They still serve as inter-plate bridges.
    """
    ctrl = samples[samples["type"] == "control"].copy()
    if ctrl.empty:
        return ctrl.assign(low=np.nan, high=np.nan, passed=np.nan)
    if control_ranges is None or control_ranges.empty:
        ctrl["low"] = np.nan
        ctrl["high"] = np.nan
    else:
        ctrl = ctrl.merge(
            control_ranges[["analyte", "sample_id", "low", "high"]],
            on=["analyte", "sample_id"],
            how="left",
        )
    has_range = ctrl["low"].notna() & ctrl["high"].notna()
    value = ctrl["conc_neat"] if "conc_neat" in ctrl else ctrl["conc_well"] * ctrl["dilution"]
    within = value.between(ctrl["low"], ctrl["high"]) & ctrl["range_flag"].eq("")
    ctrl["passed"] = np.where(has_range, within, np.nan)
    return ctrl[
        [
            "plate_id",
            "analyte",
            "sample_id",
            "dilution",
            "conc_well",
            "conc_neat",
            "cv_pct",
            "range_flag",
            "low",
            "high",
            "passed",
        ]
    ]


def plate_summary(
    fits: dict[tuple[str, str], CurveFit],
    limits: pd.DataFrame,
    controls: pd.DataFrame,
    samples: pd.DataFrame,
    cfg: QCConfig | None = None,
) -> pd.DataFrame:
    """Per-plate pass/fail with the reasons for any failure.

    A plate passes when:

    * the curve's R² ≥ ``min_r_squared``;
    * ≥ ``min_std_pass_frac`` of standard levels meet recovery and CV
      criteria (the FDA/ICH "75% of calibrators" rule);
    * a quantification range (LLOQ < ULOQ) exists;
    * ≥ ``min_control_pass_frac`` of controls *with specified ranges* fall
      inside them (in the spirit of the "4-6-X" QC rule).
    """
    cfg = cfg or QCConfig()
    rows = []
    lim = limits.set_index(["plate_id", "analyte"])
    for (pid, analyte), fit in fits.items():
        reasons = []
        lim_row = lim.loc[(pid, analyte)]
        if not fit.r_squared >= cfg.min_r_squared:
            reasons.append(f"R² {fit.r_squared:.3f} < {cfg.min_r_squared}")
        if not lim_row["std_pass_frac"] >= cfg.min_std_pass_frac:
            reasons.append(
                f"only {lim_row['n_std_levels_pass']}/{lim_row['n_std_levels']} "
                "standard levels pass"
            )
        if not np.isfinite(lim_row["lloq"]):
            reasons.append("no valid quantification range")
        c = controls[(controls["plate_id"] == pid) & (controls["analyte"] == analyte)]
        judged = c[c["passed"].notna()]
        ctrl_frac = float(judged["passed"].astype(bool).mean()) if len(judged) else np.nan
        if len(judged) and ctrl_frac < cfg.min_control_pass_frac:
            failed = ", ".join(judged.loc[~judged["passed"].astype(bool), "sample_id"])
            reasons.append(f"controls out of range: {failed}")
        s = samples[
            (samples["plate_id"] == pid)
            & (samples["analyte"] == analyte)
            & (samples["type"] == "sample")
        ]
        rows.append(
            {
                "plate_id": pid,
                "analyte": analyte,
                "model": fit.model.upper(),
                "r_squared": fit.r_squared,
                "std_levels_pass": f"{lim_row['n_std_levels_pass']}/{lim_row['n_std_levels']}",
                "lod": lim_row["lod"],
                "lloq": lim_row["lloq"],
                "uloq": lim_row["uloq"],
                "controls_pass": (
                    f"{int(judged['passed'].astype(bool).sum())}/{len(judged)}"
                    if len(judged)
                    else "n/a"
                ),
                "n_sample_groups": len(s),
                "n_flagged": int((s["flags"] != "").sum()),
                "passed": not reasons,
                "reasons": "; ".join(reasons),
            }
        )
    return pd.DataFrame(rows)
