"""Standard-curve fitting (4PL / 5PL) and back-calculation of concentrations.

Why the curve is sigmoidal
--------------------------
A sandwich ELISA well holds a finite number of capture-antibody binding
sites. At very low analyte concentration almost no sandwiches form, so the
signal sits at a floor set by background. As the concentration rises the
signal climbs roughly in proportion to the log of the concentration. At
high concentration the capture sites saturate and the signal plateaus. On a
log-concentration axis this gives an "S" shape, which is described
empirically by the **four-parameter logistic (4PL)** model::

    y = d + (a - d) / (1 + (x / c) ** b)

========  ================================================================
``a``     response at zero concentration (lower asymptote, the "floor")
``d``     response at infinite concentration (upper asymptote, saturation)
``c``     inflection point, the concentration giving a response halfway
          between ``a`` and ``d`` (EC50). The curve is steepest and most
          precise here.
``b``     Hill slope, i.e. how steep the transition is
========  ================================================================

Real curves are often *asymmetric*: they approach the top more slowly than
the bottom. The **5PL** adds an asymmetry exponent ``g``::

    y = d + (a - d) / (1 + (x / c) ** b) ** g

With ``g == 1`` the 5PL reduces to the 4PL. The extra parameter can soak up
noise when only 7-8 standard levels are available. So the models are
compared with the small-sample AIC (AICc), and 5PL is chosen only when it
wins clearly (ΔAICc > 2).

Why weighted regression
-----------------------
ELISA noise is *heteroscedastic*: the SD of replicate ODs grows with the
signal, so the %CV is roughly constant. Ordinary least squares treats every
point as equally reliable. The high-OD standards, with the largest absolute
scatter, then dominate the fit and the low end of the curve is fitted badly.
Unfortunately the low end is where most biological samples fall. Weighting
each residual by ``1 / y²`` (equivalently ``sigma ∝ y``) assumes a constant
CV and gives the low standards their proper influence. ``1 / y``
(Poisson-like) and unweighted fits are also available. Weights use the *observed*
response, and specifically the *gross* OD (net OD + blank). The constant-CV
noise acts on the total light the reader measures, background included.
Weighting by the blank-subtracted OD alone would treat the lowest standard,
whose net signal might be 0.05, as far more precise than it is, and ordinary
noise there would look like an outlier. Keeping the weights fixed means 4PL
and 5PL minimise the same objective, so their AICc values are comparable.
"""

from __future__ import annotations

import warnings
from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Literal

import numpy as np
import pandas as pd
from scipy.optimize import OptimizeWarning, curve_fit

ModelName = Literal["4pl", "5pl"]
Weighting = Literal["1/y^2", "1/y", "none"]

PARAM_NAMES: dict[str, tuple[str, ...]] = {
    "4pl": ("a", "b", "c", "d"),
    "5pl": ("a", "b", "c", "d", "g"),
}


# --------------------------------------------------------------------------- #
# Model functions and their inverses
# --------------------------------------------------------------------------- #
def four_pl(x: np.ndarray | float, a: float, b: float, c: float, d: float) -> np.ndarray:
    """4PL response at concentration ``x`` (``x >= 0``)."""
    x = np.asarray(x, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        return d + (a - d) / (1.0 + np.power(np.clip(x, 0, None) / c, b))


def five_pl(x: np.ndarray | float, a: float, b: float, c: float, d: float, g: float) -> np.ndarray:
    """5PL response at concentration ``x``. Reduces to :func:`four_pl` when ``g == 1``."""
    x = np.asarray(x, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        return d + (a - d) / np.power(1.0 + np.power(np.clip(x, 0, None) / c, b), g)


def inverse_four_pl(y: np.ndarray | float, a: float, b: float, c: float, d: float) -> np.ndarray:
    """Concentration giving response ``y``. NaN when ``y`` is not strictly between a and d."""
    return inverse_five_pl(y, a, b, c, d, 1.0)


def inverse_five_pl(
    y: np.ndarray | float, a: float, b: float, c: float, d: float, g: float
) -> np.ndarray:
    """Concentration giving response ``y`` under the 5PL.

    Solving the 5PL for ``x`` gives ``x = c * (((a - d)/(y - d))**(1/g) - 1)**(1/b)``.
    That is only real when ``y`` lies strictly between the asymptotes. Outside
    that band no finite concentration produces the signal, so NaN is returned
    and :meth:`CurveFit.back_calculate` reports which side was crossed.
    """
    y = np.asarray(y, dtype=float)
    lo, hi = min(a, d), max(a, d)
    inside = (y > lo) & (y < hi)
    out = np.full(y.shape, np.nan)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        ratio = (a - d) / (y[inside] - d)
        out[inside] = c * np.power(np.power(ratio, 1.0 / g) - 1.0, 1.0 / b)
    return out


_MODEL_FUNCS: dict[str, Callable[..., np.ndarray]] = {"4pl": four_pl, "5pl": five_pl}
_INVERSE_FUNCS: dict[str, Callable[..., np.ndarray]] = {
    "4pl": inverse_four_pl,
    "5pl": inverse_five_pl,
}


# --------------------------------------------------------------------------- #
# Fit result
# --------------------------------------------------------------------------- #
@dataclass
class CurveFit:
    """A fitted standard curve plus its diagnostics.

    Attributes:
        model: ``"4pl"`` or ``"5pl"``.
        params: Fitted parameters by name (a, b, c, d[, g]).
        stderr: Asymptotic standard errors from the scaled covariance matrix.
        r_squared: Unweighted coefficient of determination on the OD scale.
            It is reported because reviewers expect it. It says little about
            low-end accuracy, which is why per-standard % recovery exists.
        aic, aicc: Information criteria computed from the weighted residual
            sum of squares.
        conc, response, weights: Data the curve was fitted to.
        weighting: Weighting scheme used.
        comparison: When chosen via ``model="auto"``, the AICc of every
            candidate model.
    """

    model: ModelName
    params: dict[str, float]
    stderr: dict[str, float]
    r_squared: float
    aic: float
    aicc: float
    conc: np.ndarray
    response: np.ndarray
    weights: np.ndarray
    weighting: str
    comparison: dict[str, float] = field(default_factory=dict)
    weight_offset: float = 0.0

    @property
    def n(self) -> int:
        return len(self.conc)

    @property
    def lower_asymptote(self) -> float:
        return min(self.params["a"], self.params["d"])

    @property
    def upper_asymptote(self) -> float:
        return max(self.params["a"], self.params["d"])

    @property
    def increasing(self) -> bool:
        """True for sandwich-style curves (signal rises with concentration)."""
        return self.params["d"] > self.params["a"]

    def predict(self, x: np.ndarray | float) -> np.ndarray:
        """Response predicted by the fitted curve at concentration(s) ``x``."""
        return _MODEL_FUNCS[self.model](x, *self.params.values())

    def residuals(self) -> np.ndarray:
        """Observed minus fitted response at each standard."""
        return self.response - self.predict(self.conc)

    def weighted_residuals(self) -> np.ndarray:
        """Residuals multiplied by sqrt(weight). They should look like homoscedastic noise."""
        return self.residuals() * np.sqrt(self.weights)

    def back_calculate(self, y: np.ndarray | float) -> tuple[np.ndarray, np.ndarray]:
        """Convert responses into concentrations using the inverse curve.

        Returns:
            ``(conc, flag)``. ``flag`` is ``"ok"``, ``"below_curve"`` (at or
            beyond the zero-concentration asymptote), ``"above_curve"`` (at or
            beyond saturation), or ``"no_signal"`` (input NaN). Out-of-curve
            responses get NaN rather than a clipped value. An asymptote is
            only approached, never reached, so any finite number returned
            would be fabricated.
        """
        y = np.atleast_1d(np.asarray(y, dtype=float))
        conc = _INVERSE_FUNCS[self.model](y, *self.params.values())
        flag = np.full(y.shape, "ok", dtype=object)
        zero_side = self.params["a"]
        sat_side = self.params["d"]
        if self.increasing:
            flag[y <= zero_side] = "below_curve"
            flag[y >= sat_side] = "above_curve"
        else:
            flag[y >= zero_side] = "below_curve"
            flag[y <= sat_side] = "above_curve"
        flag[np.isnan(y)] = "no_signal"
        conc[flag != "ok"] = np.nan
        return conc, flag

    def summary(self) -> dict[str, float | str | int]:
        """Flat dictionary of parameters, SEs, and fit statistics (for tables)."""
        out: dict[str, float | str | int] = {"model": self.model, "n_points": self.n}
        for k, v in self.params.items():
            out[k] = v
            out[f"{k}_se"] = self.stderr.get(k, np.nan)
        out.update(r_squared=self.r_squared, aic=self.aic, aicc=self.aicc)
        return out


# --------------------------------------------------------------------------- #
# Fitting
# --------------------------------------------------------------------------- #
def compute_weights(
    response: np.ndarray, weighting: Weighting | str, offset: float = 0.0
) -> np.ndarray:
    """Weights for weighted least squares, normalised to mean 1.

    ``offset`` is added to the response before weighting. Pass the blank OD
    so weights follow the gross signal (see module docstring). A floor of 2%
    of the response range (at least 1e-3 OD) still keeps any near-zero value
    from getting almost infinite weight.
    """
    y = np.abs(np.asarray(response, dtype=float) + offset)
    floor = max(0.02 * (np.nanmax(y) - np.nanmin(y)), 1e-3)
    y = np.maximum(y, floor)
    if weighting == "1/y^2":
        w = 1.0 / y**2
    elif weighting == "1/y":
        w = 1.0 / y
    elif weighting == "none":
        w = np.ones_like(y)
    else:
        raise ValueError(f"unknown weighting {weighting!r}; use '1/y^2', '1/y' or 'none'")
    return w / w.mean()


def _initial_guess(conc: np.ndarray, y: np.ndarray) -> tuple[list[float], tuple[list, list]]:
    """Data-driven starting values and bounds for the 4PL.

    * ``a`` and ``d``: responses at the lowest and highest standard. The top
      is nudged past the maximum because the highest standard is rarely at
      full saturation.
    * ``c``: the concentration whose response is closest to the midpoint,
      interpolated on the log scale.
    * ``b``: 1, a typical Hill slope for antibody binding.

    The bounds keep the optimiser away from degenerate solutions. For
    example, a slope of 40 fits a step function through two points. The
    bounds are wide enough never to bind on a sensible curve.
    """
    order = np.argsort(conc)
    xs, ys = conc[order], y[order]
    increasing = np.median(ys[-2:]) >= np.median(ys[:2])
    span = float(np.nanmax(y) - np.nanmin(y)) or 1.0
    lo_y, hi_y = float(np.nanmin(y)), float(np.nanmax(y))
    mid = (lo_y + hi_y) / 2
    level_means = pd.Series(ys).groupby(xs).mean()
    idx = int(np.argmin(np.abs(level_means.to_numpy() - mid)))
    c0 = float(level_means.index[idx])
    if increasing:
        a0, d0 = lo_y, hi_y + 0.05 * span
        a_bounds = (lo_y - span, hi_y)
        d_bounds = (lo_y, hi_y + 3 * span)
    else:
        a0, d0 = hi_y + 0.05 * span, lo_y
        a_bounds = (lo_y, hi_y + 3 * span)
        d_bounds = (lo_y - span, hi_y)
    xmin, xmax = float(xs[xs > 0].min()), float(xs.max())
    p0 = [a0, 1.0, c0, d0]
    bounds = (
        [a_bounds[0], 0.1, xmin / 1000, d_bounds[0]],
        [a_bounds[1], 10.0, xmax * 1000, d_bounds[1]],
    )
    return p0, bounds


def _information_criteria(wrss: float, n: int, k: int) -> tuple[float, float]:
    """Gaussian AIC and AICc from a weighted residual sum of squares.

    ``AIC = n·ln(WRSS/n) + 2k``. AICc adds ``2k(k+1)/(n-k-1)``, a correction
    that matters here. With 14 standard wells and 5 parameters, plain AIC
    under-penalises the 5PL.
    """
    wrss = max(wrss, 1e-300)
    aic = n * np.log(wrss / n) + 2 * k
    aicc = aic + (2 * k * (k + 1) / (n - k - 1) if n - k - 1 > 0 else np.inf)
    return float(aic), float(aicc)


def fit_curve(
    conc: np.ndarray,
    response: np.ndarray,
    *,
    model: ModelName | Literal["auto"] = "auto",
    weighting: Weighting = "1/y^2",
    aicc_margin: float = 2.0,
    weight_offset: float = 0.0,
) -> CurveFit:
    """Fit a 4PL or 5PL standard curve by weighted least squares.

    Args:
        conc: Nominal standard concentrations, one per well. Replicate wells
            are fitted individually, not averaged. That keeps the information
            about replicate scatter in the residuals. Zero-concentration wells
            are allowed: the 4PL is defined at x = 0 (y = a), so a zero
            standard anchors the lower asymptote to data instead of leaving
            it to extrapolation.
        response: Blank-subtracted OD for each well.
        model: ``"4pl"``, ``"5pl"``, or ``"auto"``. ``"auto"`` fits both and
            keeps the 5PL only if its AICc beats the 4PL by more than
            ``aicc_margin``. A ΔAICc under 2 means the models are
            indistinguishable, so the simpler one is kept.
        weighting: ``"1/y^2"`` (default, constant-CV noise), ``"1/y"``, or ``"none"``.
        weight_offset: Added to the response when computing weights only,
            normally the plate's blank OD, so weights track the gross signal.

    Raises:
        ValueError: Fewer than 5 usable points or 4 distinct levels (4PL has
            4 parameters, so fewer points leave no degrees of freedom).
        RuntimeError: The optimiser failed to converge.
    """
    conc = np.asarray(conc, dtype=float)
    response = np.asarray(response, dtype=float)
    keep = np.isfinite(conc) & np.isfinite(response) & (conc >= 0)
    conc, response = conc[keep], response[keep]
    n_levels = len(np.unique(conc[conc > 0]))
    if len(conc) < 5 or n_levels < 4:
        raise ValueError(
            f"need >= 5 standard wells over >= 4 non-zero levels to fit a curve "
            f"(got {len(conc)} wells, {n_levels} levels)"
        )

    weights = compute_weights(response, weighting, weight_offset)
    sigma = 1.0 / np.sqrt(weights)
    p0, bounds = _initial_guess(conc, response)

    candidates: dict[str, CurveFit] = {}
    models: list[ModelName] = ["4pl", "5pl"] if model == "auto" else [model]
    for name in models:
        if name == "5pl" and len(conc) < 7:
            continue  # 5 params + AICc needs n - k - 1 > 0 with some slack
        if name == "4pl":
            start, bnds = p0, bounds
        else:
            base = (
                candidates["4pl"].params
                if "4pl" in candidates
                else dict(zip("abcd", p0, strict=True))
            )
            start = [base["a"], base["b"], base["c"], base["d"], 1.0]
            bnds = ([*bounds[0], 0.1], [*bounds[1], 10.0])
            start = list(np.clip(start, bnds[0], bnds[1]))
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", OptimizeWarning)
                popt, pcov = curve_fit(
                    _MODEL_FUNCS[name],
                    conc,
                    response,
                    p0=start,
                    sigma=sigma,
                    absolute_sigma=False,
                    bounds=bnds,
                    maxfev=20000,
                )
        except RuntimeError as exc:
            if model != "auto" or name == "4pl":
                raise RuntimeError(f"{name} fit did not converge: {exc}") from exc
            continue
        names = PARAM_NAMES[name]
        fitted = _MODEL_FUNCS[name](conc, *popt)
        resid = response - fitted
        wrss = float(np.sum(weights * resid**2))
        ss_tot = float(np.sum((response - response.mean()) ** 2))
        r2 = 1.0 - float(np.sum(resid**2)) / ss_tot if ss_tot > 0 else np.nan
        aic, aicc = _information_criteria(wrss, len(conc), len(names))
        with np.errstate(invalid="ignore"):
            se = np.sqrt(np.diag(pcov))
        candidates[name] = CurveFit(
            model=name,
            params=dict(zip(names, map(float, popt), strict=True)),
            stderr=dict(zip(names, map(float, se), strict=True)),
            r_squared=r2,
            aic=aic,
            aicc=aicc,
            conc=conc,
            response=response,
            weights=weights,
            weighting=weighting,
            weight_offset=weight_offset,
        )

    if not candidates:
        raise RuntimeError("no curve model converged")
    comparison = {k: v.aicc for k, v in candidates.items()}
    if model == "auto" and "5pl" in candidates:
        use_5pl = candidates["5pl"].aicc < candidates["4pl"].aicc - aicc_margin
        chosen = candidates["5pl" if use_5pl else "4pl"]
    else:
        chosen = next(iter(candidates.values()))
    chosen.comparison = comparison
    return chosen


# --------------------------------------------------------------------------- #
# Plate-level helpers operating on the tidy well frame
# --------------------------------------------------------------------------- #
def curve_key(row: pd.Series | dict) -> tuple[str, str]:
    """Curves are fitted per (plate_id, analyte)."""
    return (str(row["plate_id"]), str(row["analyte"]))


def fit_standard_curves(
    wells: pd.DataFrame,
    *,
    model: ModelName | Literal["auto"] = "auto",
    weighting: Weighting = "1/y^2",
    include_zero: bool = True,
) -> dict[tuple[str, str], CurveFit]:
    """Fit one curve per (plate_id, analyte) from the standard wells.

    Standard wells are used when ``od_status == "ok"`` and they are not
    flagged as outliers (column ``outlier``, if present). Overflow standards
    are dropped. Their true OD is unknown, and substituting the reader's
    maximum would bend the top of the curve.

    With ``include_zero`` (default), blank wells enter the fit as the zero
    standard (concentration 0), as in SoftMax Pro and most kit protocols.
    Without them the lower asymptote ``a`` is pure extrapolation from the
    lowest standard. On real plates it can land well above the measured
    zero, which then makes every low sample "below the curve".
    """
    types = ["standard", "blank"] if include_zero else ["standard"]
    std = wells[wells["type"].isin(types) & (wells["od_status"] == "ok")].copy()
    std.loc[std["type"] == "blank", "concentration"] = 0.0
    if "outlier" in std.columns:
        std = std[~std["outlier"].astype(bool)]
    fits: dict[tuple[str, str], CurveFit] = {}
    for (pid, analyte), grp in std.groupby(["plate_id", "analyte"], sort=False):
        try:
            fits[(str(pid), str(analyte))] = fit_curve(
                grp["concentration"].to_numpy(),
                grp["od_net"].to_numpy(),
                model=model,
                weighting=weighting,
                weight_offset=float(grp["blank_mean"].iloc[0]),
            )
        except (ValueError, RuntimeError) as exc:
            raise RuntimeError(f"plate {pid!r} / {analyte}: {exc}") from exc
    return fits


def back_calculate(wells: pd.DataFrame, fits: dict[tuple[str, str], CurveFit]) -> pd.DataFrame:
    """Add ``conc_well`` (in-well concentration, before dilution) and ``curve_flag``.

    Overflow wells are flagged ``above_curve``. A saturated reading is by
    definition above the top of the measurable range.
    """
    out = wells.copy()
    out["conc_well"] = np.nan
    out["curve_flag"] = "no_curve"
    for key, grp in out.groupby(["plate_id", "analyte"], sort=False):
        fit = fits.get((str(key[0]), str(key[1])))
        if fit is None:
            continue
        conc, flag = fit.back_calculate(grp["od_net"].to_numpy())
        flag = np.where(grp["od_status"].to_numpy() == "overflow", "above_curve", flag)
        out.loc[grp.index, "conc_well"] = conc
        out.loc[grp.index, "curve_flag"] = flag
    return out


def standard_recovery(wells: pd.DataFrame) -> pd.DataFrame:
    """Per-standard-level back-calculated mean, %CV, and % recovery.

    ``% recovery = mean back-calculated concentration / nominal × 100``.
    This is the honest test of a curve. R² can be 0.999 while the lowest
    standards back-calculate at 60% of nominal, because the low end
    contributes almost nothing to the unweighted sum of squares.

    Expects the output of :func:`back_calculate`.
    """
    std = wells[wells["type"] == "standard"].copy()
    if "outlier" in std.columns:
        std = std[~std["outlier"].astype(bool)]
    rows = []
    for (pid, analyte, level), grp in std.groupby(
        ["plate_id", "analyte", "concentration"], sort=True
    ):
        vals = grp["conc_well"].to_numpy(dtype=float)
        finite = vals[np.isfinite(vals)]
        mean = float(finite.mean()) if finite.size else np.nan
        sd = float(finite.std(ddof=1)) if finite.size > 1 else np.nan
        rows.append(
            {
                "plate_id": pid,
                "analyte": analyte,
                "concentration": float(level),
                "n_wells": len(grp),
                "n_quantified": int(finite.size),
                "mean_back_calc": mean,
                "cv_pct": 100 * sd / mean if finite.size > 1 and mean > 0 else np.nan,
                "recovery_pct": 100 * mean / float(level) if finite.size else np.nan,
            }
        )
    return pd.DataFrame(rows)
