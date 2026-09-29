"""Static figures for reports and the README (matplotlib).

The colors are one validated categorical palette, assigned to groups in a
fixed order so a group keeps its color in every figure. Text and axes use
neutral inks, never series colors. Identity is never carried by color alone:
groups are also named on the axis, and excluded wells use a distinct marker
shape as well as a color.
"""

from __future__ import annotations

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.figure import Figure

from immunoassay.curves import CurveFit

SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
SURFACE = "#fcfcfb"
CRITICAL = "#d03b3b"
BAND = "#f0efec"
SEQ_BLUE = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]

_RC = {
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 9,
    "text.color": INK,
    "axes.labelcolor": INK_2,
    "axes.edgecolor": AXIS,
    "axes.titlesize": 10,
    "axes.titleweight": "bold",
    "axes.titlecolor": INK,
    "axes.titlelocation": "left",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.6,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "xtick.labelcolor": INK_2,
    "ytick.labelcolor": INK_2,
    "legend.frameon": False,
    "legend.fontsize": 8,
    "lines.linewidth": 2,
}


def style() -> matplotlib.RcParams:
    """Context manager applying the report style: ``with style(): ...``."""
    return plt.rc_context(_RC)  # type: ignore[return-value]


def group_colors(groups: list[str]) -> dict[str, str]:
    """Fixed-order color assignment. A group's color never depends on its rank."""
    return {g: SERIES[i % len(SERIES)] for i, g in enumerate(groups)}


def plot_standard_curve(
    fit: CurveFit,
    plate_wells: pd.DataFrame,
    limits: pd.Series | None = None,
    *,
    title: str = "",
    units: str = "",
) -> Figure:
    """Standard curve (log x) with fitted model, quantification range, and weighted residuals."""
    with style():
        fig, (ax, axr) = plt.subplots(
            2, 1, figsize=(5.2, 4.6), sharex=True, gridspec_kw={"height_ratios": [3, 1.2]}
        )
        std = plate_wells[plate_wells["type"] == "standard"]
        used = std[~std.get("outlier", False).astype(bool) & (std["od_status"] == "ok")]
        excl = std[std.get("outlier", False).astype(bool)]
        xs = np.geomspace(std["concentration"].min() / 3, std["concentration"].max() * 3, 300)

        if limits is not None and np.isfinite(limits.get("lloq", np.nan)):
            ax.axvspan(limits["lloq"], limits["uloq"], color=BAND, zorder=0, lw=0)
            ax.text(
                np.sqrt(limits["lloq"] * limits["uloq"]),
                0.97,
                "quantifiable range",
                transform=ax.get_xaxis_transform(),
                ha="center",
                va="top",
                fontsize=7.5,
                color=INK_2,
            )
        if limits is not None and np.isfinite(limits.get("lod", np.nan)):
            ax.axvline(limits["lod"], color=MUTED, lw=1, ls=(0, (3, 3)))
            ax.text(
                limits["lod"],
                0.03,
                " LOD",
                transform=ax.get_xaxis_transform(),
                fontsize=7.5,
                color=INK_2,
                ha="left",
            )
        ax.plot(
            xs,
            fit.predict(xs),
            color=INK_2,
            lw=1.5,
            zorder=2,
            label=f"{fit.model.upper()} fit (R² {fit.r_squared:.4f})",
        )
        ax.scatter(
            used["concentration"],
            used["od_net"],
            s=30,
            color=SERIES[0],
            edgecolor=SURFACE,
            linewidth=1.2,
            zorder=3,
            label="standard wells",
        )
        if len(excl):
            ax.scatter(
                excl["concentration"],
                excl["od_net"],
                s=40,
                marker="x",
                color=CRITICAL,
                linewidth=1.6,
                zorder=4,
                label="excluded (outlier)",
            )
        ax.set_xscale("log")
        ax.set_ylabel("net OD")
        ax.set_title(title or "Standard curve")
        ax.legend(loc="upper left")

        resid = fit.weighted_residuals()
        axr.axhline(0, color=AXIS, lw=1)
        axr.scatter(fit.conc, resid, s=18, color=SERIES[0], edgecolor=SURFACE, linewidth=1)
        lim = max(np.abs(resid).max() * 1.3, 1e-3)
        axr.set_ylim(-lim, lim)
        axr.set_ylabel("weighted\nresidual")
        axr.set_xlabel(f"concentration{f' ({units})' if units else ''}")
        fig.tight_layout()
    return fig


def plot_plate_heatmap(plate_wells: pd.DataFrame, *, title: str = "") -> Figure:
    """Raw OD across the 8 × 12 plate. Edge effects and pipetting gradients show up here."""
    from immunoassay.io import COLUMNS, ROWS

    grid = plate_wells.pivot_table(index="row", columns="column", values="od").reindex(
        index=list(ROWS), columns=list(COLUMNS)
    )
    with style():
        fig, ax = plt.subplots(figsize=(5.6, 3.2))
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list("seq", SEQ_BLUE)
        im = ax.imshow(grid.to_numpy(dtype=float), cmap=cmap, aspect="auto")
        ax.set_xticks(range(12), [str(c) for c in COLUMNS])
        ax.set_yticks(range(8), list(ROWS))
        ax.grid(False)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.tick_params(length=0)
        out = plate_wells[plate_wells.get("outlier", False).astype(bool)]
        for _, w in out.iterrows():
            ax.scatter(
                int(w["column"]) - 1,
                ROWS.index(w["row"]),
                marker="x",
                s=40,
                color=CRITICAL,
                linewidth=1.6,
            )
        cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
        cb.outline.set_visible(False)
        cb.set_label("raw OD", color=INK_2)
        ax.set_title(title or "Plate map")
        fig.tight_layout()
    return fig


def plot_normalization(results: pd.DataFrame, analyte: str, *, units: str = "") -> Figure:
    """Bridge controls across plates before vs after inter-plate normalization.

    Each line is one control sample. Before correction the lines are
    tilted, because every plate reads all controls high or low together.
    After correction they should be flat.
    """
    ctrl = results[(results["type"] == "control") & (results["analyte"] == analyte)]
    plates = list(dict.fromkeys(ctrl["plate_id"]))
    ids = sorted(ctrl["sample_id"].unique())
    colors = {s: SERIES[i % len(SERIES)] for i, s in enumerate(ids)}
    with style():
        fig, axes = plt.subplots(1, 2, figsize=(6.4, 3.0), sharey=True)
        for ax, col, label in zip(axes, ["conc", "conc_norm"], ["Before", "After"], strict=True):
            for sid in ids:
                d = ctrl[ctrl["sample_id"] == sid].set_index("plate_id").reindex(plates)
                x = np.arange(len(plates))
                y = d[col].to_numpy(dtype=float)
                ax.plot(
                    x,
                    y,
                    color=colors[sid],
                    lw=2,
                    marker="o",
                    ms=6,
                    markeredgecolor=SURFACE,
                    markeredgewidth=1.2,
                )
                last = np.where(np.isfinite(y))[0]
                if len(last) and ax is axes[1]:
                    ax.text(
                        x[last[-1]] + 0.08, y[last[-1]], sid, va="center", fontsize=7.5, color=INK_2
                    )
            ax.set_xticks(range(len(plates)), plates)
            ax.set_xlim(-0.3, len(plates) - 0.4 + (0.6 if ax is axes[1] else 0))
            ax.set_yscale("log")
            ax.set_title(f"{label} normalization")
        axes[0].set_ylabel(f"{analyte} ({units})" if units else analyte)
        fig.tight_layout()
    return fig


def plot_group_comparison(
    data: pd.DataFrame,
    analyte: str,
    *,
    group_order: list[str] | None = None,
    pairwise: pd.DataFrame | None = None,
    units: str = "",
    log: bool = True,
) -> Figure:
    """Box plot plus every individual sample. Substituted (censored) values are hollow."""
    d = data[(data["analyte"] == analyte) & np.isfinite(data["value"])]
    groups = group_order or list(dict.fromkeys(d["group"]))
    groups = [g for g in groups if g in set(d["group"])]
    colors = group_colors(group_order or groups)
    rng = np.random.default_rng(0)
    with style():
        fig, ax = plt.subplots(figsize=(1.3 + 1.25 * len(groups), 3.4))
        for i, g in enumerate(groups):
            v = d[d["group"] == g]
            ax.boxplot(
                v["value"],
                positions=[i],
                widths=0.5,
                showfliers=False,
                patch_artist=True,
                boxprops={"facecolor": BAND, "edgecolor": AXIS, "linewidth": 1},
                medianprops={"color": INK, "linewidth": 1.5},
                whiskerprops={"color": AXIS},
                capprops={"color": AXIS},
            )
            jitter = rng.uniform(-0.13, 0.13, len(v))
            measured = v["handling"].eq("measured").to_numpy()
            ax.scatter(
                i + jitter[measured],
                v["value"][measured],
                s=26,
                color=colors[g],
                edgecolor=SURFACE,
                linewidth=1,
                zorder=3,
            )
            if (~measured).any():
                ax.scatter(
                    i + jitter[~measured],
                    v["value"][~measured],
                    s=26,
                    facecolor=SURFACE,
                    edgecolor=colors[g],
                    linewidth=1.4,
                    zorder=3,
                )
        ax.set_xticks(range(len(groups)), groups)
        ax.grid(axis="x", visible=False)
        if log:
            ax.set_yscale("log")
        ax.set_ylabel(f"{analyte} ({units})" if units else analyte)
        ax.set_title(analyte)
        if pairwise is not None and not pairwise.empty:
            pw = pairwise[pairwise["analyte"] == analyte]
            ymax = d["value"].max()
            step = 1.35 if log else 0.08 * (d["value"].max() - d["value"].min())
            level = 0
            for r in pw.itertuples():
                if r.reference not in groups or r.group not in groups:
                    continue
                i, j = groups.index(r.reference), groups.index(r.group)
                y = ymax * step ** (level + 1) if log else ymax + step * (level + 1)
                ax.plot(
                    [i, i, j, j],
                    [y / 1.05, y, y, y / 1.05] if log else [y - step / 4, y, y, y - step / 4],
                    color=MUTED,
                    lw=1,
                )
                stars = (
                    "***"
                    if r.p_holm < 0.001
                    else "**"
                    if r.p_holm < 0.01
                    else "*"
                    if r.p_holm < 0.05
                    else "ns"
                )
                ax.text((i + j) / 2, y, stars, ha="center", va="bottom", fontsize=8, color=INK_2)
                level += 1
        fig.tight_layout()
    return fig


def fig_to_png(fig: Figure, path=None, dpi: int = 150) -> bytes:
    """Render a figure to PNG bytes (and optionally save it). Closes the figure."""
    import io

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    data = buf.getvalue()
    if path is not None:
        with open(path, "wb") as fh:
            fh.write(data)
    return data
