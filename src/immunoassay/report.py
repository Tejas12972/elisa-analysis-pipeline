"""Self-contained HTML report (figures embedded as base64 PNG, no external assets).

The report follows the order a reviewer audits an ELISA run in. First, did
the plates pass? Then, are the curves sound? Which samples are flagged, and
why? Only then come the biological results. The figures are static PNGs, so
the file can be emailed, archived with a lab notebook, or printed to PDF
from any browser.
"""

from __future__ import annotations

import base64
import html
import json
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

from immunoassay import __version__
from immunoassay import plotting as plots
from immunoassay.pipeline import PipelineResult

_CSS = """
:root{--surface:#fcfcfb;--page:#f9f9f7;--ink:#0b0b0b;--ink2:#52514e;--muted:#898781;
--grid:#e1e0d9;--good:#006300;--bad:#d03b3b;--band:#f0efec;--accent:#2a78d6}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){--surface:#1a1a19;
--page:#0d0d0d;--ink:#fff;--ink2:#c3c2b7;--grid:#2c2c2a;--good:#0ca30c;--band:#383835;
--accent:#3987e5}}
:root[data-theme="dark"]{--surface:#1a1a19;--page:#0d0d0d;--ink:#fff;--ink2:#c3c2b7;
--grid:#2c2c2a;--good:#0ca30c;--band:#383835;--accent:#3987e5}
*{box-sizing:border-box}
body{margin:0;background:var(--page);color:var(--ink);
font:14px/1.5 system-ui,-apple-system,"Segoe UI",sans-serif}
main{max-width:1100px;margin:0 auto;padding:24px 16px 64px}
h1{font-size:24px;margin:0 0 4px}h2{font-size:18px;margin:40px 0 8px;
border-bottom:1px solid var(--grid);padding-bottom:6px}h3{font-size:15px;margin:24px 0 8px}
.sub{color:var(--ink2);margin:0 0 20px}
.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:12px}
.tile{background:var(--surface);border:1px solid var(--grid);border-radius:8px;padding:12px 14px}
.tile .v{font-size:24px;font-weight:600}.tile .l{color:var(--ink2);font-size:12px}
.wrap{overflow-x:auto;background:var(--surface);border:1px solid var(--grid);border-radius:8px}
table{border-collapse:collapse;width:100%;font-size:12.5px;font-variant-numeric:tabular-nums}
th,td{padding:6px 10px;text-align:left;border-bottom:1px solid var(--grid);white-space:nowrap}
th{color:var(--ink2);font-weight:600;background:var(--band)}
tr:last-child td{border-bottom:0}
.pass{color:var(--good);font-weight:600}.fail{color:var(--bad);font-weight:600}
.figs{display:flex;flex-wrap:wrap;gap:12px}
.figs img{max-width:100%;height:auto;background:#fcfcfb;border:1px solid var(--grid);
border-radius:8px}
.note{color:var(--ink2);font-size:13px;max-width:820px}
details{margin:8px 0}summary{cursor:pointer;color:var(--accent)}
pre{background:var(--surface);border:1px solid var(--grid);border-radius:8px;padding:12px;
overflow-x:auto;font-size:12px}
"""


def _img(fig) -> str:
    data = base64.b64encode(plots.fig_to_png(fig)).decode()
    return f'<img alt="figure" src="data:image/png;base64,{data}">'


def _fmt(v) -> str:
    if isinstance(v, (bool, np.bool_)):
        return '<span class="pass">✓ pass</span>' if v else '<span class="fail">✗ fail</span>'
    if isinstance(v, (float, np.floating)):
        if not np.isfinite(v):
            return "–"
        if v != 0 and (abs(v) < 0.01 or abs(v) >= 1e5):
            return f"{v:.3g}"
        return f"{v:.4g}" if abs(v) < 10 else f"{v:,.1f}"
    if v is None or (isinstance(v, str) and v == "") or v is pd.NA:
        return ""
    return html.escape(str(v))


def _table(df: pd.DataFrame, columns: list[str] | None = None, rename: dict | None = None) -> str:
    if df is None or df.empty:
        return '<p class="note">None.</p>'
    cols = [c for c in (columns or list(df.columns)) if c in df.columns]
    head = "".join(f"<th>{html.escape((rename or {}).get(c, c))}</th>" for c in cols)
    body = "".join(
        "<tr>" + "".join(f"<td>{_fmt(v)}</td>" for v in row) + "</tr>"
        for row in df[cols].itertuples(index=False)
    )
    table = f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"
    return f'<div class="wrap">{table}</div>'


def build_report(
    result: PipelineResult,
    *,
    title: str = "ELISA analysis report",
    subtitle: str = "",
    units: dict[str, str] | None = None,
) -> str:
    """Render the full HTML report for a :class:`PipelineResult`."""
    units = units or result.config.units or {}
    cfg = result.config
    qc = result.plate_qc
    samples = result.results[result.results["type"] == "sample"]
    flagged = result.samples[
        (result.samples["flags"] != "") & (result.samples["type"].isin(["sample", "control"]))
    ]
    used = result.stats.data["handling"].eq("measured").sum() if len(result.stats.data) else 0
    parts: list[str] = []
    add = parts.append

    add(f"<h1>{html.escape(title)}</h1>")
    add(
        f'<p class="sub">{html.escape(subtitle)} · generated {datetime.now():%Y-%m-%d %H:%M} · '
        f"immunoassay v{__version__}</p>"
    )
    add('<div class="tiles">')
    for value, label in [
        (f"{int(qc['passed'].sum())}/{len(qc)}", "plates passed QC"),
        (", ".join(sorted(qc["analyte"].unique())), "analytes"),
        (str(samples["sample_id"].nunique()), "samples"),
        (str(len(flagged)), "flagged replicate groups"),
        (str(used), "values measured in range"),
    ]:
        add(
            f'<div class="tile"><div class="v">{html.escape(value)}</div>'
            f'<div class="l">{label}</div></div>'
        )
    add("</div>")

    # 1. Plate QC
    add("<h2>1 · Plate acceptance</h2>")
    add(
        f'<p class="note">A plate passes when R² ≥ {cfg.qc.min_r_squared:.2f}, ≥ '
        f"{cfg.qc.min_std_pass_frac:.0%} of standard levels recover "
        f"within {cfg.qc.recovery_low_pct:.0f}–{cfg.qc.recovery_high_pct:.0f}% with CV ≤ "
        f"{cfg.qc.max_cv_pct:.0f}%, a quantifiable range exists, and ≥ "
        f"{cfg.qc.min_control_pass_frac:.0%} of "
        "controls with expected ranges fall inside them. Limits are in-well "
        "concentrations.</p>"
    )
    add(
        _table(
            qc,
            [
                "plate_id",
                "analyte",
                "passed",
                "model",
                "r_squared",
                "std_levels_pass",
                "lod",
                "lloq",
                "uloq",
                "controls_pass",
                "n_flagged",
                "reasons",
            ],
            {
                "passed": "verdict",
                "r_squared": "R²",
                "std_levels_pass": "standards OK",
                "controls_pass": "controls OK",
                "n_flagged": "flagged groups",
            },
        )
    )

    # 2. Curves
    add("<h2>2 · Standard curves</h2>")
    add(
        f'<p class="note">Weighted ({html.escape(cfg.weighting)}) least-squares fits. 5PL is '
        f"chosen over 4PL only when "
        "its AICc is lower by more than 2. Weighted residuals should scatter evenly around "
        "zero; a curved pattern means the model is wrong for the data. The plate map shows "
        "raw OD. Excluded wells are marked ✗.</p>"
    )
    curve_tab = result.curve_table()
    for (pid, analyte), fit in result.fits.items():
        plate_wells = result.wells[
            (result.wells["plate_id"] == pid) & (result.wells["analyte"] == analyte)
        ]
        lim = result.limits.set_index(["plate_id", "analyte"]).loc[(pid, analyte)]
        add(f"<h3>{html.escape(pid)} · {html.escape(analyte)}</h3>")
        add('<div class="figs">')
        add(
            _img(
                plots.plot_standard_curve(
                    fit, plate_wells, lim, title=f"{pid} · {analyte}", units=units.get(analyte, "")
                )
            )
        )
        add(_img(plots.plot_plate_heatmap(plate_wells, title=f"{pid} raw OD")))
        add("</div>")
        row = curve_tab[(curve_tab["plate_id"] == pid) & (curve_tab["analyte"] == analyte)]
        pcols = [
            c
            for c in [
                "model",
                "a",
                "a_se",
                "b",
                "b_se",
                "c",
                "c_se",
                "d",
                "d_se",
                "g",
                "g_se",
                "r_squared",
                "aicc_4pl",
                "aicc_5pl",
            ]
            if c in row
        ]
        add(_table(row, pcols, {"r_squared": "R²", "c": "c (EC50)"}))
        rec = result.recovery[
            (result.recovery["plate_id"] == pid) & (result.recovery["analyte"] == analyte)
        ]
        add("<details><summary>Standard recovery</summary>")
        add(
            _table(
                rec,
                [
                    "concentration",
                    "n_wells",
                    "n_quantified",
                    "mean_back_calc",
                    "cv_pct",
                    "recovery_pct",
                ],
            )
        )
        add("</details>")
    outl = result.wells[result.wells["outlier"].astype(bool)]
    add("<h3>Excluded wells</h3>")
    add(_table(outl, ["plate_id", "well", "type", "sample_id", "od", "outlier_reason"]))

    # 3. Controls + flags
    add("<h2>3 · Controls and flagged samples</h2>")
    add(
        _table(
            result.controls,
            [
                "plate_id",
                "analyte",
                "sample_id",
                "conc_neat",
                "low",
                "high",
                "cv_pct",
                "range_flag",
                "passed",
            ],
            {"conc_neat": "measured", "passed": "in range"},
        )
    )
    add("<h3>Flagged replicate groups</h3>")
    add(
        '<p class="note"><b>below_lloq / above_uloq</b>: outside the validated range '
        "(not reported as numbers). <b>below_lod</b>: not distinguishable from blank. "
        f"<b>cv_high</b>: replicate CV above {cfg.qc.max_cv_pct:.0f}%.</p>"
    )
    add(
        _table(
            flagged,
            [
                "plate_id",
                "analyte",
                "sample_id",
                "group",
                "dilution",
                "wells",
                "od_mean",
                "conc_well",
                "cv_pct",
                "flags",
            ],
        )
    )

    # 4. Dilution
    add("<h2>4 · Dilution linearity</h2>")
    add(
        '<p class="note">Samples run at more than one dilution should agree after correction '
        f"(±{cfg.linearity_max_deviation_pct:.0f}%). A concentration that rises with dilution "
        f"points to matrix interference "
        "or the high-dose hook effect.</p>"
    )
    add(
        _table(
            result.linearity,
            [
                "plate_id",
                "analyte",
                "sample_id",
                "dilutions",
                "neat_concs",
                "max_deviation_pct",
                "trend",
                "linear",
                "interpretation",
            ],
        )
    )

    # 5. Normalization
    add("<h2>5 · Inter-plate normalization</h2>")
    add(
        '<p class="note">Bridge controls shared by every plate estimate a multiplicative '
        "plate factor, on the log scale with the average plate as reference. Corrected "
        "concentration = raw / factor.</p>"
    )
    add(
        _table(
            result.plate_factors,
            ["plate_id", "analyte", "factor", "log_factor_se", "n_bridges", "bridge_ids"],
        )
    )
    add(
        _table(
            result.bridge_cv,
            None,
            {
                "inter_plate_cv_pct_before": "inter-plate CV % before",
                "inter_plate_cv_pct_after": "inter-plate CV % after",
            },
        )
    )
    multi = result.results.groupby("analyte")["plate_id"].nunique()
    add('<div class="figs">')
    for analyte in multi[multi > 1].index:
        add(_img(plots.plot_normalization(result.results, analyte, units=units.get(analyte, ""))))
    add("</div>")

    # 6. Stats
    add("<h2>6 · Group comparisons</h2>")
    handling = result.stats.data.groupby(["analyte", "handling"]).size().reset_index(name="n")
    add(
        '<p class="note">Tests are run on {} concentrations. Shapiro-Wilk on every group '
        "decides between parametric tests (Welch t / Welch ANOVA) and rank tests (Mann-Whitney "
        "/ Kruskal-Wallis). Omnibus p-values are Benjamini-Hochberg adjusted across analytes; "
        "pairwise p-values are Holm-adjusted within each analyte. Hollow points are censored "
        "values that were substituted.</p>".format(
            "log10-transformed" if cfg.stats.transform == "log10" else "raw"
        )
    )
    add(_table(handling, None))
    add("<h3>Omnibus tests</h3>")
    add(
        _table(
            result.stats.omnibus,
            [
                "analyte",
                "groups",
                "n",
                "test",
                "statistic",
                "p",
                "p_adj_bh",
                "significant",
                "omnibus_effect_name",
                "omnibus_effect",
                "shapiro",
            ],
            {"p_adj_bh": "p (BH)", "omnibus_effect_name": "effect", "omnibus_effect": "value"},
        )
    )
    add("<h3>Pairwise comparisons</h3>")
    add(
        _table(
            result.stats.pairwise,
            [
                "analyte",
                "reference",
                "group",
                "test",
                "effect_name",
                "estimate",
                "ci_low",
                "ci_high",
                "std_effect_name",
                "std_effect",
                "p_holm",
                "significant",
            ],
            {
                "estimate": "estimate",
                "ci_low": "95% CI low",
                "ci_high": "95% CI high",
                "std_effect_name": "std. effect",
                "std_effect": "value",
                "p_holm": "p (Holm)",
            },
        )
    )
    add('<div class="figs">')
    for analyte in result.stats.data["analyte"].unique():
        add(
            _img(
                plots.plot_group_comparison(
                    result.stats.data,
                    analyte,
                    group_order=cfg.group_order,
                    pairwise=result.stats.pairwise,
                    units=units.get(analyte, ""),
                    log=cfg.stats.transform == "log10",
                )
            )
        )
    add("</div>")

    # 7. Methods
    add("<h2>7 · Settings and warnings</h2>")
    add(f"<pre>{html.escape(json.dumps(cfg.as_dict(), indent=2, default=str))}</pre>")
    if result.warnings:
        add("<ul>" + "".join(f"<li>{html.escape(w)}</li>" for w in result.warnings) + "</ul>")

    body = "\n".join(parts)
    return (
        '<!doctype html><html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width,initial-scale=1">'
        f"<title>{html.escape(title)}</title><style>{_CSS}</style></head>"
        f"<body><main>{body}</main></body></html>"
    )


def write_report(result: PipelineResult, path: str | Path, **kwargs) -> Path:
    """Write :func:`build_report` output to ``path``."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(build_report(result, **kwargs), encoding="utf-8")
    return path
