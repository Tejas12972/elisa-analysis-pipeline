"""Build and execute the worked-example notebooks in ``notebooks/``.

The notebooks are generated from this file so that their source stays
reviewable as plain Python. Executed outputs are saved, so GitHub renders
them without running anything.

Usage::

    python scripts/build_notebooks.py            # build + execute both
    python scripts/build_notebooks.py --no-run   # build only
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import nbformat
from nbclient import NotebookClient

REPO = Path(__file__).resolve().parents[1]
NB = REPO / "notebooks"


def md(text: str):
    return nbformat.v4.new_markdown_cell(text.strip())


def code(text: str):
    return nbformat.v4.new_code_cell(text.strip())


def setup(third_party: str = "", first_party: str = "") -> str:
    """First code cell: imports, notebook-wide options, and the ``show`` helper."""
    return (
        """
import subprocess
import sys
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from IPython.display import display
"""
        + third_party
        + """
from immunoassay import plotting as P
"""
        + first_party
        + """from immunoassay.io import PlateWarning
from immunoassay.pipeline import PipelineConfig, run_pipeline

warnings.simplefilter("ignore", PlateWarning)  # unused wells on these plates are expected
pd.set_option("display.precision", 3)
pd.set_option("display.max_columns", 30)
REPO = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd()
%matplotlib inline


def show(fig):
    \"\"\"Display a figure returned by immunoassay.plotting once, then free it.\"\"\"
    display(fig)
    plt.close(fig)
"""
    )


# --------------------------------------------------------------------------- #
# Notebook 1: sleep fragmentation IL-6 (real data, biological question)
# --------------------------------------------------------------------------- #
NB1 = [
    md("""
# Serum IL-6 after acute sleep fragmentation: re-analysing a public ELISA dataset

**Data.** Nguyen, Fields & Ashley (Western Kentucky University), *Gene expression and ELISA
data*, Dryad [doi:10.5061/dryad.tdz08kq3p](https://doi.org/10.5061/dryad.tdz08kq3p), **CC0**.
The data accompany *"Inflammation from sleep fragmentation starts in the periphery rather
than brain in male mice"*. The raw workbook is kept unmodified in
`data/raw/sleep_fragmentation_il6/`, and `PROVENANCE.md` there documents every cell used.

**Design.** Male mice were exposed to **acute sleep fragmentation (ASF)** or left undisturbed
(**NSF**, non-sleep-fragmented controls). Trunk-blood serum was collected at 0, 1, 2, 6, 12 and
24 h, with about 10 mice per group, 110 samples in total. Serum IL-6 was measured with a
BioLegend ELISA MAX mouse IL-6 sandwich ELISA on four 96-well plates.

**Biological question.** IL-6 is a pleiotropic cytokine. It is produced by macrophages,
endothelium and adipose tissue, and it is the main driver of the hepatic acute-phase response.
If sleep disruption triggers *systemic* inflammation, circulating IL-6 should rise after ASF.
The time course then says how quickly this happens and how long it lasts.

**Why re-analyse?** The authors calibrated each plate with a hand-typed **linear** fit and
reported every sample as a number. That includes samples that read below the lowest standard,
and some values are set to 0 by hand. This notebook runs the raw ODs through the pipeline:
weighted 4PL/5PL curves, replicate and standard QC, LOD/LLOQ/ULOQ, and explicit handling of
censored values. It then asks which of the original conclusions survive.
"""),
    code(
        setup(
            third_party="from scipy.stats import mannwhitneyu\n"
            "from statsmodels.stats.multitest import multipletests\n",
            first_party="from immunoassay import stats as st\n",
        )
        + """
base = REPO / "data/processed/sleep_fragmentation_il6"
if not (base / "manifest.csv").exists():
    subprocess.run([sys.executable, str(REPO / "scripts/prepare_sleep_fragmentation.py")], check=True)
info = pd.read_csv(base / "sample_info.csv")
authors = pd.read_csv(base / "authors_concentrations.csv")
info.head()
"""
    ),
    md("""
## 1 · The plate layout, and a confound built into it

`scripts/prepare_sleep_fragmentation.py` rebuilds each plate's layout from the workbook's own
formulas (`M_ASF2h-1 → =AVERAGE(G30:H30)` gives wells C5/C6), so no well assignment is typed by
hand. Dilution factors (2× or 3×) come from the `2*(…)` / `3*(…)` multipliers in the same
formulas. Here is where each group was assayed:
"""),
    code("""
info.pivot_table(index="plate", columns="group", values="sample_id", aggfunc="count", fill_value=0)
"""),
    md("""
Each plate holds whole groups, so **plate is confounded with time point**. For 1, 2, 12 and
24 h, the ASF and NSF groups share a plate, so their comparison is within-plate. At **6 h,
ASF is on plate 3 and NSF on plate 2**. Any plate-to-plate difference in calibration therefore
adds directly to the 6 h treatment effect. There are no inter-plate bridge controls, so this
cannot be corrected after the fact (see the demo study in the README for what bridging buys).

The raw plate maps show the layout: standards in columns 1–2 (zero standard in H1/H2), then
duplicate sample pairs.
"""),
    code("""
cfg = PipelineConfig(normalize=False)
result = run_pipeline(base / "manifest.csv", cfg)
for pid in ["plate1", "plate2", "plate3", "plate4"]:
    show(P.plot_plate_heatmap(result.wells[result.wells.plate_id == pid], title=f"{pid} raw OD"))
"""),
    md("""
## 2 · Standard curves: weighted 4PL vs the authors' straight lines

Every plate is fitted with a 4PL and a 5PL by weighted least squares, with weights 1/OD² on the
*gross* signal. The zero standard is included as a calibrator at x = 0. The 5PL is kept only if
it lowers AICc by more than 2. Here it never does, so all four plates use the 4PL.
"""),
    code("""
curves = result.curve_table()
curves[["plate_id", "model", "a", "a_se", "b", "b_se", "c", "c_se", "d", "d_se", "r_squared",
        "aicc_4pl", "aicc_5pl"]]
"""),
    md("""
The authors' linear calibrations, as typed in the workbook, were `conc = (OD − a) / b`:

| plate | a | b | response used |
|---|---|---|---|
| 1 | 0.0032 | 0.0030 | OD − blank (H9/H10) |
| 2 | 0.0097 | 0.0144 | OD − blank (H9/H10) |
| 3 | 0.1167 | 0.0104 | raw OD |
| 4 | 0.1038 | 0.0103 | raw OD |

Both calibrations can be applied to the *standards themselves*. Back-calculating the standards
and comparing them with their nominal values (% recovery) is the direct test of a calibration.
"""),
    code("""
LINEAR = {"plate1": (0.0032, 0.0030, "net"), "plate2": (0.0097, 0.0144, "net"),
          "plate3": (0.1167, 0.0104, "raw"), "plate4": (0.1038, 0.0103, "raw")}
AUTHOR_BLANK = {"plate1": ["H9", "H10"], "plate2": ["H9", "H10"]}
plates_raw = {pid: pd.read_csv(base / f"plates/{pid}.csv", skiprows=3, index_col=0)
              for pid in LINEAR}

def linear_conc(pid, od_raw):
    a, b, kind = LINEAR[pid]
    if kind == "net":
        grid = plates_raw[pid]
        blank = np.mean([grid.loc[w[0], w[1:]] for w in AUTHOR_BLANK[pid]])
        od_raw = od_raw - blank
    return (od_raw - a) / b

std = result.wells[result.wells.type == "standard"].copy()
lvl = (std.groupby(["plate_id", "concentration"])
          .agg(od=("od", "mean")).reset_index())
lvl["linear_recovery_pct"] = [100 * linear_conc(p, o) / c
                              for p, o, c in zip(lvl.plate_id, lvl.od, lvl.concentration, strict=True)]
rec = result.recovery.rename(columns={"recovery_pct": "4pl_recovery_pct"})
comparison = rec.merge(lvl, on=["plate_id", "concentration"])
(comparison.pivot(index="concentration", columns="plate_id",
                  values=["4pl_recovery_pct", "linear_recovery_pct"]).round(0))
"""),
    code("""
with P.style():
    fig, axes = plt.subplots(1, 4, figsize=(15, 3.6), sharey=False)
    for ax, pid in zip(axes, LINEAR, strict=True):
        fit = result.fits[(pid, "IL-6")]
        s = std[(std.plate_id == pid) & ~std.outlier]
        xs = np.geomspace(3, 1500, 200)
        ax.plot(xs, fit.predict(xs), color=P.INK_2, lw=1.6, label="weighted 4PL")
        a, b, kind = LINEAR[pid]
        blank_ours = s["blank_mean"].iloc[0]
        shift = 0.0 if kind == "net" else blank_ours
        ax.plot(xs, a + b * xs - shift, color=P.SERIES[1], lw=1.6, ls="--",
                label="authors' linear")
        ax.scatter(s.concentration, s.od_net, s=26, color=P.SERIES[0], edgecolor=P.SURFACE,
                   zorder=3, label="standards")
        ax.set_xscale("log"); ax.set_ylim(-0.1, max(s.od_net) * 1.15)
        ax.set_title(pid); ax.set_xlabel("IL-6 (pg/mL)")
    axes[0].set_ylabel("net OD"); axes[0].legend(loc="upper left")
    fig.tight_layout()
"""),
    md("""
**Reading the comparison.**

* **4PL.** On plate 1 every standard recovers within 94–102%. On plates 3 and 4, every standard
  from 15.6 pg/mL up recovers within ±15%, and only the lowest level misses (140% and 118%). On
  plate 2 the three lowest levels are off (67–132%), which is why that plate fails below.
* **Linear.** The straight lines are good only over the part of the range they were
  (apparently) fitted to. Even on plate 1 the line over-recovers the lowest four standards by
  21–34%. On plates 2–4 the lines fit the bottom of the curve: they recover 88–142% at
  7.8–31 pg/mL, then collapse to 28–45% at 500 pg/mL. For these samples, which all sit at the
  bottom of the curve, that part of the error matters less than it might seem. But the fit is
  ±20–40% even there, and the lines extend below the lowest standard, where most samples fall.
* A sigmoid is the right model for a saturable binding assay. A straight line is at best a
  local approximation, and its validity range has to be demonstrated, not assumed.

## 3 · QC: which plates pass, and what range is quantifiable?
"""),
    code("""
result.plate_qc[["plate_id", "passed", "r_squared", "std_levels_pass", "lod", "lloq", "uloq", "reasons"]]
"""),
    code("""
result.wells.loc[result.wells.outlier, ["plate_id", "well", "type", "sample_id", "od", "outlier_reason"]]
"""),
    md("""
* **Plate 2 fails** acceptance: only 4 of 7 standard levels meet the criteria (the rule needs
  ≥ 75%). Its quantification range is 62.5–500 pg/mL, 8× narrower at the low end than plate 1's.
* Plates 1, 3 and 4 pass, with LLOQs of 7.8–15.6 pg/mL. The LOD (blank mean + 3 SD) is about
  3–6 pg/mL. It is estimated from only two zero-standard wells per plate, so treat it as an
  order of magnitude.
* One standard well (plate 3, G2) is excluded as a gross error by the leave-one-out residual
  screen.

## 4 · Where do the samples fall relative to the validated range?
"""),
    code("""
res = result.results.merge(info[["sample_id", "treatment", "time_h", "plate"]], on="sample_id")
res = res.merge(authors, on="sample_id")
order = sorted(res.group.unique(), key=lambda g: (int(g.split()[1][:-1]), g.split()[0] == "ASF"))
summary = (res.groupby("group")
              .agg(plate=("plate", "first"), n=("sample_id", "size"),
                   quantifiable=("censor", lambda c: int((c == "").sum())),
                   below_lloq=("censor", lambda c: int((c == "<LLOQ").sum())),
                   authors_median=("authors_conc", "median"),
                   pipeline_median_extrapolated=("conc_extrap", "median"))
              .reindex(order))
summary
"""),
    code("""
with P.style():
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    for i, t in enumerate(["NSF", "ASF"]):
        d = res[res.treatment == t]
        q = d.censor.eq("")
        ax.scatter(d.authors_conc[~q].clip(lower=0.3), d.conc_extrap[~q].clip(lower=0.3), s=26,
                   facecolor=P.SURFACE, edgecolor=P.SERIES[i], linewidth=1.3,
                   label=f"{t} (below LLOQ)")
        ax.scatter(d.authors_conc[q], d.conc_extrap[q], s=30, color=P.SERIES[i],
                   edgecolor=P.SURFACE, label=f"{t} (quantifiable)")
    lim = [0.3, 120]
    ax.plot(lim, lim, color=P.AXIS, lw=1, zorder=0)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("authors' value (pg/mL, linear fit)")
    ax.set_ylabel("4PL back-calculated (pg/mL)")
    ax.set_title("Same ODs, two calibrations")
    ax.legend(loc="upper left", fontsize=7.5)
    fig.tight_layout()
print(f"Spearman ρ = {res[['authors_conc', 'conc_extrap']].corr('spearman').iloc[0, 1]:.2f}; "
      f"{int(res.censor.eq('<LLOQ').sum())}/{len(res)} samples below LLOQ")
"""),
    md("""
**Most samples (100 of 110) are below the LLOQ.** Only ASF 6 h has several samples in the validated
range (5 of 10). Healthy mouse serum IL-6 normally sits at a few pg/mL, near or below the
sensitivity of this kit, so this is biologically expected. It does mean that most of the numbers
reported for these groups are **extrapolations below the lowest standard**, where neither
calibration has been validated.

The two calibrations agree on the ordering of samples (Spearman ρ above). They disagree on
magnitude, most of all at the bottom. Samples reading at or below the zero standard have no
defined 4PL concentration, and the plot clips them to its lower edge. The authors' linear fits
turn the same wells into values of 1–10 pg/mL, or into hand-entered zeros.

## 5 · Statistics: a strict analysis, then a labelled sensitivity analysis

For each time point with both treatments (1–24 h), ASF is compared with NSF. The pipeline
checks normality (Shapiro-Wilk on log10 values) and picks Welch's t-test or Mann-Whitney. The
five p-values are then **Benjamini-Hochberg-adjusted across time points**.

**(A) Strict, as the pipeline does by default.** Plates that fail QC are excluded. Values below
LLOQ are censored at one **study-wide** limit, so no comparison depends on which plate a sample
happened to be on.
"""),
    code("""
def by_timepoint(cfg, frame=res, qc=result.plate_qc):
    data = st.prepare_values(frame, qc, cfg)
    data = data[data.time_h > 0].copy()
    data["analyte"] = data.time_h.map(lambda t: f"{t} h")
    data["group"] = data.treatment
    out = st.compare_groups(data, cfg, group_order=["NSF", "ASF"])
    return data, out.omnibus.set_index("analyte")

strict_cfg = st.StatsConfig(value_col="conc")
strict_data, strict = by_timepoint(strict_cfg)
display(strict_data.handling.value_counts().rename("samples"))
strict[["test", "note", "p", "p_adj_bh"]]
"""),
    md("""
At validated precision **the data cannot distinguish the groups.** Plate 2 (2 h and NSF 6 h)
is excluded. After harmonisation, *every* remaining value is below the study-wide LLOQ
(46.8 pg/mL, which is 15.6 pg/mL × the 3× dilution used for some 0 h samples), so every value ties
at that limit and no comparison can be tested at all. The 6 h contrast cannot even be formed, because its
control group sits on the failed plate. This is not a failure of the pipeline. It is the correct
answer to the question "what does this assay support quantitatively?"

**(B) Sensitivity analysis: treat sub-LLOQ back-calculated values as data.** This is effectively
what the original analysis did. All plates are used. A value below LLOQ becomes its 4PL
back-calculation. A sample below the curve floor, where no concentration exists, is ranked
lowest. Because those floor substitutions are arbitrary numbers, only **rank-based** inference is
meaningful here. The comparison uses the Mann-Whitney p-value and the rank-biserial correlation
*r*. Fold changes are not used: the substituted floor values make any ratio meaningless.
"""),
    code("""
sens_cfg = st.StatsConfig(value_col="conc", exclude_failed_plates=False, censored="extrapolate")
sens_data, _ = by_timepoint(sens_cfg)
display(sens_data.handling.value_counts().rename("samples"))

rows = []
for t in (1, 2, 6, 12, 24):
    d = sens_data[sens_data.time_h == t]
    x = d.loc[d.treatment == "NSF", "value"].dropna(); y = d.loc[d.treatment == "ASF", "value"].dropna()
    u = mannwhitneyu(y, x, alternative="two-sided")
    rows.append({"time": f"{t} h", "n_NSF": len(x), "n_ASF": len(y),
                 "median_NSF": x.median(), "median_ASF": y.median(),
                 "rank_biserial_r": 2 * u.statistic / (len(x) * len(y)) - 1, "p": u.pvalue})
sens = pd.DataFrame(rows)
sens["p_BH"] = multipletests(sens.p, method="fdr_bh")[1]
sens
"""),
    code("""
with P.style():
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    times = [1, 2, 6, 12, 24]
    for i, t in enumerate(["NSF", "ASF"]):
        d = sens_data[sens_data.treatment == t]
        rng = np.random.default_rng(i)
        for k, tp in enumerate(times):
            v = d[d.time_h == tp]
            x = k + (-0.17 if t == "NSF" else 0.17) + rng.uniform(-0.06, 0.06, len(v))
            measured = v.handling.eq("measured").to_numpy()
            ax.scatter(x[measured], v.value[measured], s=28, color=P.SERIES[i], edgecolor=P.SURFACE, zorder=3)
            ax.scatter(x[~measured], v.value[~measured], s=24, facecolor=P.SURFACE,
                       edgecolor=P.SERIES[i], linewidth=1.2, zorder=3)
            ax.plot([k + (-0.17 if t == "NSF" else 0.17) - 0.12, k + (-0.17 if t == "NSF" else 0.17) + 0.12],
                    [v.value.median()] * 2, color=P.INK, lw=1.8)
        ax.scatter([], [], color=P.SERIES[i], label=t)
    top = sens_data.value.max()
    floor = sens_data.loc[sens_data.handling.str.startswith("below curve floor"), "value"].max()
    for k, r in enumerate(sens.itertuples()):
        star = "***" if r.p_BH < 0.001 else "**" if r.p_BH < 0.01 else "*" if r.p_BH < 0.05 else "ns"
        ax.text(k, top * 2.2, star, ha="center", fontsize=8, color=P.INK_2)
    ax.axhline(floor * 1.7, color=P.AXIS, lw=0.8, ls=(0, (2, 3)))
    ax.text(len(times) - 0.5, floor * 1.25, "below curve floor (ranked lowest)",
            ha="right", va="center", fontsize=7, color=P.INK_2)
    ax.set_yscale("log"); ax.set_ylim(top=top * 5)
    ax.set_xticks(range(len(times)), [f"{t} h" for t in times])
    ax.set_ylabel("serum IL-6 (pg/mL)"); ax.set_xlabel("time after onset")
    ax.set_title("IL-6 time course, sensitivity analysis (hollow = below LLOQ)")
    ax.legend(loc="lower right", bbox_to_anchor=(1, 0.2), ncol=2); ax.grid(axis="x", visible=False)
    fig.tight_layout()
"""),
    md("""
## 6 · Biological interpretation

**What survives.**

* **Ordering.** Under the sensitivity analysis, IL-6 in ASF mice ranks above time-matched NSF
  controls at **6, 12 and 24 h**. Here the rank-biserial *r* is large and the BH-adjusted p is
  below 0.05. There is no difference at 1–2 h. This matches the authors' qualitative story: a
  systemic IL-6 response that appears a few hours after fragmentation begins and persists
  through 24 h.
* **Peak.** ASF 6 h is the only group with a substantial share of samples (5/10, median
  ≈ 27 pg/mL) inside the validated quantification range, so it is the one *quantitative* result
  the assay supports. It is consistent with the authors' ≈ 27 pg/mL. IL-6 peaking hours after an
  inflammatory stimulus is typical of a cytokine cascade. Stress-axis and sympathetic activation
  drive innate-immune cells (and adipose tissue) to secrete IL-6, which then induces the hepatic
  acute-phase response. This fits the paper's conclusion that inflammation from sleep
  fragmentation starts in the periphery.

**What does not survive.**

* **Fold changes.** Fold changes such as "11× at 6 h" or "3× at 24 h" are ratios of
  sub-LLOQ extrapolations. Most control values sit at or below the zero standard, so a ratio
  against them is not a measurement.
* **The 6 h magnitude.** The 6 h comparison crosses plates, and the control plate failed QC.
  Without bridge samples, part of that difference could be calibration.
* **Strict-criteria significance.** Under the strict criteria (A), none of the differences can
  be established.

**How the experiment could be improved.**

1. Use a high-sensitivity IL-6 assay, such as an electrochemiluminescence or Simoa platform with
   an LLOQ below 1 pg/mL, or run serum neat rather than diluted, so that baseline levels are
   quantifiable.
2. **Randomise groups across plates**, never whole groups per plate, and include the same pooled
   serum as a bridge control on every plate. Then `immunoassay.normalize` can remove plate
   effects instead of leaving them confounded with treatment.
3. Pre-specify the censoring policy (substitution, or a Tobit/interval-censored model) instead of
   hand-entering zeros.

*Everything above is regenerated by `python scripts/build_notebooks.py`.*
"""),
]

# --------------------------------------------------------------------------- #
# Notebook 2: ELISAtools SoftMax plates (real QC + normalization)
# --------------------------------------------------------------------------- #
NB2 = [
    md("""
# Real SoftMax Pro plates: outlier QC vs the analyst, and kit-lot normalization

**Data.** Five sandwich-ELISA plates exported from SoftMax Pro. They ship as example data in the
R package **ELISAtools** 0.1.8 (Feng Feng, MIT license; see Feng *et al.*, bioRxiv
[10.1101/483800](https://doi.org/10.1101/483800)). The plates come in two **kit-lot batches**:
three plates in Batch 1 and two in Batch 2. Each plate has 8 standards in triplicate
(3000 → 46.9 pg/mL plus a zero), 15 unknowns in triplicate at 1:10, and the **same QC control**
in triplicate. The analyte is not disclosed. The raw files are in `data/raw/elisatools_feng2019/`.

This dataset has no biology to interpret. It is useful for two engineering questions that come
up on every real ELISA project:

1. The original analyst **masked** wells by hand in SoftMax. How does the pipeline's rule-based
   outlier QC compare?
2. Does a **kit-lot change** shift results, and can the shared QC control correct for it?
"""),
    code(
        setup()
        + """
base = REPO / "data/processed/elisatools_feng2019"
if not (base / "manifest.csv").exists():
    subprocess.run([sys.executable, str(REPO / "scripts/prepare_elisatools.py")], check=True)
manifest = pd.read_csv(base / "manifest.csv")
result = run_pipeline(base / "manifest.csv", PipelineConfig())
manifest[["plate_id", "batch"]]
"""
    ),
    md("""
The SoftMax export stores each plate as one row of 96 values per wavelength (450 and 620 nm).
`immunoassay.io.read_softmax_export` subtracts the 620 nm reference read, which corrects for
plastic and optical imperfections, and parses the `Group:` blocks that hold the real plate map.
The plate map was needed because the package's own annotation file does not match these
plates. One plate even has its top two standard rows swapped, which the group blocks record
correctly.

## 1 · Curves and plate acceptance
"""),
    code("""
display(result.curve_table()[["plate_id", "model", "a", "b", "c", "d", "r_squared"]])
result.plate_qc[["plate_id", "passed", "std_levels_pass", "lod", "lloq", "uloq", "n_flagged", "reasons"]]
"""),
    code("""
k = ("Assay_2_p1", "analyte (undisclosed)")
w = result.wells[result.wells.plate_id == k[0]]
show(P.plot_standard_curve(result.fits[k], w, result.limits.set_index(["plate_id", "analyte"]).loc[k],
                              title="Assay_2 · standard curve", units="pg/mL"))
"""),
    md("""
## 2 · Rule-based outlier QC vs the analyst's masks

These plates run in **triplicate**, so the replicate test is live. The modified z-score is
applied only to groups whose raw-OD CV already exceeds 20%. Without that gate, the MAD of three
values is simply the smaller of two deviations. The ungated rule flagged 23 wells here, only 6
of which the analyst had masked. That was the motivation for the gate.
"""),
    code("""
w = result.wells.copy()
w["analyst_masked"] = w["operator_masked"].astype(str).eq("True")
display(pd.crosstab(w.analyst_masked, w.outlier, rownames=["analyst masked"], colnames=["pipeline outlier"]))

def group_cv(row):
    g = w[(w.plate_id == row.plate_id) & (w.type == row.type) & (w.sample_id == row.sample_id)
          & (w.dilution == row.dilution)]
    return 100 * g.od.std(ddof=1) / g.od.mean()

masked = w[w.analyst_masked | w.outlier].copy()
masked["replicate_group_cv_pct"] = masked.apply(group_cv, axis=1)
masked[["plate_id", "well", "type", "sample_id", "od", "analyst_masked", "outlier",
        "replicate_group_cv_pct"]].sort_values("replicate_group_cv_pct", ascending=False)
"""),
    md("""
**Every well the pipeline removed was also removed by the analyst** (4 of 4). The analyst
removed 13 more, and **all 13** came from replicate groups whose CV was already within 20%. One
(Assay_11_and_12, A3) came from a triplicate with a CV of **1.8%**.
Masking a well in a group that already meets the precision criterion changes the mean very
little. What it does is make the reported CV look better than the measurement was. That kind of
undocumented, discretionary exclusion is what bioanalytical guidance asks labs to replace with a
**pre-specified rule**, and it is what the pipeline enforces.

## 3 · Kit-lot batch effect and bridge normalization

The QC control is the same material on every plate, so it serves as a bridge. If calibration
were perfect, its measured concentration would be identical on all five plates.
"""),
    code("""
ctrl = result.results[result.results.type == "control"].merge(manifest[["plate_id", "batch"]], on="plate_id")
f = result.plate_factors.merge(manifest[["plate_id", "batch"]], on="plate_id")
display(ctrl[["plate_id", "batch", "conc", "conc_norm"]].rename(columns={"conc": "QC control (raw)", "conc_norm": "after normalization"}))
display(f[["plate_id", "batch", "factor"]])
lot = f.groupby("batch")["factor"].apply(lambda x: np.exp(np.log(x).mean()))
print(f"Batch 2 reads {100 * (lot['Batch2'] / lot['Batch1'] - 1):+.0f}% relative to Batch 1 "
      f"(geometric mean of plate factors); inter-plate CV of the control: "
      f"{result.bridge_cv.inter_plate_cv_pct_before.iloc[0]:.1f}% raw")
"""),
    code("""
show(P.plot_normalization(result.results, "analyte (undisclosed)", units="pg/mL"))
"""),
    md("""
**Interpretation.**

* The control's between-plate CV is 22% before correction, comparable to the whole ±20–25%
  acceptance window for a single measurement. The plate factors cluster by lot: the **Batch 2
  plates read about 30% lower** than Batch 1. This is exactly the long-term lot-to-lot drift
  the ELISAtools authors set out to correct. A study that switched kit lots midway would have
  shown a spurious shift in patient values at the switch.
* **Caveat, and it matters.** With a *single* bridge sample, normalization forces the control to
  be identical on every plate. The "after" CV is 0% by construction and cannot validate the
  correction. The factor also absorbs that control's own measurement noise (8–10%
  within-plate CV on these plates). In a real study, use **several bridge levels**, as in the simulated demo: the
  pipeline then fits a two-way log-scale model, reports a standard error per plate factor, and
  the post-correction CV of each bridge becomes an honest check.
* The unknowns P1–P15 are *not* the same patients from plate to plate (the labels are template
  slots), so they cannot be used to check the correction either. `scripts/prepare_elisatools.py`
  therefore prefixes them with the plate ID.
"""),
]


def build(run: bool) -> None:
    NB.mkdir(exist_ok=True)
    for name, cells in (
        ("01_sleep_fragmentation_il6.ipynb", NB1),
        ("02_elisatools_softmax_qc.ipynb", NB2),
    ):
        nb = nbformat.v4.new_notebook(cells=cells)
        nb.metadata["kernelspec"] = {
            "name": "python3",
            "display_name": "Python 3",
            "language": "python",
        }
        if run:
            NotebookClient(
                nb, timeout=600, kernel_name="python3", resources={"metadata": {"path": str(NB)}}
            ).execute()
        nbformat.write(nb, NB / name)
        # Keep generated code cells in the repo's formatting (CI runs `ruff format --check`).
        subprocess.run([sys.executable, "-m", "ruff", "format", str(NB / name)], check=True)
        print(f"wrote notebooks/{name}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-run", action="store_true")
    build(run=not ap.parse_args().no_run)
