"""
AGRO-BHARAT — statistical analysis of bell pepper (Capsicum annuum) stem height
under Neem (Azadirachta indica) extract vs. water control.

Design: 2 groups (between-subjects) x 10 weekly measurements (within-subjects),
n = 10 plants per group. Primary test: two-way mixed (split-plot) ANOVA with
repeated measures on `week`.

Run:  python analysis/analysis.py   (from the repo root)
Outputs go to results/.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pingouin as pg
import statsmodels.formula.api as smf
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "plant_height_long.csv"
OUT = ROOT / "results"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(DATA)
report = []


def log(text=""):
    print(text)
    report.append(text)


# ---------------------------------------------------------------- descriptives
desc = (
    df.groupby(["week", "group"])["height_cm"]
    .agg(n="count", mean="mean", sd=lambda x: x.std(ddof=1), min="min", max="max")
    .reset_index()
)
t_crit = stats.t.ppf(0.975, desc["n"] - 1)
desc["ci95_low"] = desc["mean"] - t_crit * desc["sd"] / np.sqrt(desc["n"])
desc["ci95_high"] = desc["mean"] + t_crit * desc["sd"] / np.sqrt(desc["n"])
desc.round(3).to_csv(OUT / "descriptive_stats.csv", index=False)

log("=" * 70)
log("1. DESCRIPTIVE STATISTICS (cm)")
log("=" * 70)
log(desc.round(2).to_string(index=False))

# ---------------------------------------------------------- mixed ANOVA (main)
aov = pg.mixed_anova(data=df, dv="height_cm", within="week",
                     between="group", subject="plant_id", correction=True,
                     effsize="np2")
aov.round(4).to_csv(OUT / "mixed_anova.csv", index=False)

log("\n" + "=" * 70)
log("2. TWO-WAY MIXED ANOVA (between: group, within: week)")
log("=" * 70)
log(aov.round(4).to_string(index=False))

# Sphericity check on the within factor
spher = pg.sphericity(data=df, dv="height_cm", within="week", subject="plant_id")
log(f"\nMauchly sphericity (week): W = {spher.W:.4f}, p = {spher.pval:.4g}"
    f" -> {'violated, use GG-corrected p' if spher.pval < .05 else 'not violated'}")

# Homogeneity of variances between groups at each week
log("\nLevene test (group variances, per week):")
for w, sub in df.groupby("week"):
    lv = stats.levene(*[g["height_cm"] for _, g in sub.groupby("group")])
    log(f"  week {w:>2}: W = {lv.statistic:.3f}, p = {lv.pvalue:.3f}")

# Normality of residuals (Shapiro-Wilk) per group x week
log("\nShapiro-Wilk normality (per group x week, p-values):")
sw = df.groupby(["group", "week"])["height_cm"].apply(lambda x: stats.shapiro(x).pvalue)
log(sw.unstack("week").round(3).to_string())

# ------------------------------------------------------------- post-hoc tests
posthoc = pg.pairwise_tests(data=df, dv="height_cm", within="week",
                            between="group", subject="plant_id",
                            padjust="holm", effsize="hedges")
inter = posthoc[posthoc["Contrast"] == "week * group"]
inter = inter[["week", "A", "B", "T", "dof", "p_unc", "p_corr", "hedges"]]
inter.round(4).to_csv(OUT / "posthoc_group_by_week.csv", index=False)

log("\n" + "=" * 70)
log("3. POST-HOC: Neem vs Control at each week (Welch t, Holm-adjusted)")
log("=" * 70)
log(inter.round(4).to_string(index=False))

# ------------------------------------------ growth-rate model (mixed effects)
# Random intercept + random slope per plant; tests whether Neem changes slope.
df["week_c"] = df["week"] - 1
df["neem"] = (df["group"] == "neem").astype(int)
lmm = smf.mixedlm("height_cm ~ week_c * neem", df, groups=df["plant_id"],
                  re_formula="~week_c").fit(reml=True)
log("\n" + "=" * 70)
log("4. LINEAR MIXED MODEL: height ~ week * group, random slope per plant")
log("=" * 70)
log(str(lmm.summary()))
b = lmm.params
log(f"Control growth rate : {b['week_c']:.3f} cm/week")
log(f"Neem growth rate    : {b['week_c'] + b['week_c:neem']:.3f} cm/week "
    f"(extra {b['week_c:neem']:.3f} cm/week, p = {lmm.pvalues['week_c:neem']:.4f})")

# Per-plant slopes, compared with a simple t-test (robustness check)
slopes = (df.groupby(["plant_id", "group"])
            .apply(lambda g: np.polyfit(g["week"], g["height_cm"], 1)[0])
            .rename("slope_cm_per_week").reset_index())
slopes.round(4).to_csv(OUT / "per_plant_growth_slopes.csv", index=False)
s_c = slopes.query("group=='control'")["slope_cm_per_week"]
s_n = slopes.query("group=='neem'")["slope_cm_per_week"]
tt = stats.ttest_ind(s_n, s_c, equal_var=False)
mw = stats.mannwhitneyu(s_n, s_c, alternative="two-sided")
log(f"\nPer-plant slopes: control {s_c.mean():.3f} ± {s_c.std():.3f}, "
    f"neem {s_n.mean():.3f} ± {s_n.std():.3f} cm/week")
log(f"  Welch t = {tt.statistic:.3f}, p = {tt.pvalue:.4f}; "
    f"Mann-Whitney U = {mw.statistic:.1f}, p = {mw.pvalue:.4f}")

(OUT / "statistical_report.txt").write_text("\n".join(report) + "\n")

# -------------------------------------------------------------------- figures
COL = {"control": "#C98A1B", "neem": "#2B7A9E"}
LAB = {"control": "Control (water)", "neem": "Neem extract"}

fig, ax = plt.subplots(figsize=(7, 4.5))
for g, sub in df.groupby("group"):
    for _, p in sub.groupby("plant_id"):
        ax.plot(p["week"], p["height_cm"], color=COL[g], alpha=0.18, lw=1)
    d = desc[desc["group"] == g]
    ax.fill_between(d["week"], d["ci95_low"], d["ci95_high"], color=COL[g], alpha=0.2)
    ax.plot(d["week"], d["mean"], color=COL[g], lw=2.5, marker="o", label=LAB[g])
ax.set(xlabel="Week", ylabel="Stem height (cm)", xticks=range(1, 11),
       title="Plant growth curves by treatment (mean ± 95% CI)")
ax.spines[["top", "right"]].set_visible(False)
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(FIG / "growth_curves.png", dpi=200)

fig, ax = plt.subplots(figsize=(8, 4.5))
w = 0.38
for i, g in enumerate(["control", "neem"]):
    d = desc[desc["group"] == g]
    ax.bar(d["week"] + (i - 0.5) * w, d["mean"], w, color=COL[g], label=LAB[g],
           yerr=d["sd"], capsize=2, error_kw=dict(lw=0.8))
ax.set(xlabel="Week", ylabel="Mean stem height (cm)", xticks=range(1, 11),
       title="Mean plant height per week (± SD)")
ax.spines[["top", "right"]].set_visible(False)
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(FIG / "weekly_means.png", dpi=200)

fig, ax = plt.subplots(figsize=(4.5, 4.5))
ax.boxplot([s_c, s_n], tick_labels=[LAB["control"], LAB["neem"]], widths=0.5)
for i, s in enumerate([s_c, s_n], 1):
    ax.scatter(np.full(len(s), i) + np.random.default_rng(0).uniform(-.08, .08, len(s)),
               s, color=list(COL.values())[i - 1], zorder=3)
ax.set(ylabel="Growth rate (cm/week)", title="Per-plant growth slopes")
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(FIG / "growth_slopes.png", dpi=200)

print("\nDone. Results written to results/.")
