"""SpaPath pi versus the true pi (350 gene sets, 7 seeds): one scatter per effect
distribution coloured by K, and the median and interquartile range of the true and
estimated pi in each group.
Usage: python plot_simulation_pi_accuracy.py SIM_ROOT OUT.png
SIM_ROOT contains tables_input/per_pathway_all.csv (tables_input_v3.R)."""
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

root = Path(sys.argv[1])
rows = list(csv.DictReader(open(root / "tables_input" / "per_pathway_all.csv", encoding="utf-8")))
assert len(rows) == 350
GROUPS = [("K0", "0", "#9e9e9e"), ("K1_10", "1\u201310", "#a6d543"), ("K11_30", "11\u201330", "#21a685"),
          ("K31_60", "31\u201360", "#2e6e8e"), ("K61_100", "61\u2013100", "#440d54")]
PROFILES = [("broad_weak", "Broad, weak effects", "#4292c6"), ("narrow_strong", "Narrow, strong effects", "#08306b")]

fig, axes = plt.subplots(1, 3, figsize=(11, 3.7))
for ax, (pkey, pname, _) in zip(axes[:2], PROFILES):
    d = [r for r in rows if r["effect_profile"] == pkey]
    for gkey, glab, col in GROUPS:
        g = [r for r in d if r["response_group"] == gkey]
        ax.scatter([float(r["pi_true"]) for r in g], [float(r["spapath_pi"]) for r in g], s=16,
                   color=col, edgecolor="#555555", linewidth=.3, label=glab)
    mae = np.mean([abs(float(r["spapath_pi"]) - float(r["pi_true"])) for r in d])
    ax.plot([.4, .95], [.4, .95], "k--", lw=.8)
    ax.text(.42, .91, f"MAE = {mae:.3f}", fontsize=9)
    ax.set_xlim(.4, .95); ax.set_ylim(.4, .95)
    ax.set_xlabel("True \u03c0"); ax.set_ylabel("SpaPath \u03c0"); ax.set_title(pname)
    ax.spines[["top", "right"]].set_visible(False)
axes[1].legend(title="Responding genes K", frameon=False, loc="lower right", fontsize=8, title_fontsize=8)

ax = axes[2]
xs = np.arange(len(GROUPS))
for i, (pkey, pname, col) in enumerate(PROFILES):
    for j, (gkey, _, _) in enumerate(GROUPS):
        g = [r for r in rows if r["effect_profile"] == pkey and r["response_group"] == gkey]
        for k, (field, marker, face) in enumerate((("pi_true", "D", "white"), ("spapath_pi", "o", col))):
            v = np.array([float(r[field]) for r in g])
            q1, med, q3 = np.percentile(v, [25, 50, 75])
            x = j + (i * 2 + k - 1.5) * .12
            ax.plot([x, x], [q1, q3], color=col, lw=1)
            ax.scatter([x], [med], marker=marker, s=34, facecolor=face, edgecolor=col, linewidth=1.2, zorder=3,
                       label=f"{pname}, {'true' if k == 0 else 'SpaPath'} \u03c0" if j == 0 else None)
ax.set_xticks(xs); ax.set_xticklabels([g[1] for g in GROUPS])
ax.set_ylim(.4, .95); ax.set_xlabel("Responding genes K"); ax.set_ylabel("\u03c0 (median and interquartile range)")
ax.set_title("True and estimated \u03c0 by group")
ax.legend(frameon=False, fontsize=7.5, loc="upper left")
ax.spines[["top", "right"]].set_visible(False)
fig.suptitle("SpaPath \u03c0 versus the true \u03c0 across 350 simulated gene sets (7 seeds)", fontsize=11)
fig.tight_layout()
out = Path(sys.argv[2])
fig.savefig(out, dpi=200)
allmae = np.mean([abs(float(r["spapath_pi"]) - float(r["pi_true"])) for r in rows])
print(out, "overall MAE", round(allmae, 4))
for pkey, pname, _ in PROFILES:
    print(pname, round(np.mean([abs(float(r["spapath_pi"]) - float(r["pi_true"])) for r in rows if r["effect_profile"] == pkey]), 4))
