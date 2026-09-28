"""Pillar 2 Step 9 — journal-style figures (SVG+PNG), reusing frozen plotting_style.
Figures:
  fig_pillar2_gsea_ladder      : PT AKI-vs-normal GSEA NES ladder (EV vs controls) - specificity failure
  fig_pillar2_permutation      : observed vs matched null, EV/ribosomal/mito across PT contrasts
  fig_pillar2_module_violin    : PT EV module score donor distributions (AKI/CKD/normal)
  fig_pillar2_composition_sev  : PT abundance (no shift) + EV-score vs injury-score scatter (severity confound)
"""
import os, sys, json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = ["Liberation Sans","Arimo","DejaVu Sans"]
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt
sys.path.insert(0, "/mnt/results/ev_aki_ckd_project/00_preregistration/code")
from plotting_style import set_journal_style, save_figure, COLORS, SINGLE_COL, ONEHALF_COL, DOUBLE_COL, MM
set_journal_style()

OUT = "/mnt/results/ev_aki_ckd_project/02_singlecell_specificity/data"
FIGSUB = "02_singlecell_specificity/figures"
gsea = pd.read_csv(f"{OUT}/sender_ev_enrichment_gsea.tsv", sep="\t")
perm = pd.read_csv(f"{OUT}/permutation_null_sender.tsv", sep="\t")
mod = pd.read_csv(f"{OUT}/pt_ev_module_score_donor.tsv", sep="\t")
sev = pd.read_csv(f"{OUT}/severity_within_PT.tsv", sep="\t")
ab = pd.read_csv(f"{OUT}/pt_abundance_by_donor.tsv", sep="\t")

LADDER = ["EV_TIER1_PRIMARY","EV_TIER2_CARGO","NONEV_SECRETED_BACKGROUND",
          "CTRL_LYSO_ENDOSOMAL","CTRL_MITOCHONDRIAL","CTRL_RIBOSOMAL"]
LAB = {"EV_TIER1_PRIMARY":"EV Tier-1\n(machinery)","EV_TIER2_CARGO":"EV Tier-2\n(cargo)",
       "NONEV_SECRETED_BACKGROUND":"Non-EV\nsecreted","CTRL_LYSO_ENDOSOMAL":"Lyso/\nendosomal",
       "CTRL_MITOCHONDRIAL":"Mito\n(ctrl)","CTRL_RIBOSOMAL":"Ribosomal\n(ctrl)"}

# ---------- Fig 1: GSEA NES ladder, PT AKI-vs-normal ----------
sub = gsea[(gsea.compartment=="PT")&(gsea.contrast=="AKI_vs_normal")].set_index("term").loc[LADDER]
fig, ax = plt.subplots(figsize=(ONEHALF_COL, 62*MM))
cols = [COLORS["ev"] if t=="EV_TIER1_PRIMARY" else (COLORS["background"] if t.startswith("CTRL") else "#7F7F7F") for t in LADDER]
cols = []
for t in LADDER:
    if t=="EV_TIER1_PRIMARY": cols.append(COLORS["ev"])
    elif t in ("CTRL_RIBOSOMAL","CTRL_MITOCHONDRIAL"): cols.append("#CC5500")  # negative controls
    else: cols.append(COLORS["background"])
bars = ax.bar(range(len(LADDER)), sub["nes"].values, color=cols, edgecolor="black", linewidth=0.5)
ax.axhline(0, color="black", lw=0.6)
ax.axhline(1.0, color="grey", ls="--", lw=0.6)
nmin = float(sub["nes"].min()); nmax = float(sub["nes"].max())
ax.set_ylim(nmin-0.55, nmax+0.45)   # headroom so labels never collide with axis/ticks
for i,(t,row) in enumerate(sub.iterrows()):
    star = "*" if row["fdr"]<0.05 else ""
    if row["nes"]>=0:
        ax.text(i, row["nes"]+0.06, f"{row['nes']:.2f}{star}", ha="center", va="bottom", fontsize=6)
    else:
        # place just below the bar bottom, with clear gap; axis lower limit gives room
        ax.text(i, row["nes"]-0.08, f"{row['nes']:.2f}{star}", ha="center", va="top", fontsize=6)
ax.set_xticks(range(len(LADDER))); ax.set_xticklabels([LAB[t] for t in LADDER], fontsize=6)
ax.set_ylabel("GSEA NES"); ax.set_title("PT (injured proximal tubule): AKI vs normal", fontsize=8)
ax.text(0.5, -0.30, "EV Tier-1 (green) is weak and NON-SPECIFIC: ribosomal & non-EV secreted controls (orange/grey) enrich more.\n* FDR<0.05.  Dashed line = |NES|=1 effect-size floor.",
        transform=ax.transAxes, ha="center", va="top", fontsize=5.2, color="#333333")
fig.tight_layout(rect=[0,0.02,1,1])
save_figure(fig, "fig_pillar2_gsea_ladder", subdir=FIGSUB, dpi=600); plt.close(fig)

# ---------- Fig 2: permutation OBSERVED vs matched NULL (redesigned for clarity) ----------
# Encoding: matched-null 95% CI drawn as a grey whisker band at the NULL location;
# OBSERVED value drawn as a distinct filled diamond in the set colour. The reader sees
# whether the observed enrichment sits ABOVE its own matched null (the specificity test),
# which is exactly what the one-sided p_greater measures. No bar-to-observed overlap.
from matplotlib.lines import Line2D
contrasts = ["AKI_vs_normal","CKD_vs_normal","AKI_vs_CKD"]
terms = ["EV_TIER1_PRIMARY","CTRL_RIBOSOMAL","CTRL_MITOCHONDRIAL"]
tcol = {"EV_TIER1_PRIMARY":COLORS["ev"],"CTRL_RIBOSOMAL":"#CC5500","CTRL_MITOCHONDRIAL":"#E69F00"}
NULL_GREY = "#9E9E9E"
fig, axes = plt.subplots(1, 3, figsize=(DOUBLE_COL, 58*MM), sharey=False)
for j,cn in enumerate(contrasts):
    ax = axes[j]; x = np.arange(len(terms))
    rows_cn = [perm[(perm.compartment=="PT")&(perm.contrast==cn)&(perm.term==t)].iloc[0] for t in terms]
    all_vals = ([r["observed"] for r in rows_cn]
                + [r["null_mean"]+1.96*r["null_sd"] for r in rows_cn]
                + [r["null_mean"]-1.96*r["null_sd"] for r in rows_cn])
    ymax, ymin = max(all_vals), min(all_vals)
    pad = (ymax-ymin)*0.18
    ax.set_ylim(ymin-pad*0.8, ymax+pad*1.8)
    lab_y = ymax + pad*0.7     # fixed label row above everything
    dx = 0.16  # horizontal offset so observed diamond never overlaps the null whisker
    for i,(t,r) in enumerate(zip(terms, rows_cn)):
        lo, hi = r["null_mean"]-1.96*r["null_sd"], r["null_mean"]+1.96*r["null_sd"]
        xn = i - dx   # null whisker slightly left
        xo = i + dx   # observed diamond slightly right
        # matched-null 95% CI as a grey vertical whisker with caps + null-mean tick
        ax.plot([xn, xn], [lo, hi], color=NULL_GREY, lw=3.0, solid_capstyle="round", zorder=1, alpha=0.85)
        ax.plot([xn-0.08, xn+0.08], [r["null_mean"], r["null_mean"]], color=NULL_GREY, lw=1.0, zorder=2)
        # observed as a filled diamond in set colour, black-edged
        ax.scatter([xo], [r["observed"]], marker="D", s=46, color=tcol[t],
                   edgecolor="black", linewidth=0.6, zorder=3)
        # thin connector to make the observed-vs-null comparison explicit
        ax.plot([xn, xo], [r["observed"], r["observed"]], color=tcol[t], lw=0.5, ls=":", zorder=1, alpha=0.6)
        p = r["p_greater"]
        ax.text(i, lab_y, f"p={p:.0e}" if p<0.01 else f"p={p:.2f}", ha="center", va="bottom", fontsize=5.5)
    ax.axhline(0, color="black", lw=0.5, zorder=0)
    ax.set_xlim(-0.6, len(terms)-0.4)
    ax.set_xticks(x); ax.set_xticklabels(["EV\nTier-1","Ribo\n(ctrl)","Mito\n(ctrl)"], fontsize=6)
    ax.set_title(cn.replace("_"," "), fontsize=7.5)
    if j==0: ax.set_ylabel("EV/control set score")
# single shared legend BELOW the panels (outside plot area) to avoid p-label collisions
leg = [Line2D([0],[0], marker="D", color="none", markerfacecolor="#555555",
              markeredgecolor="black", markersize=7, label="observed set score"),
       Line2D([0],[0], color=NULL_GREY, lw=3, label="matched-null 95% CI")]
fig.legend(handles=leg, loc="lower center", ncol=2, fontsize=6, frameon=False,
           handletextpad=0.5, columnspacing=1.6, bbox_to_anchor=(0.5, -0.02))
fig.suptitle("PT matched-background permutation (B=10,000): observed vs matched null (EV specificity test)", y=1.02, fontsize=8.5)
fig.tight_layout(rect=[0,0.05,1,0.93])
save_figure(fig, "fig_pillar2_permutation", subdir=FIGSUB, dpi=600); plt.close(fig)

# ---------- Fig 3: PT EV module-score violins ----------
order = ["normal","CKD","AKI"]; pal = {"normal":"#B9B9B9","CKD":COLORS["ckd"],"AKI":COLORS["aki"]}
fig, ax = plt.subplots(figsize=(SINGLE_COL, 68*MM))
data = [mod.loc[mod.disease_code==g,"module_score"].dropna().values for g in order]
parts = ax.violinplot(data, showmeans=False, showmedians=True, widths=0.8)
for pc,g in zip(parts["bodies"], order):
    pc.set_facecolor(pal[g]); pc.set_alpha(0.6); pc.set_edgecolor("black"); pc.set_linewidth(0.5)
for key in ["cmedians","cbars","cmins","cmaxes"]:
    if key in parts: parts[key].set_color("black"); parts[key].set_linewidth(0.8)
for i,g in enumerate(order):
    v = mod.loc[mod.disease_code==g,"module_score"].dropna().values
    ax.scatter(np.random.normal(i+1,0.05,len(v)), v, s=4, color="black", alpha=0.35, zorder=3)
ax.set_xticks(range(1,4)); ax.set_xticklabels([f"{g}\n(n={ (mod.disease_code==g).sum() })" for g in order], fontsize=6.5)
ax.set_ylabel("PT EV Tier-1 module score\n(donor-level, z-standardized)")
ax.set_title("Per-donor PT EV module score", fontsize=8)
ax.text(0.5,-0.22,"Supporting/non-gating. Donor-consistent shift, but see severity confound (next).",
        transform=ax.transAxes, ha="center", va="top", fontsize=5.2, color="#333333")
fig.tight_layout(rect=[0,0.02,1,1])
save_figure(fig, "fig_pillar2_module_violin", subdir=FIGSUB, dpi=600); plt.close(fig)

# ---------- Fig 4: composition (no shift) + severity confound scatter ----------
fig, axes = plt.subplots(1, 2, figsize=(DOUBLE_COL, 62*MM))
# 4a abundance
ax = axes[0]
for i,g in enumerate(order):
    v = ab.loc[ab.disease_code==g,"pt_frac"].dropna().values
    ax.scatter(np.random.normal(i,0.06,len(v)), v, s=5, color=pal[g], alpha=0.6, edgecolor="none")
    ax.hlines(np.median(v), i-0.2, i+0.2, color="black", lw=1.2)
ax.set_xticks(range(3)); ax.set_xticklabels(order, fontsize=6.5)
ax.set_ylabel("PT fraction of donor's cells"); ax.set_title("PT abundance: no significant shift\n(AKI vs normal p=0.19)", fontsize=7.5)
# 4b severity confound
ax = axes[1]
for g in order:
    s = sev[sev.disease_code==g]
    ax.scatter(s["injury_score"], s["module_score"], s=9, color=pal[g], alpha=0.65, edgecolor="none", label=g)
# regression line all donors
from numpy.polynomial import polynomial as P
x = sev["injury_score"].values; y = sev["module_score"].values
b = np.polyfit(x,y,1); xs=np.linspace(x.min(),x.max(),50); ax.plot(xs, np.polyval(b,xs), color="black", lw=1, ls="--")
from scipy.stats import spearmanr
rho_all,_ = spearmanr(x,y); rho_aki,_ = spearmanr(sev.loc[sev.disease_code=="AKI","injury_score"], sev.loc[sev.disease_code=="AKI","module_score"])
ax.text(0.04,0.96,f"all ρ={rho_all:.2f}\nwithin-AKI ρ={rho_aki:.2f}", transform=ax.transAxes, va="top", fontsize=6)
ax.set_xlabel("PT injury-severity panel score"); ax.set_ylabel("PT EV Tier-1 module score")
ax.set_title("EV score entangled with injury severity\n(replicates Pillar-1 confound)", fontsize=7.5)
ax.legend(fontsize=5.5, frameon=False, loc="lower right")
fig.tight_layout()
save_figure(fig, "fig_pillar2_composition_severity", subdir=FIGSUB, dpi=600); plt.close(fig)

print("Saved 4 Pillar-2 figures to "+FIGSUB)
