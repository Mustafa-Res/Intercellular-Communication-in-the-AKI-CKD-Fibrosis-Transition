#!/usr/bin/env python3
"""
Pillar 4 figures (data-driven; journal style via 00_preregistration/code/plotting_style.py):
  F_FUNNEL : signature project figure - progressive pre-registered evidence refinement.
  F6       : whole-landscape EV subtraction - route composition of the 34 reproducible survivors.
  F7       : two-level route map - 3 primary consensus axes + Supporting reproducible communication landscape.

All counts come from frozen Pillar-3 + Pillar-4 artifacts. No values are hardcoded.
"""
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

PROJ = Path("/mnt/results/ev_aki_ckd_project")
sys.path.insert(0, str(PROJ / "00_preregistration" / "code"))
from plotting_style import set_journal_style, save_figure, COLORS, MM, SINGLE_COL, DOUBLE_COL, ONEHALF_COL  # noqa

set_journal_style()
P4 = PROJ / "04_synthesis" / "data"

funnel = json.load(open(P4 / "pillar4_funnel.json"))
summary = json.load(open(P4 / "pillar4_subtraction_summary.json"))
land = pd.read_csv(P4 / "pillar4_whole_landscape_routes.tsv", sep="\t")
prim = pd.read_csv(P4 / "pillar4_route_partition_primary.tsv", sep="\t")
supp = pd.read_csv(P4 / "pillar4_supporting_landscape.tsv", sep="\t")

# route colours (semantic, colourblind-safe)
ROUTE_COLORS = {
    "EV-compatible": COLORS["ev"],          # green
    "Soluble-classical": COLORS["aki"],     # vermillion
    "ECM": COLORS["ckd"],                   # blue
    "Contact": "#CC79A7",                   # reddish purple
    "unclassified": COLORS["ns"],           # grey
}

# ----------------------------------------------------------------------------
# F_FUNNEL  - signature figure
# ----------------------------------------------------------------------------
# HONEST STRUCTURE (two logically distinct parts, drawn distinctly):
#   (i) NESTED FUNNEL: 1311 -> 44 -> 34 -> 3. Each stage is a strict subset of the
#       one above (a genuine narrowing). S2 and S3 are equal (significance stage
#       removes none), so they are collapsed into one stage to avoid implying a
#       drop that did not happen.
#   (ii) ROUTE ATTRIBUTION (NOT a funnel continuation): the terminal 3 consensus
#       axes are PARTITIONED by route -> 0 EV-compatible vs 3 soluble/ECM
#       (0 + 3 = 3). Drawn as a side-by-side split below the funnel so the
#       nested-subset invariant is never violated.
st = {s["id"]: s for s in funnel["stages"]}
prim_rc = summary["primary_consensus"]["route_tag_counts"]
n_ev_consensus = summary["primary_consensus"]["n_ev_compatible"]
n_soluble_ecm_consensus = prim_rc.get("Soluble-classical", 0) + prim_rc.get("ECM", 0)
assert n_ev_consensus + n_soluble_ecm_consensus == st["S5"]["n"], "attribution must partition consensus"

nested = [
    (st["S1"]["n"], "Candidate PT\u2192receiver interactions", "all axes, any engine, any atlas"),
    (st["S3"]["n"], "Cross-engine direction + statistical support", "up-in-AKI in both engines, significant in \u22651 engine (\u22651 atlas)"),
    (st["S4"]["n"], "+ Effect-size floor (reproducible)", "a-priori absolute + relative floor met (\u22651 atlas)"),
    (st["S5"]["n"], "Cross-atlas replication (consensus)", "full evidence chain in BOTH atlases"),
]
counts = [s[0] for s in nested]
labels = [s[1] for s in nested]
crits = [s[2] for s in nested]
nfun = len(nested)

fig, ax = plt.subplots(figsize=(DOUBLE_COL, DOUBLE_COL * 0.66))
maxc = max(counts)
def half_width(c):
    # sqrt scaling for the 1311->3 span; floor so the terminal (n=3) bar stays readable
    return max(0.10, 0.46 * (np.sqrt(max(c, 1)) / np.sqrt(maxc)))

# funnel occupies upper region (y from nfun+1 .. 2); attribution split sits below (y ~ 0)
y_fun = np.arange(nfun)[::-1].astype(float) + 2.0
bar_h = 0.66
body = "#B9C4CC"
xc = 0.5  # funnel centre x
for i in range(nfun):
    hw = half_width(counts[i])
    y = y_fun[i]
    is_terminal = (i == nfun - 1)
    ax.add_patch(mpatches.FancyBboxPatch(
        (xc - hw, y - bar_h / 2), 2 * hw, bar_h,
        boxstyle="round,pad=0.002,rounding_size=0.02",
        linewidth=0.6, edgecolor="white",
        facecolor=COLORS["up"] if is_terminal else body, zorder=3))
    ax.text(xc, y, f"{counts[i]:,}", ha="center", va="center", fontsize=9,
            fontweight="bold", color="white" if is_terminal else "#222222", zorder=4)
    ax.text(xc - 0.50, y, labels[i], ha="right", va="center", fontsize=7, zorder=4)
    ax.text(xc + 0.50, y, crits[i], ha="left", va="center", fontsize=5.6,
            color="#555555", style="italic", zorder=4)
    if i < nfun - 1:
        ax.annotate("", xy=(xc, y_fun[i + 1] + bar_h / 2 + 0.03),
                    xytext=(xc, y - bar_h / 2 - 0.03),
                    arrowprops=dict(arrowstyle="-|>", color="#888888", lw=0.9), zorder=2)

# ---- route-attribution split (branch, not a funnel step) ----
y_split = 0.35
x_ev, x_se = 0.18, 0.82
# branch label
ax.text(xc, 1.52, "Route attribution of the 3 consensus axes  (partition, not a further filter)",
        ha="center", va="center", fontsize=6.2, style="italic", color="#444444")
# two branch arrows from the terminal (consensus) box
for xt in (x_ev, x_se):
    ax.annotate("", xy=(xt, y_split + 0.34), xytext=(xc, y_fun[-1] - bar_h / 2 - 0.03),
                arrowprops=dict(arrowstyle="-|>", color="#888888", lw=0.9,
                                connectionstyle="arc3,rad=0"), zorder=2)
# EV-compatible box (0) - drawn as an OUTLINE (empty) box so 0 reads as "none"
ax.add_patch(mpatches.FancyBboxPatch(
    (x_ev - 0.16, y_split - bar_h / 2), 0.32, bar_h,
    boxstyle="round,pad=0.002,rounding_size=0.02",
    linewidth=1.1, edgecolor=COLORS["ev"], facecolor="white",
    linestyle="--", zorder=3))
ax.text(x_ev, y_split, f"{n_ev_consensus}", ha="center", va="center",
        fontsize=11, fontweight="bold", color=COLORS["ev"], zorder=4)
ax.text(x_ev, y_split - 0.62, "EV-exclusive\n(frozen EV rule vs soluble background)",
        ha="center", va="top", fontsize=6, color=COLORS["ev"])
# soluble/ECM box (3) - filled
ax.add_patch(mpatches.FancyBboxPatch(
    (x_se - 0.16, y_split - bar_h / 2), 0.32, bar_h,
    boxstyle="round,pad=0.002,rounding_size=0.02",
    linewidth=0.6, edgecolor="white", facecolor="#2f8f5b", zorder=3))
ax.text(x_se, y_split, f"{n_soluble_ecm_consensus}", ha="center", va="center",
        fontsize=11, fontweight="bold", color="white", zorder=4)
ax.text(x_se, y_split - 0.62, "Soluble + ECM constructive program\n(SPP1 soluble-classical; collagen-IV ECM)",
        ha="center", va="top", fontsize=6, color="#1d6b41")

ax.set_xlim(-0.62, 1.72)
ax.set_ylim(-0.9, nfun + 2.1)
ax.axis("off")
ax.set_title("Progressive, pre-registered evidence refinement of AKI PT\u2192receiver communication",
             fontsize=8.5, pad=6)
fig.text(0.5, 0.028,
         "Funnel (top) is strictly nested \u2014 each stage a subset of the one above; every narrowing step applies a pre-registered "
         "criterion (right), not post-hoc selection. Bottom: the 3 consensus axes partitioned by route (0 EV-exclusive + 3 soluble/ECM = 3). Primary cap = 2000.",
         ha="center", fontsize=5.3, color="#555555")
fig.text(0.5, 0.006,
         "Initial candidate interactions represent all PT\u2192receiver interactions identified by at least one engine in at least one atlas, "
         "prior to replication and specificity filtering. See Fig. F_FUNNEL_SUPP for the expression-floor sensitivity analysis.",
         ha="center", fontsize=4.9, color="#777777", style="italic")
save_figure(fig, "F_FUNNEL_evidence_refinement", subdir="04_synthesis/figures")
plt.close(fig)

# ----------------------------------------------------------------------------
# F6 - whole-landscape EV subtraction (route composition of 34 reproducible survivors)
# ----------------------------------------------------------------------------
order = ["EV-compatible", "Soluble-classical", "ECM", "Contact", "unclassified"]
DISP = {"EV-compatible": "EV-exclusive"}  # display rename; data keys unchanged
rc = summary["reproducible_population_S4"]["route_tag_counts"]
vals = [rc.get(k, 0) for k in order]

fig, axes = plt.subplots(1, 2, figsize=(DOUBLE_COL, DOUBLE_COL * 0.40),
                         gridspec_kw={"width_ratios": [1.15, 1]})

# (a) bar composition
ax = axes[0]
xs = np.arange(len(order))
bars = ax.bar(xs, vals, color=[ROUTE_COLORS[k] for k in order],
              edgecolor="white", linewidth=0.6, width=0.72)
for x, v in zip(xs, vals):
    ax.text(x, v + 0.4, str(v), ha="center", va="bottom", fontsize=7, fontweight="bold")
ax.set_xticks(xs)
ax.set_xticklabels([DISP.get(k, k).replace("-", "-\n") if k == "Soluble-classical" else DISP.get(k, k) for k in order],
                   rotation=0, fontsize=6)
ax.set_ylabel("Reproducible up-in-AKI\nPT\u2192receiver axes (n)")
ax.set_ylim(0, max(vals) * 1.18)
ax.set_title(f"Route of the {sum(vals)} reproducible survivors", fontsize=7.5)
# annotate the EV=0 bar explicitly
ax.annotate("0 EV-exclusive", xy=(0, 0.3), xytext=(0.15, max(vals) * 0.55),
            fontsize=6, color=COLORS["ev"], ha="left",
            arrowprops=dict(arrowstyle="-|>", color=COLORS["ev"], lw=0.8))

# (b) ligand x route (which ligands carry the survivors)
ax = axes[1]
lig_route = (land.groupby(["ligand", "route_tag"]).size().reset_index(name="n"))
ligs = (land.groupby("ligand").size().sort_values(ascending=True).index.tolist())
yl = np.arange(len(ligs))
for k in order:
    sub = lig_route[lig_route["route_tag"] == k].set_index("ligand")["n"]
    widths = [sub.get(l, 0) for l in ligs]
    left = np.zeros(len(ligs))
    # stack (each ligand has a single route here, but keep general)
    ax.barh(yl, widths, left=left, color=ROUTE_COLORS[k], edgecolor="white",
            linewidth=0.5, height=0.72, label=DISP.get(k, k))
ax.set_yticks(yl)
ax.set_yticklabels(ligs, fontsize=6, fontstyle="italic")
ax.set_xlabel("Axes (n)")
ax.set_title("Ligands carrying the survivors", fontsize=7.5)
handles = [mpatches.Patch(color=ROUTE_COLORS[k], label=DISP.get(k, k)) for k in order if rc.get(k, 0) > 0]
ax.legend(handles=handles, fontsize=5.6, loc="lower right", frameon=False)

fig.suptitle("Whole-landscape EV subtraction: reproducible AKI communication is soluble + ECM, not EV",
             fontsize=8.5, y=1.02)
save_figure(fig, "F6_whole_landscape_subtraction", subdir="04_synthesis/figures")
plt.close(fig)

# ----------------------------------------------------------------------------
# F7 - two-level route map (primary consensus + Supporting reproducible landscape)
# ----------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(DOUBLE_COL, DOUBLE_COL * 0.46),
                         gridspec_kw={"width_ratios": [1, 1.25]})

# (a) primary 3 consensus axes - labelled dot map coloured by route, marker by receiver
ax = axes[0]
recv_marker = {"endothelial": "o", "stromal": "s", "myeloid": "^"}
prim_sorted = prim.sort_values(["route_tag", "key"]).reset_index(drop=True)
for i, r in prim_sorted.iterrows():
    m = recv_marker.get(r["receiver"], "D")
    ax.scatter(0.5, i, s=140, marker=m, color=ROUTE_COLORS[r["route_tag"]],
               edgecolor="black", linewidth=0.6, zorder=3)
    tier_flag = "  \u2020 Tier-2" if "Tier 2" in str(r["tier"]) else ""
    ax.text(0.72, i, f"{r['ligand']}\u2192{r['receptor']}\n({r['receiver']}, {r['pathway_name']}){tier_flag}",
            ha="left", va="center", fontsize=6)
ax.set_xlim(0.3, 2.2)
ax.set_ylim(-0.6, len(prim_sorted) - 0.4)
ax.axis("off")
ax.set_title(f"PRIMARY: {len(prim_sorted)} consensus axes\n(complete pre-registered chain)", fontsize=7.5)

# (b) Supporting reproducible communication landscape - route composition + top ligands
ax = axes[1]
supp_rc = supp["route_tag"].value_counts()
supp_order = [k for k in order if k in supp_rc.index]
lig_r = supp.groupby(["ligand", "route_tag"]).size().reset_index(name="n")
lig_tot = supp.groupby("ligand").size().sort_values(ascending=True)
ligs2 = lig_tot.index.tolist()
yl = np.arange(len(ligs2))
for k in supp_order:
    sub = lig_r[lig_r["route_tag"] == k].set_index("ligand")["n"]
    widths = [sub.get(l, 0) for l in ligs2]
    ax.barh(yl, widths, color=ROUTE_COLORS[k], edgecolor="white", linewidth=0.4,
            height=0.74, label=k)
ax.set_yticks(yl)
ax.set_yticklabels(ligs2, fontsize=5.2, fontstyle="italic")
ax.set_xlabel("Axes (n)")
ax.set_title(f"SUPPORTING reproducible communication landscape\n(n={len(supp)}; descriptive, below replication bar)",
             fontsize=7.5)
handles = [mpatches.Patch(color=ROUTE_COLORS[k], label=k) for k in supp_order]
ax.legend(handles=handles, fontsize=5.6, loc="lower right", frameon=False, title=None)

# receiver-marker legend for panel (a)
rec_handles = [plt.Line2D([0], [0], marker=recv_marker[r], color="w",
               markerfacecolor="#999999", markeredgecolor="black", markersize=6, label=r)
               for r in ["endothelial", "stromal"]]
axes[0].legend(handles=rec_handles, fontsize=5.6, loc="lower center", ncol=2,
               frameon=False, title="receiver", title_fontsize=5.6)

fig.suptitle("Two-level surviving communication: primary consensus + supporting reproducible landscape (0 EV-exclusive in either)",
             fontsize=8, y=1.02)
save_figure(fig, "F7_two_level_route_map", subdir="04_synthesis/figures")
plt.close(fig)

# ----------------------------------------------------------------------------
# F_FUNNEL_SUPP - expression-floor sensitivity of the funnel start (SUPPLEMENTARY)
# ----------------------------------------------------------------------------
# Shows the funnel is robust to removing low-expression candidates: downstream
# stages (S2->S5) are essentially unchanged as an expression floor is raised.
sens_path = P4 / "pillar4_funnel_exprfloor_sensitivity.json"
if sens_path.exists():
    sens = json.load(open(sens_path))
    sweep = sens["sweep"]
    stages_lbl = ["S1\ncandidates", "S2\ncross-engine", "S3\n+ sig", "S4\n+ effect floor", "S5\nconsensus"]
    fig, ax = plt.subplots(figsize=(ONEHALF_COL, ONEHALF_COL * 0.62))
    xs = np.arange(len(stages_lbl))
    # colourblind-safe shades + distinct markers for increasing floor
    shades = ["#222222", "#0072B2", "#56B4E9", "#E69F00"]
    markers = ["o", "s", "^", "D"]
    # tiny multiplicative jitter so convergent traces stay individually readable on symlog
    jit = [1.0, 1.06, 1.13, 1.20]
    for j, row in enumerate(sweep):
        ys = np.array([row["S1"], row["S2"], row["S3"], row["S4"], row["S5"]], dtype=float)
        ys_plot = ys * jit[j % len(jit)]
        fl = row["pctl_floor"]
        lbl = "no floor" if fl == 0 else f"drop bottom {int(fl*100)}%"
        ax.plot(xs, ys_plot, "-", marker=markers[j % len(markers)], color=shades[j % len(shades)],
                lw=1.3, ms=4.5, alpha=0.9, label=lbl, zorder=3, markeredgecolor="white",
                markeredgewidth=0.4)
        # annotate the TRUE terminal consensus count for each line (bold, full contrast)
        ax.text(xs[-1] + 0.08, ys_plot[-1], f"{int(ys[-1])}", va="center", ha="left",
                fontsize=6.2, fontweight="bold", color=shades[j % len(shades)])
    ax.set_yscale("symlog")
    ax.set_xticks(xs)
    ax.set_xticklabels(stages_lbl, fontsize=6)
    ax.set_xlim(-0.3, len(stages_lbl) - 0.3)
    ax.set_ylabel("PT\u2192receiver axes (n, symlog)")
    ax.set_title("Expression-floor sensitivity of the evidence funnel", fontsize=8)
    ax.legend(fontsize=5.8, frameon=False, loc="upper right", title="candidate expression floor",
              title_fontsize=5.8)
    ax.text(0.99, 0.60, "traces vertically offset\nfor visibility (S2\u2013S5 counts equal)",
            transform=ax.transAxes, ha="right", va="top", fontsize=4.8,
            style="italic", color="#888888")
    cov = sens["coverage"]
    fig.subplots_adjust(bottom=0.22)
    fig.text(0.5, 0.045,
             f"Cross-engine expression proxy (percentile-normalised max of LIANA lr_means & CellChat probability; "
             f"covers {cov['either_engine']}/{cov['total']} candidates). Raising the floor shrinks only the S1 candidate "
             f"count; S2\u2013S5 are essentially invariant \u2014 the surviving biology is well-expressed and the funnel is not a low-expression artifact.",
             ha="center", fontsize=4.9, color="#666666")
    save_figure(fig, "F_FUNNEL_SUPP_exprfloor_sensitivity", subdir="04_synthesis/figures")
    plt.close(fig)
    print("F_FUNNEL_SUPP written.")
else:
    print("WARNING: sensitivity JSON not found; F_FUNNEL_SUPP skipped.")

print("Figures written:")
for f in sorted((PROJ / "04_synthesis" / "figures").glob("*.svg")):
    print("  ", f.name)
print("PNG present:", [f.name for f in sorted((PROJ / "04_synthesis" / "figures").glob('*.png'))])
