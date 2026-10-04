"""
Tier-2 expanded figures (E-H) for the pillar-synthesis visualization set.
All numbers pulled from frozen artifacts (no new inference); journal style via plotting_style.
Figures written to 04_synthesis/figures/ as SVG (editable) + PNG (600 dpi).

E. Receiver x ligand bubble map of the 34 reproducible survivors (route-colored).
F. SPP1 / COLLAGEN pathway-coherence bars (dprob + up-rank per atlas).
G. Transplant-injury meta forest (top concordant genes; injury-axis CONTEXT, not EV-specific).
H. Two-atlas effect concordance scatter (LIANA magnitude-rank improvement; 3 consensus highlighted).
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe

import plotting_style as ps
import synthesis_data as sdata

ps.set_journal_style()
C = ps.COLORS
SUBDIR = "04_synthesis/figures"

ROUTE_COLORS = {
    "EV-compatible": C["ev"],
    "Soluble-classical": C["aki"],
    "ECM": C["ckd"],
    "Contact": "#CC79A7",
    "unclassified": C["ns"],
}
CONSENSUS_KEYS = {
    "receiver_endothelial|SPP1|ITGA5_ITGB1",
    "receiver_stromal|SPP1|ITGAV_ITGB1",
    "receiver_stromal|COL4A1|ITGA1_ITGB1",
}


# ---------------------------------------------------------------------------
# Figure E: receiver x ligand bubble map of the 34 reproducible survivors
# ---------------------------------------------------------------------------
def fig_E_bubble_map(d):
    r = d["p4_routes34"].copy()  # 34 rows: key, receiver, ligand, receptor, route_tag, ev_compatible, is_consensus
    # bubble size = # axes for that (receiver, ligand) pair; color = route
    grp = (r.groupby(["receiver", "ligand", "route_tag"])
             .agg(n_axes=("key", "size"), any_consensus=("is_consensus", "any"))
             .reset_index())

    # order ligands by total axes (desc), receivers by total axes (desc)
    lig_order = (r.groupby("ligand").size().sort_values(ascending=False)).index.tolist()
    rec_order = (r.groupby("receiver").size().sort_values(ascending=False)).index.tolist()
    lig_ix = {l: i for i, l in enumerate(lig_order)}
    rec_ix = {rc: i for i, rc in enumerate(rec_order)}

    fig, ax = plt.subplots(figsize=(ps.ONEHALF_COL, ps.ONEHALF_COL * 0.72))
    for _, g in grp.iterrows():
        x = lig_ix[g["ligand"]]; y = rec_ix[g["receiver"]]
        size = 40 + g["n_axes"] * 80
        ax.scatter(x, y, s=size, color=ROUTE_COLORS.get(g["route_tag"], C["ns"]),
                   alpha=0.85, edgecolor="black", linewidth=(1.4 if g["any_consensus"] else 0.5),
                   zorder=3)
        ax.text(x, y, str(int(g["n_axes"])), ha="center", va="center", fontsize=5.4,
                color="white", fontweight="bold", zorder=4,
                path_effects=[pe.withStroke(linewidth=1.1, foreground="black")])

    ax.set_xticks(range(len(lig_order))); ax.set_xticklabels(lig_order, rotation=45, ha="right", fontsize=6.0)
    ax.set_yticks(range(len(rec_order))); ax.set_yticklabels(rec_order, fontsize=6.5)
    ax.set_xlim(-0.6, len(lig_order) - 0.2); ax.set_ylim(-0.7, len(rec_order) - 0.3)
    ax.set_xlabel("Ligand"); ax.set_ylabel("Receiver compartment")
    ax.set_title("Supporting reproducible communication landscape\n(34 survivors; bubble = # axes, thick ring = consensus)",
                 fontsize=8, fontweight="bold")
    ax.grid(True, linestyle=":", linewidth=0.4, alpha=0.5, zorder=0)

    # legends: route color + consensus ring + a size guide
    route_handles = [Line2D([0], [0], marker="o", linestyle="none", markersize=6,
                            markerfacecolor=ROUTE_COLORS[k], markeredgecolor="black",
                            markeredgewidth=0.4, label=k)
                     for k in ["ECM", "Soluble-classical"] if k in set(grp["route_tag"])]
    ring = Line2D([0], [0], marker="o", linestyle="none", markersize=7, markerfacecolor="#DDDDDD",
                  markeredgecolor="black", markeredgewidth=1.4, label="contains consensus axis")
    leg1 = ax.legend(handles=route_handles + [ring], loc="upper right", fontsize=5.6,
                     title="route", title_fontsize=5.8, framealpha=0.9)
    ax.add_artist(leg1)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Figure F: SPP1 / COLLAGEN pathway coherence across the two atlases
# ---------------------------------------------------------------------------
def fig_F_pathway_coherence(d):
    pc = d["p3_pathway"]["consensus_pathways"]
    pathways = ["SPP1", "COLLAGEN"]
    ds_order = ["a12ccb9b", "dea717d4"]
    ds_names = {"a12ccb9b": "snRNA (a12)", "dea717d4": "scRNA (dea)"}
    ds_colors = {"a12ccb9b": C["ckd"], "dea717d4": C["aki"]}

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(ps.ONEHALF_COL, ps.ONEHALF_COL * 0.5))
    x = np.arange(len(pathways)); wid = 0.38

    # Panel 1: aggregate pathway-level dprob (CellChat), up in AKI both atlases
    for j, ds in enumerate(ds_order):
        vals = [pc[pw][ds]["dprob"] for pw in pathways]
        ax1.bar(x + (j - 0.5) * wid, vals, wid, color=ds_colors[ds], label=ds_names[ds])
    ax1.set_xticks(x); ax1.set_xticklabels(pathways, fontsize=6.5)
    ax1.set_ylabel("Pathway $\\Delta$prob (AKI - normal)")
    ax1.set_title("A  Aggregate pathway signal up in AKI", fontsize=7.2, fontweight="bold", loc="left")
    ax1.legend(fontsize=5.4, loc="upper right")
    ax1.axhline(0, color="#666666", linewidth=0.5)

    # Panel 2: rank of the pathway among up-in-AKI pathways (1 = strongest); lower is better.
    # Bars grow DOWNWARD from the top after inversion; place the #rank label just inside the
    # bar top (small positive offset in data units) so it never clips at the y-limit.
    max_rank = max(pc[pw][ds]["rank_up_in_AKI"] for pw in pathways for ds in ds_order)
    for j, ds in enumerate(ds_order):
        vals = [pc[pw][ds]["rank_up_in_AKI"] for pw in pathways]
        ax2.bar(x + (j - 0.5) * wid, vals, wid, color=ds_colors[ds], label=ds_names[ds])
        for i, v in enumerate(vals):
            ax2.text(x[i] + (j - 0.5) * wid, v - 0.18, f"#{int(v)}", ha="center", va="top",
                     fontsize=5.8, color="white", fontweight="bold")
    ax2.set_xticks(x); ax2.set_xticklabels(pathways, fontsize=6.5)
    ax2.set_ylabel("Rank among up-in-AKI pathways")
    ax2.set_ylim(max_rank + 0.6, 0)  # inverted (rank 1 at top) with headroom at both ends
    ax2.set_title("B  Both consensus pathways rank near the top", fontsize=7.2, fontweight="bold", loc="left")
    # no legend in Panel B (Panel A legend already identifies the two atlases)

    fig.suptitle("Pathway-level coherence of the consensus signal (both atlases, both up in AKI)",
                 fontsize=8.5, fontweight="bold", y=1.02)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Figure G: transplant-injury meta forest (injury-axis CONTEXT, not EV-specific)
# ---------------------------------------------------------------------------
def fig_G_meta_forest(d, n_each=8):
    m = d["p1_transplant_meta"].copy()
    up = m[m["pooled_logFC"] > 0].nsmallest(n_each, "adj_pval").copy()
    dn = m[m["pooled_logFC"] < 0].nsmallest(n_each, "adj_pval").copy()
    up["dir"] = "up"; dn["dir"] = "down"
    sel = pd.concat([up.sort_values("pooled_logFC"), dn.sort_values("pooled_logFC")], ignore_index=True)
    sel = sel.sort_values("pooled_logFC").reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(ps.SINGLE_COL, ps.SINGLE_COL * 1.15))
    y = np.arange(len(sel))
    for i, row in sel.iterrows():
        col = C["up"] if row["pooled_logFC"] > 0 else C["down"]
        ax.plot([row["ci_lb"], row["ci_ub"]], [i, i], color=col, linewidth=1.4, zorder=2)
        ax.scatter(row["pooled_logFC"], i, s=22, color=col, zorder=3, edgecolor="white", linewidth=0.4)
    ax.axvline(0, color="#666666", linewidth=0.6, linestyle="-")
    ax.set_yticks(y); ax.set_yticklabels(sel["gene"], fontsize=5.8)
    # widen x so no CI touches the frame; give one row of vertical padding so the last
    # gene row is never crowded by an annotation
    xlo = sel["ci_lb"].min(); xhi = sel["ci_ub"].max()
    xpad = (xhi - xlo) * 0.14
    ax.set_xlim(xlo - xpad, xhi + xpad)
    ax.set_ylim(-1.4, len(sel) - 0.4)
    ax.set_xlabel("Pooled log$_2$FC (random-effects meta, transplant-injury cohorts)")
    ax.set_title("Transplant-injury meta-signature context\n(top concordant genes; injury axis, not EV-specific)",
                 fontsize=7.6, fontweight="bold")
    # I2 note in the reserved bottom-margin strip (below all gene rows), left-aligned
    i2max = sel["I2"].max()
    ax.text(xlo - xpad*0.6, -1.05, f"all shown: $I^2$ = {i2max:.0f}%  (homogeneous pooling)",
            fontsize=5.4, va="center", ha="left",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#CCCCCC", linewidth=0.5))
    # direction legend (upper-right; down-regulated genes are at the bottom-left so this is clear)
    handles = [Line2D([0], [0], marker="o", linestyle="-", color=C["up"], markersize=5, label="up in injury"),
               Line2D([0], [0], marker="o", linestyle="-", color=C["down"], markersize=5, label="down in injury")]
    ax.legend(handles=handles, loc="center right", fontsize=5.4, framealpha=0.95)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Figure H: two-atlas effect concordance (LIANA magnitude-rank improvement)
# ---------------------------------------------------------------------------
def fig_H_concordance(d):
    pa = d["p3_per_axis"].copy()
    # the 44 cross-engine direction-concordant axes (S2): dir_both in >=1 dataset
    s2_keys = pa.loc[pa["dir_both"], "key"].unique().tolist()
    sub = pa[pa["key"].isin(s2_keys)].copy()
    # LIANA magnitude-rank improvement is the cross-engine-complete effect metric
    piv = sub.pivot_table(index="key", columns="dataset", values="rank_improvement", aggfunc="first")
    both = piv.dropna().copy()                       # axes scored by LIANA in both atlases
    n_off = piv.shape[0] - both.shape[0]             # concordant axes scored in only one atlas
    # membership tiers
    passes = pa.loc[pa["passes_in_dataset"]]
    s4_keys = set(passes["key"].unique())
    nds = passes.groupby("key")["dataset"].nunique()
    s5_keys = set(nds[nds >= 2].index)

    fig, ax = plt.subplots(figsize=(ps.SINGLE_COL, ps.SINGLE_COL))
    xa, ya = "a12ccb9b", "dea717d4"
    # background (concordant-only) points
    for k, row in both.iterrows():
        if k in s5_keys:
            continue
        is_s4 = k in s4_keys
        ax.scatter(row[xa], row[ya], s=(26 if is_s4 else 16),
                   color=(C["ns"] if not is_s4 else "#7A7A7A"),
                   alpha=0.8, edgecolor="white", linewidth=0.3, zorder=2)
    # consensus points highlighted, with manually placed labels + leader lines into open
    # regions so the two near-coincident SPP1 points do not collide. Offsets in DATA units.
    cons_labels = {
        "receiver_endothelial|SPP1|ITGA5_ITGB1": "SPP1\u2192ITGA5_ITGB1",
        "receiver_stromal|SPP1|ITGAV_ITGB1": "SPP1\u2192ITGAV_ITGB1",
        "receiver_stromal|COL4A1|ITGA1_ITGB1": "COL4A1\u2192ITGA1_ITGB1",
    }
    # label anchor positions (data coords) chosen to sit in sparse areas of the cloud
    label_xy = {
        "receiver_stromal|COL4A1|ITGA1_ITGB1": (-95, 255),   # upper-left open region
        "receiver_endothelial|SPP1|ITGA5_ITGB1": (210, 55),  # right, mid
        "receiver_stromal|SPP1|ITGAV_ITGB1": (170, -140),    # lower-right open region
    }
    for k in s5_keys:
        if k not in both.index:
            continue
        px, py = both.loc[k, xa], both.loc[k, ya]
        ax.scatter(px, py, s=70, color=C["ev"], edgecolor="black", linewidth=1.0, zorder=5)
        lx, ly = label_xy.get(k, (px + 20, py + 20))
        ha = "left" if lx >= px else "right"
        ax.annotate(cons_labels[k], xy=(px, py), xytext=(lx, ly),
                    fontsize=5.2, color="#0A5C3A", fontweight="bold", zorder=7,
                    ha=ha, va="center",
                    bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                              edgecolor=C["ev"], linewidth=0.6, alpha=0.95),
                    arrowprops=dict(arrowstyle="-", color=C["ev"], linewidth=0.7,
                                    shrinkA=1, shrinkB=3))

    # 1:1 reference and quadrant guides
    lo = min(both[xa].min(), both[ya].min())
    hi = max(both[xa].max(), both[ya].max())
    pad = (hi - lo) * 0.08
    ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], linestyle="--", color="#999999",
            linewidth=0.8, zorder=1, label="1:1 (equal effect)")
    ax.set_xlim(lo - pad, hi + pad); ax.set_ylim(lo - pad, hi + pad)
    ax.set_xlabel("LIANA magnitude-rank improvement\nsnRNA atlas (a12)")
    ax.set_ylabel("LIANA magnitude-rank improvement\nscRNA atlas (dea)")
    ax.set_title("Two-atlas effect concordance\n(3 consensus axes reproduce in both atlases)",
                 fontsize=7.6, fontweight="bold")

    handles = [
        Line2D([0], [0], marker="o", linestyle="none", markersize=7, markerfacecolor=C["ev"],
               markeredgecolor="black", markeredgewidth=0.8, label="consensus (pass both atlases)"),
        Line2D([0], [0], marker="o", linestyle="none", markersize=5, markerfacecolor="#7A7A7A",
               markeredgecolor="white", label="single-atlas full pass"),
        Line2D([0], [0], marker="o", linestyle="none", markersize=4, markerfacecolor=C["ns"],
               markeredgecolor="white", label="concordant only"),
        Line2D([0], [0], linestyle="--", color="#999999", label="1:1 (equal effect)"),
    ]
    ax.legend(handles=handles, loc="lower left", fontsize=5.0, framealpha=0.9)
    ax.text(0.98, 0.98,
            f"n = {both.shape[0]} axes scored by LIANA in both atlases\n"
            f"({n_off} of 44 concordant axes scored in one atlas only, off-panel)",
            transform=ax.transAxes, fontsize=4.8, va="top", ha="right",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#CCCCCC", linewidth=0.5))
    fig.tight_layout()
    return fig


BUILDERS = {
    "FE_receiver_ligand_bubble_map": fig_E_bubble_map,
    "FF_pathway_coherence": fig_F_pathway_coherence,
    "FG_transplant_injury_forest": fig_G_meta_forest,
    "FH_two_atlas_concordance": fig_H_concordance,
}


def main():
    d = sdata.load_all()
    sdata.verify(d)
    written = {}
    for name, fn in BUILDERS.items():
        fig = fn(d)
        paths = ps.save_figure(fig, name, subdir=SUBDIR)
        plt.close(fig)
        written[name] = paths
        print(f"wrote {name}: {paths['svg'].split('/')[-1]} + png")
    return written


if __name__ == "__main__":
    main()
