"""
Tier-1 expanded figures (A-D) for the pillar-synthesis visualization set.
All numbers pulled from the frozen synthesis_data layer; journal style via plotting_style.
Figures written to 04_synthesis/figures/ as SVG (editable) + PNG (600 dpi).
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Patch
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


# ---------------------------------------------------------------------------
# Figure A: cross-pillar convergence (the scientific-signature figure, "Figure 2")
# ---------------------------------------------------------------------------
def fig_A_convergence(d):
    sub = d["p4_subtraction"]
    g3 = d["p3_gate_summary"]
    p1 = d["p1_null"].set_index(["context", "term"])
    p2 = d["p2_null"]
    p2_pt = p2[(p2["compartment"] == "PT") & (p2["contrast"] == "AKI_vs_normal")].set_index("term")

    # per-pillar row content (all values from frozen files)
    rows = [
        dict(tag="Pillar 1", modality="Bulk transcriptomics",
             test="EV gene-set enrichment vs matched\ngeneric (ribosomal) permutation null",
             chip=[("EV_TIER1 z", p1.loc[("native_GSE139061", "EV_TIER1_PRIMARY"), "z"]),
                   ("ribosomal z", p1.loc[("native_GSE139061", "CTRL_RIBOSOMAL"), "z"])],
             outcome="EV signal sits inside the null;\ngeneric control clears it \u2192 not EV-specific"),
        dict(tag="Pillar 2", modality="Single-cell / single-nucleus",
             test="PT sender EV-module score vs matched\ngeneric (ribosomal) permutation null",
             chip=[("EV_TIER1 z", p2_pt.loc["EV_TIER1_PRIMARY", "z"]),
                   ("ribosomal z", p2_pt.loc["CTRL_RIBOSOMAL", "z"])],
             outcome="EV module elevated, but generic\ncontrol ~4\u00d7 stronger \u2192 not EV-specific"),
        dict(tag="Pillar 3", modality="Communication inference\n(CellChat + LIANA, 2 atlases)",
             test="Cross-engine, cross-atlas reproducible\nconsensus + descriptive EV overlay",
             chip=[("consensus axes", g3["n_consensus_axes"]),
                   ("EV-exclusive", sub["primary_consensus"]["n_ev_compatible"])],
             outcome="3 robust consensus axes;\n0 EV-exclusive (descriptive)"),
        dict(tag="Pillar 4", modality="Constructive synthesis",
             test="Route partition of the reproducible\ncommunication landscape",
             chip=[("survivors", sub["reproducible_population_S4"]["n"]),
                   ("EV-exclusive", sub["reproducible_population_S4"]["n_ev_compatible"])],
             outcome="34/34 survivors soluble + ECM;\n0 EV-exclusive"),
    ]

    fig, ax = plt.subplots(figsize=(ps.DOUBLE_COL, ps.DOUBLE_COL * 0.66))
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

    n = len(rows)
    # reserve top for headers (100->93) and bottom for terminus (0->17); rows live in [20, 91]
    top, bot = 91.0, 20.0
    row_h = (top - bot) / n
    y0 = top
    box_x, box_w = 1.5, 24          # pillar/modality box
    test_x, test_w = 27, 25         # test box
    chip_x, chip_w = 54, 12         # numeric chip
    out_x, out_w = 68.5, 30.0       # outcome box (wider, clear gap from chip, room for text)

    for i, r in enumerate(rows):
        yc = y0 - (i + 0.5) * row_h
        h = row_h * 0.80
        # pillar/modality
        bb = FancyBboxPatch((box_x, yc - h/2), box_w, h, boxstyle="round,pad=0.2,rounding_size=1.2",
                            linewidth=1.0, edgecolor="#333333", facecolor="#F2F2F2")
        ax.add_patch(bb)
        ax.text(box_x + box_w/2, yc + h*0.22, r["tag"], ha="center", va="center",
                fontsize=8, fontweight="bold")
        ax.text(box_x + box_w/2, yc - h*0.22, r["modality"], ha="center", va="center", fontsize=6.3)
        # test
        tb = FancyBboxPatch((test_x, yc - h/2), test_w, h, boxstyle="round,pad=0.2,rounding_size=1.2",
                            linewidth=0.8, edgecolor="#888888", facecolor="white")
        ax.add_patch(tb)
        ax.text(test_x + test_w/2, yc, r["test"], ha="center", va="center", fontsize=6.0)
        # numeric chip (two mini-bars comparing signal vs control z, or counts)
        (lab1, v1), (lab2, v2) = r["chip"]
        cb = FancyBboxPatch((chip_x, yc - h/2), chip_w, h, boxstyle="round,pad=0.2,rounding_size=1.0",
                            linewidth=0.8, edgecolor="#888888", facecolor="#FbFbFb")
        ax.add_patch(cb)
        # normalize the two values to a shared scale within the chip for the mini bars
        vmax = max(abs(v1), abs(v2), 1e-9)
        barw = chip_w * 0.30
        # signal bar (first) colored EV-green, control bar grey; bars grow upward from a common baseline
        base = yc - h*0.30
        for j, (lab, v, col) in enumerate([(lab1, v1, C["ev"]), (lab2, v2, C["ns"])]):
            bx = chip_x + chip_w*0.14 + j * (barw + chip_w*0.12)
            frac = abs(v) / vmax
            bh = h * 0.44 * frac
            ax.add_patch(plt.Rectangle((bx, base), barw, bh, facecolor=col, edgecolor="none"))
            ax.text(bx + barw/2, base + bh + 0.4, f"{v:.1f}" if abs(v) < 100 else f"{v:.0f}",
                    ha="center", va="bottom", fontsize=5.0)
            ax.text(bx + barw/2, base - 0.6, lab.split()[0], ha="center", va="top", fontsize=4.6, color="#555555")
        # outcome
        ob = FancyBboxPatch((out_x, yc - h/2), out_w, h, boxstyle="round,pad=0.2,rounding_size=1.2",
                            linewidth=0.8, edgecolor="#B36A00", facecolor="#FFF3E6")
        ax.add_patch(ob)
        ax.text(out_x + out_w/2, yc, r["outcome"], ha="center", va="center", fontsize=5.6)
        # down arrow to next pillar
        if i < n - 1:
            ax.add_patch(FancyArrowPatch((box_x + box_w/2, yc - h/2 - 0.3),
                                         (box_x + box_w/2, yc - row_h + h/2 + 0.3),
                                         arrowstyle="-|>", mutation_scale=8, linewidth=1.0, color="#555555"))

    # convergence terminus node (sits in the reserved bottom margin, clear of the last row)
    term_y = 6.0
    term_h = 8.0
    term_x, term_w = 22, 56
    tb = FancyBboxPatch((term_x, term_y - term_h/2), term_w, term_h,
                        boxstyle="round,pad=0.3,rounding_size=2",
                        linewidth=1.4, edgecolor=C["aki"], facecolor="#FDECEC")
    ax.add_patch(tb)
    ax.text(50, term_y, "PT communication architecture converges on\nSPP1 (soluble) + extracellular-matrix signaling  (non-EV)",
            ha="center", va="center", fontsize=7.2, fontweight="bold")
    # single clean centered funnel arrow from the bottom of the row stack into the terminus
    last_row_bottom = y0 - n * row_h  # bottom edge of the last row band
    ax.add_patch(FancyArrowPatch((50, last_row_bottom + 0.5),
                                 (50, term_y + term_h/2 + 0.4),
                                 arrowstyle="-|>", mutation_scale=13,
                                 linewidth=1.6, color=C["aki"]))

    # column headers
    for hx, htxt in [(box_x + box_w/2, "Level / modality"), (test_x + test_w/2, "Specificity test"),
                     (chip_x + chip_w/2, "Signal vs\ngeneric control"), (out_x + out_w/2, "Outcome")]:
        ax.text(hx, 99.2, htxt, ha="center", va="top", fontsize=6.5, fontweight="bold", color="#333333")

    ax.set_title("Specificity-guided refinement: convergence across four analysis levels",
                 fontsize=9, fontweight="bold", pad=6)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Figure B: specificity permutation-null panel (P1 + P2)
# ---------------------------------------------------------------------------
def fig_B_specificity_nulls(d):
    p1 = d["p1_null"].copy()
    p2 = d["p2_null"].copy()

    # assemble a tidy frame of (context_label, term, z, p) for EV_TIER1 + CTRL_RIBOSOMAL
    def keep(df, terms):
        return df[df["term"].isin(terms)]
    terms = ["EV_TIER1_PRIMARY", "CTRL_RIBOSOMAL"]
    rows = []
    for _, r in keep(p1, terms).iterrows():
        ctx = "P1 " + ("native bulk" if r["context"].startswith("native") else "transplant meta")
        rows.append(dict(panel=ctx, term=r["term"], z=r["z"], p=r["p_greater"]))
    p2sel = p2[(p2["compartment"] == "PT") & (p2["contrast"].isin(["AKI_vs_normal"]))]
    for _, r in keep(p2sel, terms).iterrows():
        rows.append(dict(panel="P2 PT single-cell\n(AKI vs normal)", term=r["term"], z=r["z"], p=r["p_greater"]))
    import pandas as pd
    D = pd.DataFrame(rows)
    panels = ["P1 native bulk", "P1 transplant meta", "P2 PT single-cell\n(AKI vs normal)"]

    fig, ax = plt.subplots(figsize=(ps.ONEHALF_COL, ps.ONEHALF_COL * 0.66))
    ylabels = []
    yticks = []
    y = 0
    gap = 1.0
    term_style = {"EV_TIER1_PRIMARY": dict(color=C["ev"], marker="o", label="EV Tier-1 (hypothesis)"),
                  "CTRL_RIBOSOMAL": dict(color=C["ns"], marker="s", label="Ribosomal (generic control)")}
    seen_labels = set()
    for panel in panels:
        sub = D[D["panel"] == panel]
        for term in ["EV_TIER1_PRIMARY", "CTRL_RIBOSOMAL"]:
            r = sub[sub["term"] == term]
            if r.empty:
                y += 1; continue
            zz = float(r["z"].iloc[0]); pp = float(r["p"].iloc[0])
            st = term_style[term]
            lab = st["label"] if st["label"] not in seen_labels else None
            seen_labels.add(st["label"])
            ax.plot([0, zz], [y, y], color=st["color"], linewidth=1.6, zorder=2, alpha=0.5)
            ax.scatter([zz], [y], color=st["color"], marker=st["marker"], s=42, zorder=3,
                       edgecolor="white", linewidth=0.5, label=lab)
            sig = "*" if pp < 0.05 else "n.s."
            ax.text(zz + 0.35, y, f"z={zz:.2f}  {sig}", va="center", ha="left", fontsize=5.6)
            ylabels.append(term.replace("_", " ").replace("EV TIER1 PRIMARY", "EV Tier-1").replace("CTRL RIBOSOMAL", "Ribosomal ctrl"))
            yticks.append(y)
            y += 1
        y += gap  # gap between panels

    # significance threshold band (|z| for p~0.05 one-sided ~1.645)
    ax.axvspan(-1.645, 1.645, color="#DDDDDD", alpha=0.55, zorder=0, label="null band (|z|<1.645)")
    ax.axvline(0, color="#666666", linewidth=0.6, zorder=1)

    ax.set_yticks(yticks); ax.set_yticklabels(ylabels, fontsize=5.6)
    ax.invert_yaxis()
    ax.set_xlabel("Enrichment z vs permutation null (higher = more enriched)")
    # extra right-hand room so the far-right z-label (ribosomal z=13.23) is never crowded
    ax.set_xlim(-2.5, D["z"].max() * 1.32)
    # add panel bracket labels on right
    ax.set_title("EV signal vs matched generic control across permutation nulls\n(specificity failure made visual)",
                 fontsize=8, fontweight="bold")
    # dedup legend and place it OUTSIDE the plot area (below the axis) to avoid
    # colliding with the z=13.23 data label in the lower-right of the panel
    h, l = ax.get_legend_handles_labels()
    uniq = dict(zip(l, h))
    ax.legend(uniq.values(), uniq.keys(), loc="upper center",
              bbox_to_anchor=(0.5, -0.16), ncol=3, fontsize=5.6,
              frameon=False, handletextpad=0.4, columnspacing=1.2)
    # panel group annotations
    y = 0
    for panel in panels:
        n_here = 2
        ymid = y + (n_here - 1) / 2
        ax.text(-2.55, ymid, panel, ha="right", va="center", fontsize=6.0, fontweight="bold",
                color="#333333", rotation=0)
        y += n_here + gap
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Figure C: per-donor EV module vs injury score (P2)
# ---------------------------------------------------------------------------
def fig_C_ev_vs_injury(d):
    s = d["p2_sev"].copy()
    corr = d["p2_sev_corr"].set_index("subset")
    disease_col = {"normal": C["ns"], "AKI": C["aki"], "CKD": C["ckd"]}
    fig, ax = plt.subplots(figsize=(ps.SINGLE_COL, ps.SINGLE_COL * 0.85))
    for dis, g in s.groupby("disease_code"):
        ax.scatter(g["injury_score"], g["module_score"], s=16, alpha=0.8,
                   color=disease_col.get(dis, "#000000"), edgecolor="white", linewidth=0.3,
                   label=f"{dis} (n={len(g)})")
    # overall regression line
    x = s["injury_score"].values; y = s["module_score"].values
    b, a = np.polyfit(x, y, 1)
    xs = np.linspace(x.min(), x.max(), 50)
    ax.plot(xs, a + b*xs, color="#333333", linewidth=1.0, linestyle="--", zorder=1)
    rho = corr.loc["all_donors", "spearman_rho"]; pval = corr.loc["all_donors", "spearman_p"]
    ax.text(0.03, 0.97, f"Spearman $\\rho$ = {rho:.2f}\n$p$ = {pval:.1e}  (n={int(corr.loc['all_donors','n'])})",
            transform=ax.transAxes, va="top", ha="left", fontsize=6.2,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#CCCCCC", linewidth=0.6))
    ax.set_xlabel("Per-donor PT injury score")
    ax.set_ylabel("Per-donor PT EV-module score")
    ax.set_title("EV-module activity tracks PT injury severity\n(not an EV-specific axis)",
                 fontsize=8, fontweight="bold")
    ax.legend(loc="lower right", fontsize=5.6)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Figure D: consensus-axis robustness multipanel (P3)
# ---------------------------------------------------------------------------
def fig_D_robustness(d):
    axes_order = sdata.CONSENSUS_AXES  # (key, label, receiver, pathway)
    labels = [a[1] for a in axes_order]
    short = [l.replace("_", "\n") for l in labels]
    lodo = d["p3_lodo"]["verdict"]
    cap = d["p3_cap"]["axes"]
    ckd = {a["key"]: a for a in d["p3_ckd"]["per_axis"]}
    scr = d["p3_scramble"]

    capdisp = {
        "receiver_endothelial|SPP1|ITGA5_ITGB1": "SPP1->ITGA5_ITGB1 (endothelial)",
        "receiver_stromal|SPP1|ITGAV_ITGB1": "SPP1->ITGAV_ITGB1 (stromal)",
        "receiver_stromal|COL4A1|ITGA1_ITGB1": "COL4A1->ITGA1_ITGB1 (stromal)",
    }
    x = np.arange(len(labels))
    wid = 0.38
    ds_colors = {"a12ccb9b": C["ckd"], "dea717d4": C["aki"]}
    ds_names = {"a12ccb9b": "snRNA (a12)", "dea717d4": "scRNA (dea)"}

    fig, axs = plt.subplots(2, 2, figsize=(ps.DOUBLE_COL, ps.DOUBLE_COL * 0.62))

    # Panel 1: LODO frac_updir
    ax = axs[0, 0]
    for j, ds in enumerate(["a12ccb9b", "dea717d4"]):
        vals = [lodo[k]["per_dataset_frac"][ds] for k, *_ in axes_order]
        ax.bar(x + (j - 0.5)*wid, vals, wid, color=ds_colors[ds], label=ds_names[ds])
    ax.axhline(0.5, color="#888888", linestyle=":", linewidth=0.8)
    ax.set_ylim(0, 1.05); ax.set_ylabel("frac. donor-drop iters up")
    ax.set_xticks(x); ax.set_xticklabels(short, fontsize=5.2)
    ax.set_title("A  Leave-one-donor-out stability", fontsize=7.5, fontweight="bold", loc="left")
    ax.legend(fontsize=5.2, loc="lower right")

    # Panel 2: cap-4000 CellChat relative increase (dprob/prob_normal, fractional) per dataset, vs +20% floor.
    # cc_rel_increase is the FRACTIONAL increase: 0.1904 -> +19.04%, 1.0781 -> +107.8%.
    # COL4A1 in dea717d4 has prob_normal=0 => 'inf' (string); cap its bar and mark with an infinity glyph.
    ax = axs[0, 1]
    floor_pct = 20.0  # a-priori rel_increase floor = +20%
    def _relpct(v):
        try:
            f = float(v)
            return (np.inf if not np.isfinite(f) else f * 100.0)
        except (TypeError, ValueError):
            return np.inf
    all_finite = []
    for k, *_ in axes_order:
        for ds in ["a12ccb9b", "dea717d4"]:
            v = _relpct(cap[capdisp[k]]["per_dataset"][ds]["cc_rel_increase"])
            if np.isfinite(v):
                all_finite.append(v)
    cap_top = max(all_finite) * 1.15  # display ceiling for inf bars
    inf_drawn = False
    for j, ds in enumerate(["a12ccb9b", "dea717d4"]):
        raw = [_relpct(cap[capdisp[k]]["per_dataset"][ds]["cc_rel_increase"]) for k, *_ in axes_order]
        # finite values as solid bars; infinite entries (prob. in normal = 0) get no solid bar
        fin_vals = [v if np.isfinite(v) else 0.0 for v in raw]
        ax.bar(x + (j - 0.5)*wid, fin_vals, wid, color=ds_colors[ds], label=ds_names[ds])
        for i, v in enumerate(raw):
            if not np.isfinite(v):
                # open hatched bar to the display ceiling + infinity glyph: reads as
                # "unbounded", never as a finite value
                ax.bar(x[i] + (j - 0.5)*wid, cap_top, wid, facecolor="white",
                       edgecolor=ds_colors[ds], hatch="////", linewidth=1.0, zorder=3)
                ax.text(x[i] + (j - 0.5)*wid, cap_top, "$\\infty$", ha="center", va="bottom",
                        fontsize=7, fontweight="bold", color=ds_colors[ds])
                ax.text(x[i] + (j - 0.5)*wid, cap_top * 0.45, "prob. in\nnormal = 0",
                        ha="center", va="center", fontsize=4.4, color=ds_colors[ds], zorder=4,
                        bbox=dict(boxstyle="round,pad=0.15", facecolor="white",
                                  edgecolor="none", alpha=0.85))
                inf_drawn = True
    ax.axhline(floor_pct, color=C["up"], linestyle="--", linewidth=1.0, label="a-priori floor (+20%)")
    ax.set_ylabel("CellChat relative increase (%)")
    ax.set_xticks(x); ax.set_xticklabels(short, fontsize=5.2)
    # generous headroom so the infinity glyph (atop the rightmost/tallest bar) clears the
    # top spine AND the legend (moved to upper-left over the shorter left-hand bars)
    ax.set_ylim(0, cap_top * 1.32)
    ax.set_title("B  Cap-4000 effect size vs floor", fontsize=7.5, fontweight="bold", loc="left")
    handles, labels_ = ax.get_legend_handles_labels()
    if inf_drawn:
        handles.append(Patch(facecolor="white", edgecolor="#666666", hatch="////",
                             linewidth=1.0, label="\u221e (prob. in normal = 0)"))
    ax.legend(handles=handles, fontsize=5.0, loc="upper left", framealpha=0.9)
    # annotate any axis/atlas that drops below the floor (COL4A1 in snRNA a12)
    for i, (k, *_ ) in enumerate(axes_order):
        v = _relpct(cap[capdisp[k]]["per_dataset"]["a12ccb9b"]["cc_rel_increase"])
        if not cap[capdisp[k]]["per_dataset"]["a12ccb9b"]["passes_in_dataset"] and np.isfinite(v):
            ax.annotate("below floor\n(snRNA)", xy=(i - 0.5*wid, v), xytext=(i - 0.5*wid, floor_pct + cap_top*0.12),
                        fontsize=4.8, ha="center", color=C["up"],
                        arrowprops=dict(arrowstyle="->", color=C["up"], linewidth=0.7))

    # Panel 3: scramble negative control (consensus count collapses)
    ax = axs[1, 0]
    seeds = scr["seeds"]
    labels3 = ["primary"] + [f"scramble\nseed {s['seed']}" for s in seeds]
    vals3 = [seeds[0]["n_primary_consensus"]] + [s["n_consensus_scrambled"] for s in seeds]
    cols3 = [C["ev"]] + [C["ns"]] * len(seeds)
    ax.bar(np.arange(len(vals3)), vals3, color=cols3, width=0.6)
    for i, v in enumerate(vals3):
        ax.text(i, v + 0.06, str(v), ha="center", va="bottom", fontsize=6.5, fontweight="bold")
    ax.set_xticks(np.arange(len(labels3))); ax.set_xticklabels(labels3, fontsize=5.6)
    ax.set_ylabel("# consensus axes"); ax.set_ylim(0, 3.6)
    ax.set_title("C  Cell-label scramble control", fontsize=7.5, fontweight="bold", loc="left")

    # Panel 4: CKD persistence
    ax = axs[1, 1]
    cc_up = [ckd[k]["n_datasets_cc_up"] for k, *_ in axes_order]
    li_up = [ckd[k]["n_datasets_liana_up"] for k, *_ in axes_order]
    ax.bar(x - 0.5*wid, cc_up, wid, color=C["ckd"], label="CellChat up (n datasets)")
    ax.bar(x + 0.5*wid, li_up, wid, color=C["aki"], label="LIANA up (n datasets)")
    ax.set_yticks([0, 1, 2]); ax.set_ylabel("# atlases up in CKD")
    ax.set_xticks(x); ax.set_xticklabels(short, fontsize=5.2)
    ax.set_title("D  Persistence into CKD", fontsize=7.5, fontweight="bold", loc="left")
    ax.legend(fontsize=5.0, loc="upper left")
    for i, (k, *_ ) in enumerate(axes_order):
        if ckd[k]["persists_both_datasets_both_engines"]:
            ax.text(i, 2.05, "persists\nboth engines", ha="center", va="bottom", fontsize=4.8, color=C["up"])

    fig.suptitle("Robustness of the three reproducible consensus axes", fontsize=9, fontweight="bold", y=1.0)
    fig.tight_layout()
    return fig


BUILDERS = {
    "FA_cross_pillar_convergence": fig_A_convergence,
    "FB_specificity_null_panel": fig_B_specificity_nulls,
    "FC_ev_module_vs_injury": fig_C_ev_vs_injury,
    "FD_consensus_axis_robustness": fig_D_robustness,
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
