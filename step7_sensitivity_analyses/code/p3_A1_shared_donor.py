"""A1: shared-donor exclusion sensitivity (reviewer M1).

Drops the 10 AKI donors shared between the two atlases (same Census collection)
from each atlas's AKI arm in turn, re-runs BOTH engines on the reduced AKI shard,
and re-evaluates the 3 consensus axes against the FIXED normal side under the
full pre-registered per-dataset pass criteria (direction-both + sig>=1 + floor-both).

Mirrors pillar3_differential_consensus.py logic exactly; reuses pillar3_robustness
engine runners. Reduced arms: a12 AKI 16->6 donors; dea AKI 14->4 donors.
"""
import sys, os, json
import numpy as np
import pandas as pd

sys.path.insert(0, "/mnt/shared-workspace/shared/pillar3")
import pillar3_robustness as pr

ENGINE_DIR = "/workspace/p3_primary"   # primary engine TSVs (normal side fixed)
FORR_ROOT = "/workspace/p3_primary"
OUT = "/workspace/p3_A1"
os.makedirs(OUT, exist_ok=True)

CC_SIG_P = 0.05; CC_FLOOR_ABS = 0.05; CC_FLOOR_REL = 0.20
LIANA_FLOOR_RANKPOS = 50; LIANA_FLOOR_RRA = 0.15

shared = json.load(open("/workspace/shared_donors.json"))["AKI"]["shared"]
print(f"shared AKI donors excluded: {len(shared)}")

cons = pd.read_csv(f"{ENGINE_DIR}/consensus/consensus_landscape_AKIvsnormal.tsv", sep="\t")
keys = list(cons["key"])
print("consensus axes:", keys)

results = {}
for ds in ["a12ccb9b", "dea717d4"]:
    a = pr.load_bundle_adata(f"{FORR_ROOT}/forR_{ds}__AKI")
    a_r = a[~a.obs["donor_id"].astype(str).isin(shared)].copy()
    n_left = a_r.obs["donor_id"].nunique()
    comp_counts = a_r.obs["compartment"].value_counts().to_dict()
    print(f"\n[{ds}] AKI reduced: {a_r.n_obs} cells, {n_left} donors left; compartments {comp_counts}")

    # fixed normal side
    cc_norm = pr.cc_pt_probs(pd.read_csv(f"{ENGINE_DIR}/cellchat_LR_{ds}__normal.tsv", sep="\t"))
    li_norm = pr.liana_pt_positions(pd.read_csv(f"{ENGINE_DIR}/liana_{ds}__normal.tsv", sep="\t"))
    cc_norm_map = cc_norm.set_index("key")["prob"].to_dict()
    li_norm_pos = li_norm.set_index("key")["mag_position"].to_dict()
    li_norm_rra = li_norm.set_index("key")["magnitude_rank"].to_dict()
    worst_norm_pos = (max(li_norm_pos.values()) + 1) if li_norm_pos else 1

    # engines on reduced AKI
    tag = f"{ds}__AKI_noShared"
    lr = pr.run_liana_on_adata(a_r, seed=1)
    lr.to_csv(f"{OUT}/liana_{tag}.tsv", sep="\t", index=False)
    sfr = f"{OUT}/forR_{tag}"
    pr.write_forR_from_adata(a_r, sfr)
    cc = pr.run_cellchat_bundle(sfr, OUT, tag)
    if cc is not None:
        cc.to_csv(f"{OUT}/cellchat_LR_{tag}.tsv", sep="\t", index=False)

    cc_p = pr.cc_pt_probs(cc).set_index("key")
    cc_pv = pr.cc_pt_probs(cc).set_index("key")["pval"].to_dict() if cc is not None else {}
    li_pos = pr.liana_pt_positions(lr).set_index("key")

    axes = {}
    for k in keys:
        prob_aki = float(cc_p["prob"].get(k, 0.0)); prob_nrm = float(cc_norm_map.get(k, 0.0))
        dprob = prob_aki - prob_nrm
        rel = (dprob / prob_nrm) if prob_nrm > 0 else (float("inf") if dprob > 0 else 0.0)
        pval = float(cc_pv.get(k, 1.0))
        cc_up = dprob > 0
        cc_sig = pval < CC_SIG_P
        cc_floor = (dprob >= CC_FLOOR_ABS) or (rel >= CC_FLOOR_REL)
        if k in li_pos.index:
            pos_aki = float(li_pos.loc[k, "mag_position"]); rra_aki = float(li_pos.loc[k, "magnitude_rank"])
        else:
            pos_aki = worst_norm_pos + 1; rra_aki = 1.0
        pos_nrm = float(li_norm_pos.get(k, worst_norm_pos)); rra_nrm = float(li_norm_rra.get(k, 1.0))
        rank_impr = pos_nrm - pos_aki; rra_impr = rra_nrm - rra_aki
        li_up = (rank_impr > 0) and (rra_impr > 0)
        li_floor = (rank_impr >= LIANA_FLOOR_RANKPOS) or (rra_impr >= LIANA_FLOOR_RRA)
        li_sig = li_floor
        passes = bool(cc_up and li_up and (cc_sig or li_sig) and cc_floor and li_floor)
        axes[k] = dict(prob_AKI=prob_aki, prob_normal=prob_nrm, dprob=dprob, rel_increase=rel,
                       pval_AKI=pval, cc_up=bool(cc_up), cc_sig=bool(cc_sig), cc_floor=bool(cc_floor),
                       rank_improvement=rank_impr, rra_improvement=rra_impr,
                       liana_up=bool(li_up), liana_floor=bool(li_floor),
                       passes_full_criteria=passes)
        print(f"  {k}: dprob={dprob:.4f} rel={rel:.2f} p={pval:.3g} | rankImpr={rank_impr:.0f} rraImpr={rra_impr:.3f} | pass={passes}")
    results[ds] = dict(n_donors_left=int(n_left), n_cells=int(a_r.n_obs), compartments=comp_counts, axes=axes)

verdict = {k: all(results[ds]["axes"][k]["passes_full_criteria"] for ds in results) for k in keys}
out = dict(test="A1_shared_donor_exclusion", shared_donors=shared, per_dataset=results,
           consensus_survives_exclusion=verdict,
           all_axes_survive=all(verdict.values()))
with open(f"{OUT}/A1_result.json", "w") as fh:
    json.dump(out, fh, indent=2)
print("\nA1 verdict (axis -> survives full criteria in BOTH reduced atlases):")
print(json.dumps(verdict, indent=2))
print("A1_DONE")
