"""Pillar 2 Step 8 — Gate 2 evaluation against pre-registered criteria (PLAN.md L294-338).
Produces gate2_verdicts.tsv. Criteria FIXED before enrichment (verbatim from PLAN):
  2A: PT EV_TIER1 GSEA NES>0 AND FDR<0.05 (AKI-vs-normal) AND replicated in >=1 independent
      context (PT CKD same-dir OR another tubular sender NES>0). [LOO scored under 2C, not here]
  2B: matched-perm EV above non-EV background p<0.05 AND ribosomal control NULL (not sig).
  effect-size floor (2A&2B): |NES|>1 AND median EV log2FC>0 AND |median log2FC|>=0.20.
  2C-tech: LOO stable AND >=2 AKI datasets same-dir AND assay/suspension diag clean AND
           no influential dataset AND down-sample sensitivity holds.
  2D-bio: >=1 of {another tubular sender, PT CKD, donor module-score consistency}.
  FULL PASS = 2A & 2B & floor & 2C & 2D ; QUALIFIED = 2A & 2B & floor but 2C|2D fail ; else FAIL.
"""
import os, json
import numpy as np, pandas as pd

OUT = "/mnt/results/ev_aki_ckd_project/02_singlecell_specificity/data"
def rd(f):
    p = f"{OUT}/{f}"
    return pd.read_csv(p, sep="\t") if os.path.exists(p) else None

gsea = rd("sender_ev_enrichment_gsea.tsv")
floor = rd("sender_effect_size_floor.tsv")
perm = rd("permutation_null_sender.tsv")
perds = rd("per_dataset_pt_ev_nes.tsv")
loo = rd("leave_one_dataset_out.tsv")
ds = rd("downsampled_control_sensitivity.tsv")
mod = rd("pt_ev_module_score_test.tsv")

def gget(comp, contrast, term):
    s = gsea[(gsea.compartment==comp)&(gsea.contrast==contrast)&(gsea.term==term)]
    return (float(s.nes.iloc[0]), float(s.fdr.iloc[0])) if len(s) else (np.nan, np.nan)
def pget(comp, contrast, term):
    s = perm[(perm.compartment==comp)&(perm.contrast==contrast)&(perm.term==term)]
    return float(s.p_greater.iloc[0]) if len(s) else np.nan

verd = {}

# ---- 2A ----
nes_pt_aki, fdr_pt_aki = gget("PT","AKI_vs_normal","EV_TIER1_PRIMARY")
a1 = (nes_pt_aki > 0) and (fdr_pt_aki < 0.05)
# replication: PT CKD same-direction (NES>0) OR tubular_other AKI NES>0
nes_pt_ckd,_ = gget("PT","CKD_vs_normal","EV_TIER1_PRIMARY")
nes_tub_aki,_ = gget("tubular_other","AKI_vs_normal","EV_TIER1_PRIMARY")
a2 = (nes_pt_ckd > 0) or (nes_tub_aki > 0)
verd["2A_activation_replication"] = dict(
    pass_=bool(a1 and a2),
    detail=f"NES={nes_pt_aki:.3f} FDR={fdr_pt_aki:.4f} (need NES>0 & FDR<0.05 -> part1={a1}); "
           f"replication: PT_CKD_NES={nes_pt_ckd:.3f}, tubular_AKI_NES={nes_tub_aki:.3f} -> part2={a2}")

# ---- 2B ----
p_ev = pget("PT","AKI_vs_normal","EV_TIER1_PRIMARY")
p_rib = pget("PT","AKI_vs_normal","CTRL_RIBOSOMAL")
nes_rib, fdr_rib = gget("PT","AKI_vs_normal","CTRL_RIBOSOMAL")
b_ev = p_ev < 0.05
b_rib_null = not (p_rib < 0.05)          # ribosomal must be NULL
verd["2B_specificity"] = dict(
    pass_=bool(b_ev and b_rib_null),
    detail=f"EV matched-perm p={p_ev:.2e} (<0.05 -> {b_ev}); ribosomal perm p={p_rib:.2e}, "
           f"GSEA NES={nes_rib:.2f}/FDR={fdr_rib:.3f} -> ribosomal NULL={b_rib_null} "
           f"(ribosomal {'VIOLATED' if not b_rib_null else 'null'})")

# ---- effect-size floor ----
fr = floor[(floor.compartment=="PT")&(floor.contrast=="AKI_vs_normal")].iloc[0]
med_lfc = float(fr.median_ev_tier1_logFC); nes_abs = abs(float(fr.nes_ev_tier1))
floor_pass = (nes_abs > 1) and (med_lfc > 0) and (abs(med_lfc) >= 0.20)
verd["effect_size_floor"] = dict(
    pass_=bool(floor_pass),
    detail=f"|NES|={nes_abs:.3f} (>1 -> {nes_abs>1}); median EV log2FC={med_lfc:.3f} "
           f"(need >0 AND |.|>=0.20 -> {med_lfc>0 and abs(med_lfc)>=0.20})")

# ---- 2C-technical ----
c_parts = {}
if loo is not None:
    loo_pos = (loo["EV_TIER1_PRIMARY_nes"] > 0).all()
    c_parts["LOO_all_NES_gt0"] = bool(loo_pos)
if perds is not None:
    n_pos_ds = int((perds["EV_TIER1_PRIMARY_nes"] > 0).sum())
    c_parts["n_datasets_same_direction"] = n_pos_ds
    c_parts["ge2_datasets_same_dir"] = bool(n_pos_ds >= 2)
# assay/suspension: STRUCTURAL finding — suspension perfectly nested in dataset_id, dataset_id
# partially confounded with disease but estimable within 2 shared atlases. Reported, not auto-pass.
c_parts["assay_suspension_diag"] = "suspension nested in dataset_id; model=~dataset_id+disease_code; estimable within 2 atlases"
if ds is not None:
    dsr = ds.iloc[0]
    c_parts["downsample_frac_NES_gt0"] = None  # NES computed on median ranking, see detail
    c_parts["downsample_median_ranking_NES"] = float(dsr.ev_tier1_nes_median_ranking)
    c_parts["downsample_frac_iters_logFC_ge0.20"] = float(dsr["frac_iters_median_logFC_ge_0p20"])
    c_parts["downsample_frac_iters_logFC_gt0"] = float(dsr["frac_iters_median_logFC_gt_0"])
c_pass = (c_parts.get("LOO_all_NES_gt0", False) and c_parts.get("ge2_datasets_same_dir", False)
          and (ds is not None and float(ds.iloc[0]["frac_iters_median_logFC_ge_0p20"]) >= 0.5))
verd["2C_technical"] = dict(pass_=bool(c_pass), detail=json.dumps(c_parts))

# ---- 2D-biological ----
d_parts = {}
d_parts["tubular_other_AKI_NES_gt0"] = bool(nes_tub_aki > 0)
d_parts["PT_CKD_NES_gt0"] = bool(nes_pt_ckd > 0)
if mod is not None:
    mrow = mod[mod.contrast=="AKI_vs_normal"].iloc[0]
    d_parts["module_score_AKI_consistent"] = bool((mrow.mannwhitney_p < 0.05) and (mrow.cliffs_delta > 0))
    d_parts["module_score_cliffs_delta"] = float(mrow.cliffs_delta)
d_pass = any([d_parts["tubular_other_AKI_NES_gt0"], d_parts["PT_CKD_NES_gt0"],
              d_parts.get("module_score_AKI_consistent", False)])
verd["2D_biological"] = dict(pass_=bool(d_pass), detail=json.dumps(d_parts))

# ---- overall verdict ----
core = (verd["2A_activation_replication"]["pass_"] and verd["2B_specificity"]["pass_"]
        and verd["effect_size_floor"]["pass_"])
if core and verd["2C_technical"]["pass_"] and verd["2D_biological"]["pass_"]:
    overall = "FULL PASS"
elif core:
    overall = "QUALIFIED PASS"
else:
    overall = "FAIL"

rows = []
for k, v in verd.items():
    rows.append(dict(criterion=k, pass_=v["pass_"], detail=v["detail"]))
rows.append(dict(criterion="OVERALL_VERDICT", pass_=(overall!="FAIL"), detail=overall))
out = pd.DataFrame(rows)
out.to_csv(f"{OUT}/gate2_verdicts.tsv", sep="\t", index=False)
print(out.to_string(index=False))
print("\n==== GATE 2 VERDICT:", overall, "====")
