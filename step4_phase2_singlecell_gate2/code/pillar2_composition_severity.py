"""Pillar 2 Step 6 — composition/abundance + severity-within-PT + donor EV module score.
All from existing pseudobulk + meta (no re-pull). Module score is per-donor pseudobulk-based
(mean standardized EV_TIER1 CPM) — SUPPORTING/NON-GATING. Per-cell ssGSEA refinement deferred
(non-gating; moot given gate outcome) and noted in memo.

(A) PT abundance by donor: fraction of a donor's cells that are PT, AKI vs normal vs CKD.
    Reported SEPARATELY from DE -> licenses "program (if any) persists after conditioning
    on cell identity" since DE is PT-vs-PT.
(B) Severity-within-PT: per-donor PT EV module score vs frozen injury-severity panel score,
    WITHIN PT (Pillar-1 Analysis D analogue; bulk was rho=0.70).
(C) Donor EV module score: AKI-vs-normal & CKD-vs-normal, Mann-Whitney + Cliff's delta.
"""
import os, json
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

DURABLE = "/mnt/shared-workspace/shared/pillar2"
OUT = "/mnt/results/ev_aki_ckd_project/02_singlecell_specificity/data"
os.makedirs(OUT, exist_ok=True)
GS = json.load(open("/mnt/results/ev_aki_ckd_project/00_preregistration/feasibility/ev_genesets_frozen_v1.3.json"))
EV1 = set(GS["EV_TIER1_PRIMARY"])
INJ = set(json.load(open("/mnt/results/ev_aki_ckd_project/01_bulk_signature/data/injury_panel_frozen.json"))["panel"])
gene_map = pd.read_csv(f"{DURABLE}/gene_map.tsv", sep="\t").set_index("feature_id")["feature_name"]

# ---- load full pseudobulk + meta ----
counts = pd.read_parquet(f"{DURABLE}/pseudobulk_counts.parquet").set_index("sample_id")
meta = pd.read_csv(f"{DURABLE}/pseudobulk_meta.tsv", sep="\t")
if "sample_id" not in meta.columns:
    meta = meta.rename(columns={"index": "sample_id"})
for c in ["disease","dataset_id","suspension_type","donor_id","compartment"]:
    meta[c] = meta[c].astype(str)
meta = meta.set_index("sample_id").loc[counts.index]
dcode = {"normal":"normal","chronic kidney disease":"CKD","acute kidney failure":"AKI"}
meta["disease_code"] = meta["disease"].map(dcode)

# ============ (A) PT abundance by donor ============
# total cells per donor across ALL compartments (from n_cells), and PT cells
# note: only compartments passing >=25-cell floor were emitted as samples, so use ALL
# compartment rows (incl 'other') for a fair denominator.
donor_tot = meta.groupby(["donor_id","disease_code"])["n_cells"].sum().rename("total_cells")
pt_cells = (meta[meta["compartment"]=="PT"].groupby(["donor_id","disease_code"])["n_cells"]
            .sum().rename("pt_cells"))
ab = pd.concat([donor_tot, pt_cells], axis=1).reset_index()
ab["pt_cells"] = ab["pt_cells"].fillna(0)
ab["pt_frac"] = ab["pt_cells"] / ab["total_cells"]
ab.to_csv(f"{OUT}/pt_abundance_by_donor.tsv", sep="\t", index=False)
# test AKI vs normal, CKD vs normal
comp_rows = []
for grp in ["AKI","CKD"]:
    a = ab.loc[ab.disease_code==grp, "pt_frac"].dropna()
    n = ab.loc[ab.disease_code=="normal", "pt_frac"].dropna()
    if len(a)>=3 and len(n)>=3:
        u,p = mannwhitneyu(a, n, alternative="two-sided")
        # Cliff's delta
        delta = (np.sum(a.values[:,None] > n.values) - np.sum(a.values[:,None] < n.values))/(len(a)*len(n))
        comp_rows.append(dict(contrast=f"{grp}_vs_normal", metric="pt_frac",
                              median_grp=float(a.median()), median_normal=float(n.median()),
                              mannwhitney_p=float(p), cliffs_delta=float(delta),
                              n_grp=len(a), n_normal=len(n)))
comp_df = pd.DataFrame(comp_rows)
comp_df.to_csv(f"{OUT}/composition_shift.tsv", sep="\t", index=False)
print("=== (A) PT abundance shift ===")
print(comp_df.to_string(index=False))

# ============ helper: per-donor pseudobulk module score within a compartment ============
def donor_module_score(compartment, gene_set):
    """mean standardized (z across donors) log1p-CPM of gene-set members, per donor."""
    sel = meta[(meta.compartment==compartment) & (meta.disease_code.isin(["AKI","CKD","normal"]))
               & (meta.n_cells>=25)]
    c = counts.loc[sel.index]                     # samples x genes(feature_id)
    # collapse suspensions within donor (sum), then CPM
    c = c.groupby(sel["donor_id"].values).sum()
    cpm = c.div(c.sum(axis=1)+1e-9, axis=0)*1e6
    lcpm = np.log1p(cpm)
    # map to symbols
    lcpm.columns = lcpm.columns.map(gene_map)
    lcpm = lcpm.loc[:, ~lcpm.columns.isna()]
    lcpm = lcpm.T.groupby(level=0).mean().T          # collapse dup symbols
    gs = [g for g in gene_set if g in lcpm.columns]
    sub = lcpm[gs]
    z = (sub - sub.mean(axis=0)) / (sub.std(axis=0)+1e-9)
    score = z.mean(axis=1)
    # donor disease
    dd = sel.drop_duplicates("donor_id").set_index("donor_id")["disease_code"]
    return pd.DataFrame({"donor_id":score.index, "module_score":score.values,
                         "disease_code":dd.loc[score.index].values, "n_genes":len(gs)})

# ============ (C) donor EV module score, PT ============
ev_score = donor_module_score("PT", EV1)
ev_score.to_csv(f"{OUT}/pt_ev_module_score_donor.tsv", sep="\t", index=False)
mod_rows = []
for grp in ["AKI","CKD"]:
    a = ev_score.loc[ev_score.disease_code==grp, "module_score"].dropna()
    n = ev_score.loc[ev_score.disease_code=="normal", "module_score"].dropna()
    u,p = mannwhitneyu(a,n,alternative="two-sided")
    delta = (np.sum(a.values[:,None]>n.values)-np.sum(a.values[:,None]<n.values))/(len(a)*len(n))
    mod_rows.append(dict(contrast=f"{grp}_vs_normal", median_grp=float(a.median()),
                         median_normal=float(n.median()), mannwhitney_p=float(p),
                         cliffs_delta=float(delta), n_grp=len(a), n_normal=len(n)))
mod_df = pd.DataFrame(mod_rows)
mod_df.to_csv(f"{OUT}/pt_ev_module_score_test.tsv", sep="\t", index=False)
print("\n=== (C) PT EV module score (donor-level) ===")
print(mod_df.to_string(index=False))

# ============ (B) severity-within-PT ============
inj_score = donor_module_score("PT", INJ)   # injury severity proxy per donor, within PT
sev = ev_score.merge(inj_score[["donor_id","module_score"]].rename(columns={"module_score":"injury_score"}),
                     on="donor_id")
sev.to_csv(f"{OUT}/severity_within_PT.tsv", sep="\t", index=False)
# correlation across all donors + within AKI
rho_all, p_all = spearmanr(sev["module_score"], sev["injury_score"])
sev_rows = [dict(subset="all_donors", spearman_rho=float(rho_all), spearman_p=float(p_all), n=len(sev))]
for grp in ["AKI","CKD","normal"]:
    s = sev[sev.disease_code==grp]
    if len(s)>=5:
        r,p = spearmanr(s["module_score"], s["injury_score"])
        sev_rows.append(dict(subset=grp, spearman_rho=float(r), spearman_p=float(p), n=len(s)))
sev_corr = pd.DataFrame(sev_rows)
sev_corr.to_csv(f"{OUT}/severity_within_PT_correlation.tsv", sep="\t", index=False)
print("\n=== (B) severity-within-PT: EV module score vs injury-panel score ===")
print(sev_corr.to_string(index=False))
print("\nStep 6 complete ->", OUT)
