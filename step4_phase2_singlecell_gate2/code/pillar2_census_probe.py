"""Pillar 2 — probe Census for 10x 3' v3 kidney datasets (AKI/CKD/normal).
Metadata-only; identifies datasets + cell counts to plan the full pull and pilot.
"""
import cellxgene_census
import pandas as pd

pd.set_option("display.max_colwidth", 90)
pd.set_option("display.width", 200)

CENSUS_VERSION = "2025-01-30"
DISEASES = ["normal", "chronic kidney disease", "acute kidney failure"]
ASSAY = "10x 3' v3"  # exact Census label

vf = (
    "tissue_general == 'kidney' and is_primary_data == True "
    f"and disease in {DISEASES!r} and assay == {ASSAY!r}"
)

with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
    obs = (
        census["census_data"]["homo_sapiens"]
        .obs.read(
            value_filter=vf,
            column_names=["disease", "donor_id", "cell_type", "dataset_id", "assay", "suspension_type"],
        )
        .concat()
        .to_pandas()
    )

for c in obs.columns:
    obs[c] = obs[c].astype(str)

print("value_filter:", vf)
print("Total 10x 3' v3 cells (3 diseases):", len(obs))
print()
print("=== cells per (disease x dataset_id) ===")
tab = obs.groupby(["disease", "dataset_id"]).size().to_frame("n_cells")
tab["n_donors"] = obs.groupby(["disease", "dataset_id"])["donor_id"].nunique()
print(tab.to_string())
print()
print("=== per-disease donor totals (v3 only) ===")
print(obs.groupby("disease")["donor_id"].nunique())
print()
print("=== suspension_type breakdown (v3) ===")
print(obs.groupby(["disease", "suspension_type"])["donor_id"].nunique())

obs[["disease", "dataset_id", "donor_id", "suspension_type"]].drop_duplicates().to_csv(
    "/workspace/v3_donor_dataset_map.csv", index=False
)
print("\nsaved donor/dataset map -> /workspace/v3_donor_dataset_map.csv")
