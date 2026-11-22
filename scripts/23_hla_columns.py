#!/usr/bin/env python3
"""Minor 6: per-feature HLA columns for Supplementary Table S1.
- whether the feature contains a chromosome-6 (25-35 Mb) instrument and how many
- distance of the nearest HLA instrument to the nearest HLA boundary (Mb)
- BH family (scale)
"""
import pickle, json
import pandas as pd

ROOT = "/home/user/mycobiome_joint_mr"
LO, HI = 25_000_000, 35_000_000

harm = pickle.load(open(f"{ROOT}/results/harmonized_forward.pkl", "rb"))
mapping = json.load(open(f"{ROOT}/data/gwas_microbe/fungi_genus/mapping.json"))

rows = []
for pheno in sorted(mapping.keys(), key=lambda x: int(x)):
    df = harm[(pheno, "expl1e5")]
    keep = df[df["keep"]].copy()
    hla = keep[(keep["chrom"] == 6) & (keep["pos"] >= LO) & (keep["pos"] <= HI)]
    n_hla = len(hla)
    if n_hla > 0:
        # distance of every HLA instrument to the nearest boundary, then take the smallest
        d = (hla["pos"] - LO).where(hla["pos"] - LO < HI - hla["pos"], HI - hla["pos"])
        i = d.idxmin()
        nearest_dist_mb = d.loc[i] / 1e6
        nearest_rsid = hla.loc[i, "rsid"]
        nearest_pos = int(hla.loc[i, "pos"])
        has = f"Yes ({n_hla})"
        dist = f"{nearest_dist_mb:.2f}"
        nrsid = nearest_rsid
        npos = nearest_pos
    else:
        has = "No"
        dist = "–"
        nrsid = "–"
        npos = "–"
    scale = mapping[pheno].split("_")[0]
    family = "Continuous (58)" if scale == "Continuous" else "Binary (57)"
    rows.append({
        "pheno": int(pheno),
        "bh_family": family,
        "chr6_25_35_iv": has,
        "n_hla_iv": n_hla,
        "nearest_boundary_dist_mb": dist,
        "nearest_hla_rsid": nrsid,
        "nearest_hla_pos": npos,
    })

out = pd.DataFrame(rows)
out.to_csv(f"{ROOT}/results/hla_table_columns.csv", index=False)
print("rows:", len(out))
print("features with HLA IV:", (out["n_hla_iv"] > 0).sum())
print(out[out["n_hla_iv"] > 0].to_string(index=False))
