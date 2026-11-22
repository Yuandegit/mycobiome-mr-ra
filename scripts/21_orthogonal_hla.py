"""
M6 response: orthogonal HLA pleiotropy control. Using the same BBJ platform and
Japanese population, query the driver variant rs111394117 (NOTCH4, chr6:32187213)
on three other HLA-driven autoimmune diseases (SLE, T1D, UC) via the PheWeb
region API, and compare with its RA effect. A consistent effect across diseases
locates the signal to shared HLA autoimmune pleiotropy rather than RA or the
mycobiome. Self-contained Python, no R/PLINK.
"""
import sys, os, json, urllib.parse, urllib.request, pickle
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, pandas as pd
from scipy import stats
from mr_core import ROOT

BASE = "https://pheweb.jp/api/region/{pheno}/lz-results/?filter={flt}"
DRIVER_POS = 32187213
HALF = 1500


def fetch(pheno):
    flt = ("chromosome in  '6' and position ge %d and position le %d"
           % (DRIVER_POS - HALF, DRIVER_POS + HALF))
    url = BASE.format(pheno=pheno, flt=urllib.parse.quote(flt))
    with urllib.request.urlopen(url, timeout=60) as r:
        d = json.load(r)["data"]
    i = next(i for i, p in enumerate(d["position"]) if p == DRIVER_POS)
    b, se = d["beta"][i], d["sebeta"][i]
    p = 2 * stats.norm.sf(abs(b / se)) if se > 0 else np.nan
    return {"pheno": pheno, "rsid": d["rsid"][i], "gene": d["nearest_genes"][i],
            "beta": b, "se": se, "p": p, "af": d["af"][i]}


# RA effect from the local harmonised data (already on the BBJ scale)
H = pickle.load(open(os.path.join(ROOT, "results/harmonized_forward.pkl"), "rb"))
ra = H[("109", "expl1e5")]
ra_drv = ra[ra.rsid == "rs111394117"].iloc[0]

rows = [{"pheno": "RA", "rsid": "rs111394117", "gene": "NOTCH4",
         "beta": ra_drv.by, "se": ra_drv.sy,
         "p": 2 * stats.norm.sf(abs(ra_drv.by / ra_drv.sy)), "af": ra_drv.maf}]
for ph in ["SLE", "T1D", "UC"]:
    rows.append(fetch(ph))

out = pd.DataFrame(rows)
out["direction"] = np.where(out.beta < 0, "negative", "positive")
path = os.path.join(ROOT, "results/orthogonal_hla_control.csv")
out.to_csv(path, index=False)
print(out.to_string(index=False))
print("\nwrote", path)
print("sign concordance with RA:", (np.sign(out.beta) == np.sign(out.iloc[0].beta)).all())
