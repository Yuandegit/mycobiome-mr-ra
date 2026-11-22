#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build the OSF/Zenodo registration data package:
stage reproducible files, write README and SHA-256 manifest, create a zip."""
import hashlib, zipfile, csv
from pathlib import Path
import shutil

ROOT = Path("/home/user/mycobiome_joint_mr")
PKG_NAME = "mycobiome_mr_registration_package"
STAGE = ROOT / PKG_NAME
OUT_ZIP = ROOT / f"{PKG_NAME}.zip"

# (source relative to ROOT, destination relative to package, include rule)
INCLUDE_DIRS = [
    ("scripts", "scripts"),
    ("figures", "figures"),
    ("manuscript", "manuscript"),
    ("revision", "revision"),
    ("preregistration", "preregistration"),
]
INCLUDE_FILES = [
    ("data/gwas_microbe/fungi_genus/mapping.json", "data/mapping.json"),
    ("data/fungi_phenotypes.tsv", "data/fungi_phenotypes.tsv"),
    ("data/genus_phenocodes.txt", "data/genus_phenocodes.txt"),
]

def keep_file(p: Path):
    if p.suffix in {".pyc", ".log"}:
        return False
    if "__pycache__" in p.parts:
        return False
    return True

if STAGE.exists():
    shutil.rmtree(STAGE)
STAGE.mkdir(parents=True)

# stage directories
for src, dst in INCLUDE_DIRS:
    for f in (ROOT / src).rglob("*"):
        if f.is_file() and keep_file(f):
            rel = f.relative_to(ROOT / src)
            target = STAGE / dst / rel
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(f, target)

# stage results (exclude logs)
for f in (ROOT / "results").rglob("*"):
    if f.is_file() and keep_file(f):
        rel = f.relative_to(ROOT / "results")
        target = STAGE / "results" / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(f, target)

# stage individual data files
for src, dst in INCLUDE_FILES:
    target = STAGE / dst
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(ROOT / src, target)

README = """# Gut mycobiome and rheumatoid arthritis: two-sample Mendelian randomization

Reproducible data package accompanying the preregistration and the manuscript
"Mendelian randomization finds no robust causal evidence linking the gut
mycobiome to rheumatoid arthritis and highlights HLA-region pleiotropy".

## Contents
- `preregistration/`  OSF preregistration (analysis plan and decision rules).
- `scripts/`           All analysis and figure scripts (Python 3).
- `results/`           Harmonised datasets (including `harmonized_forward.pkl`),
                       per-feature MR estimates, sensitivity and calibration outputs.
- `figures/`           Main figures (PNG and PDF).
- `manuscript/`        Method/results/table chapter files.
- `revision/`          Point-by-point response and revised sections.
- `data/`              Trait mapping, phenotype list and phenocodes.
- `MANIFEST.csv`       File list with size and SHA-256 checksums.

## Reproducing the analysis
1. Python 3.10+ with numpy, pandas, scipy, scikit-learn, statsmodels, matplotlib.
2. The harmonised per-SNP data in `results/harmonized_forward.pkl` is the direct
   input for the MR, sensitivity, winner's-curse, transferability and robust-estimator
   scripts. Run the numbered scripts in order (08 -> 23).
3. Figures are produced by `fig2_overall.py`, `fig3_hla.py` and the figure scripts.

## Source data (not redistributed; obtained from public/access-controlled sources)
- Exposure host-mycobiome GWAS (n=7,350): Westlake University omics portal,
  DOI 10.15302/vita.2026.04.0024.
- Outcome BioBank Japan RA (5,348 cases / 173,268 controls):
  https://humandbs.dbcls.jp under the BBJ access terms.
- 1000 Genomes Phase 3 East Asian reference (public).
The derived harmonised summary data needed for every reported value are included,
so the original restricted files are not required to reproduce the MR results.

## Notes
- This is a secondary analysis of de-identified summary statistics; no new
  human-subjects data were collected.
- Analyses added during peer review are recorded as deviations in the
  preregistration (Section 9) and labelled sensitivity analyses.
"""
(STAGE / "README.md").write_text(README, encoding="utf-8")

# manifest with checksums
files = sorted(p for p in STAGE.rglob("*") if p.is_file())
with open(STAGE / "MANIFEST.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["path", "size_bytes", "sha256"])
    for p in files:
        h = hashlib.sha256(p.read_bytes()).hexdigest()
        w.writerow([str(p.relative_to(STAGE)), p.stat().st_size, h])

# zip (rebuild to include the final manifest)
if OUT_ZIP.exists():
    OUT_ZIP.unlink()
with zipfile.ZipFile(OUT_ZIP, "w", zipfile.ZIP_DEFLATED) as z:
    for p in sorted(STAGE.rglob("*")):
        if p.is_file():
            z.write(p, Path(PKG_NAME) / p.relative_to(STAGE))

n = len(list(STAGE.rglob("*")))
print("staged files:", n)
print("package size:", f"{OUT_ZIP.stat().st_size/1e6:.2f} MB")
print("zip:", OUT_ZIP)
