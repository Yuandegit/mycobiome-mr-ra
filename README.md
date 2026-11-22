# Gut mycobiome and rheumatoid arthritis: two-sample Mendelian randomization

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24433532.svg)](https://doi.org/10.5281/zenodo.24433532)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Reproducible code and derived data accompanying the manuscript
**"Mendelian randomization finds no robust causal evidence linking the gut
mycobiome to rheumatoid arthritis and highlights HLA-region pleiotropy."**

The study uses two-sample summary-data Mendelian randomization (MR) to test
whether genetically predicted variation in 115 gut fungal (mycobiome) genera
causally affects rheumatoid arthritis (RA) risk, and it documents why the
nominal signals are better explained by HLA-region horizontal pleiotropy than
by a fungal causal effect.

- **OSF preregistration:** https://osf.io/rfpk2
- **Archived data package (OSF storage):** `mycobiome_mr_registration_package.zip`
- **Archived snapshot DOI (Zenodo):** [10.5281/zenodo.24433532](https://doi.org/10.5281/zenodo.24433532) (release v1.0.0)

## Repository layout

| Path | Contents |
|---|---|
| `preregistration/` | OSF preregistration, analysis plan and decision rules (EN + CN) |
| `scripts/` | All analysis and figure scripts (Python 3) and the shared `mr_core.py` |
| `results/` | Harmonised per-SNP data (`harmonized_forward.pkl`), per-feature MR estimates, sensitivity, calibration and robust-estimator outputs |
| `figures/` | Main figures (PNG and PDF) |
| `manuscript/` | Methods, results and table chapter files |
| `revision/` | Point-by-point response to reviewers and revised sections |
| `data/` | Trait mapping, phenotype list and genus phenocodes |
| `MANIFEST.csv` | File list with size and SHA-256 checksums |

## Reproducing the analysis

```bash
# 1. Create an environment and install dependencies (Python 3.10+)
python -m venv .venv
# Windows:
.venv\Scripts\activate
# macOS/Linux:
# source .venv/bin/activate
pip install -r requirements.txt

# 2. Run the numbered pipeline from the repository root.
#    results/harmonized_forward.pkl is the direct input for every reported
#    estimate, so the restricted source GWAS files are not required.
python scripts/08_forward_mr.py
python scripts/09_sensitivity_hla.py
python scripts/11_steiger.py
python scripts/12_mr_presso.py
python scripts/16_verify_r2_method.py
python scripts/17_build_tables.py
python scripts/18_robust_estimators.py
python scripts/19_robust_estimators_all.py
python scripts/20_winners_curse.py
python scripts/21_orthogonal_hla.py
python scripts/22_transferability.py
python scripts/23_hla_columns.py

# 3. Render the figures
python scripts/fig2_overall.py
python scripts/fig3_hla.py
python scripts/mhc_boundary_mr.py
```

Scripts 03–06 cover downloading and processing the original public source
data (including reading reference VCFs with `pysam`); they are not needed to
reproduce the MR results from the included harmonised data.

## Source data

The original restricted-access source files are not redistributed:

- **Exposure host–mycobiome GWAS (n = 7,350):** Westlake University omics
  portal, DOI [10.15302/vita.2026.04.04.0024](https://doi.org/10.15302/vita.2026.04.04.0024).
- **Outcome BioBank Japan RA (5,348 cases / 173,268 controls):**
  https://humandbs.dbcls.jp under the BBJ access terms.
- **1000 Genomes Phase 3 East-Asian reference (n = 504):** public.

All derived harmonised summary data needed to reproduce every reported value
are included in `results/`, so the original restricted files are not required.

## Ethics and notes

- This is a secondary analysis of de-identified public summary statistics;
  no new human-subjects data were collected and no individual-level data were
  accessed.
- Analyses added during peer review are recorded as documented deviations in
  the preregistration and labelled sensitivity analyses.

## License

Code is released under the [MIT License](LICENSE). The corresponding OSF
registration and its stored materials are dedicated to the public domain under
CC0 1.0.

## Citation

Please cite the archived Zenodo snapshot as:

> Yuan, S. (2026). *Yuandegit/mycobiome-mr-ra: Reproducible two-sample MR of gut mycobiome and rheumatoid arthritis* (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.24433532

The manuscript citation will be added once the paper is accepted.
