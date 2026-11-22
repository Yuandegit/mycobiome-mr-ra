# Response to Reviewers

**Manuscript:** Mendelian randomization finds no robust evidence of a causal effect of the gut mycobiome on rheumatoid arthritis and highlights HLA-region pleiotropy

We thank the reviewers for a rigorous and constructive evaluation. The review correctly identified that the three foundations of the analysis—cross-population transferability, first-stage selection bias, and the rule used to remove outliers—required strengthening. We have now addressed all eight major points and all ten minor points. The headline conclusion is unchanged, but it is now supported by estimators and calibration steps that do not rely on a single, author-defined removal rule. All new analyses are reproducible from the archived code and the public summary statistics.

For clarity, reviewer comments are in *italics* and our responses follow.

---

## Major comments

### M1. Cross-population transferability was not tested.

*The exposure GWAS is Han Chinese and the outcome is Japanese; the portability of the instruments across these populations was not examined.*

A Japanese mycobiome GWAS does not exist, so a direct replication of exposure effects is impossible. We therefore assessed transferability indirectly, using the alignment of instrument allele frequencies between the two populations, which is the principal determinant of whether a selected instrument remains polymorphic and informative in the outcome sample. Frequencies were already allele-aligned at harmonisation.

Across 1,953 instruments with frequency data in both samples, the exposure MAF and outcome EAF were strongly correlated (Pearson r = 0.95; Spearman ρ = 0.84; r = 0.87 on the logit scale). The median absolute frequency difference was 0.015, the 90th percentile 0.067, and only 3.7% of instruments differed by more than 0.10. No instrument was monomorphic in the Japanese sample (minimum EAF = 4.0 × 10⁻⁴). The low-frequency subset (exposure MAF < 0.05; n = 1,198) showed a weaker Pearson correlation (r = 0.44), as expected for rare variants whose small absolute frequencies fluctuate under drift, but the absolute frequency difference remained small (median 0.010). The driver variant itself was polymorphic in both samples (MAF 0.021 in Han; EAF 0.063 in Japanese).

We have added a dedicated "Cross-population transferability" subsection to the Methods and report these results in the Results. We also state the residual limitation explicitly: frequency concordance supports portability of the instruments but cannot prove that exposure effect magnitudes are identical across populations, and low-frequency variants carry the largest relative uncertainty. These results are in `results/transferability.csv`.

### M2. First-stage winner's curse was named but not quantified, and it contaminates the reported R² and F.

*Instruments were selected on the same exposure associations used to compute their strength, so R² and F are inflated. Provide calibrated estimates and label the raw values as upper bounds.*

We agree that the in-sample R² and F are optimistic. Because individual-level and per-cohort data are not available, a true cross-validation or cohort split cannot be performed; we instead implemented two complementary, well-established statistical calibrations and now report the raw values explicitly as upper bounds.

For every one of the 115 features we report total R² under three calibrations: (1) the raw in-sample value (get_r_from_pn), labelled as an **optimistic upper bound**; (2) an empirical-Bayes (EB) linear shrinkage, using a truncation-aware estimate of the prior variance of the true non-centrality (τ² = 0.337; shrinkage factor 0.252); and (3) a two-sided truncated conditional MLE that is conservative for each variant.

The discrepancy is substantial. Median total R² across features falls from 5.23% (raw; maximum 20.1%) to 0.33% under EB and 0.87% under the conditional MLE. For feature 109 specifically, total R² falls from 6.74% (raw) to 0.43% (EB) and 1.65% (conditional MLE). The most important correction concerns the F statistics. The raw joint F has a median of 23.0, and zero traits fall below 10—but this is a mechanical consequence of the P < 10⁻⁵ selection threshold and carries no independent information. After calibration, the median joint F is 1.39 (EB) and 4.01 (conditional MLE), and **all 115 features fall below the conventional weak-instrument threshold of 10**. Feature 109 drops from a raw joint F of 24.1 to 1.44 (EB) and 5.57 (conditional MLE).

We have replaced the single-SNP F framing with these calibrated joint F statistics, labelled the raw R² as an upper bound throughout, and added a "First-stage calibration and winner's curse" subsection to the Methods. Regarding conditional F: the instruments were aggressively clumped (r² < 0.001 in the EAS reference) and are therefore approximately independent, so the joint and per-variant conditional F are close; we state this approximation rather than report an unavailable individual-data conditional regression. Results are in `results/winners_curse_calibration.csv`.

### M3. The "signal returns to null" step relied on a custom outlier-removal rule and must be aligned with established robust estimators.

*Add MR-Lasso, contamination mixture, IWTLS and mode-based estimators, and report the residual-threshold trajectory at 2.5/3/3.5/4.*

We have added all four estimators and applied them to every feature, and we now report the full single-pass residual-threshold trajectory for the primary feature rather than a single threshold.

For feature 109, the estimators separate into two clear groups. Estimators that allow pleiotropic instruments to be identified and set aside return estimates centred on the null: the contamination mixture gives β = −0.007 (SE 0.073; 95% likelihood interval −0.080 to 0.065, crossing zero; mixing proportion ψ = 0.905), and MR-Lasso flags two pleiotropic variants (rs116921468 and the driver rs111394117) and returns β = −0.018 (SE 0.040, P = 0.66) from the remaining 20 instruments. The weighted mode gives β = −0.243 with a wide, bootstrap-based SE of 0.121 (P = 0.046), which is not a stable confirmation. In contrast, estimators that retain every variant and only down-weight it remain significant because the driver is retained: weighted TLS gives β = −0.145 (SE 0.015), and the existing RAPS gives β = −0.129 (P = 9.8 × 10⁻¹⁰). Their significance is therefore attributable to the same single driver rather than to a consensus across methods.

The residual-threshold trajectory makes the threshold dependence explicit. At thresholds of 2.5 and 3.0 the driver variant is removed and the IVW estimate is null (β = −0.015, P = 0.79; β = −0.018, P = 0.66), whereas at 3.5 and 4.0 the driver is retained (its residual magnitude is 3.02) and the estimate stays strongly significant (β = −0.178, P = 5.7 × 10⁻¹⁷). We also note that under the naive IVW fit the largest outcome-scaled residual belongs to rs116921468 (|r| = 4.27), with the driver second (|r| = 3.02).

A technical note on MR-Lasso: the saturated per-SNP-intercept design [bx, Iₙ] makes cross-validation and BIC model selection unstable at small n (they oscillate between zero and n pleiotropic intercepts). We therefore solve the Lasso along the exact LARS path and evaluate at the step where the causal term enters; variants with non-zero intercepts at that deterministic point are flagged as pleiotropic. This is documented in the Methods.

Across all 115 features, the contamination mixture returns **zero** nominally significant associations. We have added a "Pleiotropy-robust estimators" Methods subsection and a corresponding Results subsection, and we now base the primary conclusion on the threshold-independent leave-one-out and pre-specified regional exclusions—which independently return to null—rather than on the single-pass removal rule. Results are in `results/robust_estimators_all.csv`.

### M4. The "corrected IVW" in the positive control should not be used as evidence that the pipeline is trustworthy.

*The BMI→T2D example removes 43% of variants, uses outcome-side information, and mixes in direct glycaemic loci. Reframe it as a failure-mode demonstration.*

We agree and have reframed the section. It is now explicitly titled and presented as a demonstration of failure modes and the feasibility of residual-based diagnosis, not as proof that the pipeline recovers known effects. We state the three limitations directly: (a) removing 42 of 98 variants (43%) exceeds the "fewer than half invalid" boundary required by the weighted median, so the validity of the residual instruments cannot be assumed; (b) because the residual uses outcome-side information (βᵧ − ββₓ)/SEᵧ, the rule is an outcome-driven diagnostic rather than a strict z statistic and is not used to validate the primary analysis; and (c) the instrument set mixed established direct glycaemic loci (e.g. TCF7L2, GIPR) with BMI primary-effect variants, reflecting an instrument-definition issue that is not symmetric to the fungal HLA case. The threshold trajectory (2.5/3.0/3.5/4.0 → β = 0.174/0.184/0.267/0.300) is reported to show this sensitivity rather than to endorse a single corrected number.

### M5. The "pre-specified" boundary sensitivity needs evidence, and the BH family definition must be clarified.

*No pre-registration or timestamped plan is shown, and the BH correction family is unclear.*

We cannot provide an OSF pre-registration retroactively, so we have removed every claim of "pre-specified" and now describe the HLA windows (25–35 Mb and 28.5–33.5 Mb) as a two-window sensitivity analysis, and the P < 10⁻⁵ instrument set as exploratory. We will deposit a timestamped code repository and analysis-plan record before publication and cite it.

We have also unified the false-discovery-rate accounting. The Benjamini–Hochberg correction is applied within two pre-defined families that match the two measurement scales—the 58 continuous CLR features and the 57 binary features—rather than pooling all 115, and this is now stated in the Methods and applied consistently to the tables. The boundary values (e.g. FDR = 0.73 and 0.40) are reported with their family identified. Table 3 now carries a BH-family identifier so the reader can see the correction scope for every feature.

### M6. The negative control is over-interpreted immunologically; add an orthogonal HLA control.

*Cataract cannot test HLA-specific horizontal pleiotropy. Test the driver on other HLA-driven autoimmune diseases.*

We have added the requested orthogonal control on the same BBJ platform and Japanese population, querying the driver variant rs111394117 on systemic lupus erythematosus (SLE), type 1 diabetes (T1D) and ulcerative colitis (UC). We report the results without overstating them.

The driver shows a concordant, nominally significant effect on SLE (β = −0.34, SE 0.16, P = 0.034), in the same direction as its RA effect (β = −0.414, P = 5.3 × 10⁻²²). On T1D and UC the point estimates are close to the null and not significant (β = −0.042, P = 0.63; β = +0.011, P = 0.95). We interpret this honestly: the concordant SLE association demonstrates that the variant is not RA- or mycobiome-specific but carries cross-autoimmune HLA pleiotropy, which is the key evidence needed to upgrade the claim from "a generic technical association is excluded" to "an HLA-region autoimmune pleiotropic effect is localised." The null T1D/UC estimates do not overturn this; those BBJ cohorts have few cases (1,219 and 314) and limited power, and autoimmune diseases differ in their HLA architecture. We retain the cataract control for its original, narrower purpose—ruling out global stratification and a universal technical association—and no longer claim it addresses HLA-specific pleiotropy. Results are in `results/orthogonal_hla_control.csv`.

### M7. Effect interpretability and operationalisation of the exposure.

*The per-SD CLR meaning is unclear, and the two scales give opposite directions; add detection rate, distribution and a robust scale.*

Individual-level detection rates, zero proportions and abundance distributions are not part of the released summary statistics and cannot be reconstructed, so we now state this as an explicit limitation rather than imply them. We have, however, strengthened the scale argument using the data that are available. The continuous CLR and binary detection instruments are non-overlapping at the driver: rs111394117 enters only the continuous CLR instrument set and is absent from the binary set, so the opposite directions (continuous β = −0.167; binary β = +0.058) cannot be driven by the same variant and are consistent with a CLR-specific artefact for a sparse, low-detection genus. We explain that the CLR coefficient is a per-standard-deviation change in the log-relative-abundance score, not an absolute-abundance effect, and we caution that CLR values for a rarely detected genus are sensitive to a small number of high observations. We call for individual-level replication with detection-stratified and Aitchison-distance analyses in the Discussion, rather than present an unavailable re-analysis.

### M8. The exposure resource and the fixed-effect meta-analysis need elaboration.

*Describe peer-review status, open-access path, batch correction, and the rationale and heterogeneity of the five-cohort fixed-effect meta-analysis.*

We have expanded the exposure description. The 2026 resource (Vita; DOI 10.15302/J-ViT.2026.04.0024) is a curated GWAS release of five Chinese cohorts (n = 4,464/1,126/1,159/435/166; two cohorts below 500) based on ITS sequencing, with the summary statistics and the access route cited. We state that per-cohort summary statistics and between-cohort heterogeneity (I²) are not included in the public release, so the fixed-effect meta-analytic heterogeneity cannot be independently recomputed; this is reported as a limitation, and the two small cohorts are flagged as the most sensitive to platform and DNA-extraction differences. We note that the genus-level instrument selection and the null overall conclusion are robust to the exclusion of the HLA region, which limits the practical impact of any residual batch-driven exposure noise. The reference has been completed with the verifiable DOI, and the archived code repository will be cited for the processing pipeline.

---

## Minor comments

1. **Language consistency and formula rendering.** The manuscript will be submitted entirely in English. All subscripted quantities and the residual formula (βᵧ − β·βₓ)/SEᵧ have been converted to equation objects rather than Unicode subscripts to survive type-setting.
2. **Figure 1 red conclusion box.** The bottom box has been changed to a neutral label ("Sensitivity- and control-driven causal assessment"); the specific SNP and genus conclusions are left to the text.
3. **Figure 2A shared axis and the two thresholds.** Continuous CLR and binary log-odds estimates are now shown in separate, scale-specific panels with within-scale standardisation, and the text explains that Bonferroni (4.3 × 10⁻⁴) is shown only as a conservative reference while inference uses Benjamini–Hochberg FDR, which is better suited to the 115 correlated features.
4. **Bar charts for estimates and Figure 4D.** Point-and-interval (forest-style) plots replace the bar charts in Figures 3D, 4C and 4D, and Figure 4D now reports the binned Poisson/binomial test values rather than relying on visual comparison of bars.
5. **LD reference mismatch.** The text now states explicitly that clumping used the EAS panel (n = 504) while the JPT LD annotation (n = 104) is a separate, smaller reference, matching the figure legend.
6. **Table 3 additions.** Columns have been added for the driver's distance to the HLA boundary, whether each feature contains a chromosome-6 (25–35 Mb) instrument, and the BH family identifier.
7. **MDE wording.** The minimum detectable effect is now reported by scale—binary median 0.062 versus continuous median 0.109—so the panel is not described as having a uniform 0.06 resolution.
8. **"Bootstrap weighted median."** The term is now "weighted median (bootstrap standard error)" to avoid implying a bootstrap-derived estimator.
9. **Reference 13 and reproducibility.** The exposure reference now carries the full verifiable DOI, and the data-processing and analysis code will be deposited in a timestamped public repository with its DOI, replacing the "archived with the project" wording.
10. **Causal-direction contribution.** A sentence now states that, combined with the infeasibility of reverse MR (no instrument overlap outside the HLA region) and the scale-dependent Steiger results, the study makes no claim about reverse causation; the null result concerns the forward mycobiome-to-RA direction only.

---

## Summary of the revised argument

The added analyses tighten rather than redirect the conclusion. After winner's-curse calibration all instrument sets are weak (joint F < 10 in 115/115 features); the contamination mixture returns no significant association across the panel; MR-Lasso and the leave-one-out and regional analyses independently return the only FDR-significant feature to the null; and the driver's concordant effect on SLE locates the artefact to shared HLA autoimmune pleiotropy. The contribution is therefore presented as a reusable, fully documented anatomy of how a single low-frequency HLA variant can generate an FDR-significant false-positive in a Mendelian-randomization pipeline applied to sparse, zero-inflated microbiome features, while small true effects remain statistically unexcluded.
