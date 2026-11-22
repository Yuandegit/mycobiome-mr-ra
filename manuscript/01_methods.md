# Materials and Methods

## Study design and data sources

We performed a two-sample, summary-data Mendelian randomization (MR) analysis to test whether
genetically predicted variation in the gut mycobiome influences the risk of rheumatoid arthritis
(RA). The study was built entirely on publicly available genome-wide association study (GWAS)
summary statistics, so no new participant recruitment or experimental work was involved; both
source studies had obtained appropriate ethical approvals and informed consent, and the data used
here contain no personally identifiable information.

Genetic associations with the gut mycobiome were taken from the 2026 host–mycobiome GWAS
(Vita; DOI 10.15302/vita.2026.04.0024), which meta-analysed five independent Chinese cohorts
(the Tianjin Health Study and Brain Cohort, n=4,464; the Zhoushan Maternal and Child Health
Cohort, n=1,126; the Guangzhou Nutrition and Health Study, n=1,159; the WeBirth cohort, n=435;
and the WePN-DA cohort, n=166) using fixed-effects inverse-variance meta-analysis, giving a
total of 7,350 Han Chinese participants. Fungal abundance was profiled by internal transcribed
spacer (ITS) sequencing. Continuous taxonomic features were transformed by the centred
log-ratio (CLR) to remove compositional dependence and then standardised to unit standard
deviation, and were analysed with a mixed linear model (GCTA-MLMA); presence/absence features
were analysed by logistic regression (PLINK2). The consortium reported a study-wide significance
threshold of P<2.2×10⁻¹⁰. We focused on the genus level and retained 115 genus-level trait
measures, comprising 58 continuous (CLR abundance) and 57 binary (detection) features.

The RA outcome was drawn from the BioBank Japan (BBJ) GWAS, which included 5,348 cases and
173,268 controls (total N=178,616; case fraction 0.0299) of Japanese ancestry. Although the
exposure and outcome samples differ in nationality, both are East Asian, which reduces but does
not eliminate the cross-ancestry transfer problem that is more severe for European–Asian
comparisons. All genomic coordinates
were on the human genome reference GRCh37/hg19. For linkage-disequilibrium (LD) estimation we
used the 1000 Genomes Project Phase 3 East Asian (EAS) panel (n=504), and the same BBJ platform
provided the control datasets described below (type 2 diabetes, cataract and body mass index).

The interpretation of MR rests on three assumptions about the instruments: they are genuinely
associated with the exposure (relevance); they are independent of confounders of the
exposure–outcome relationship (independence); and they influence the outcome only through the
exposure (the exclusion-restriction assumption). The relevance assumption is addressed through
instrument-strength screening. The exclusion-restriction assumption is the most demanding in
this setting, because variants in the HLA region can act on both fungal measurement and
autoimmune disease through separate biological pathways; we therefore did not rely on a single
estimator and built a programme of pleiotropy tests, regional exclusion and external controls to
test this condition directly. Because the exposure and outcome GWAS are derived from
separate participants, sample overlap is not a concern; the two populations nevertheless differ
in nationality, so residual population stratification cannot be ruled out and is examined in the
sensitivity and control analyses below.

## Instrument selection and clumping

For each fungal trait we selected candidate variants at two association thresholds. The
exploratory set used P<1×10⁻⁵, a relaxed threshold deliberately chosen because microbiome GWAS
are comparatively underpowered and few individual fungal variants reach conventional
genome-wide significance; the strict set used the conventional P<5×10⁻⁸. Candidate variants were
pruned to independent instruments using the EAS reference panel rather than a distance-only rule.
We applied greedy clumping, implemented in Python with genotype correlation computed directly
from the EAS reference panel through pysam rather than approximated by physical distance, that
ranks variants by association P value and retains a variant only if it has r²<0.10 with every
previously retained variant within a 500-kb window. The EAS panel combines the Han Chinese
(CHB, CHS), Japanese (JPT), Vietnamese (KHV) and Dai (CDX) samples, so it spans both ancestral
backgrounds represented in the exposure and outcome. To verify that inferences did not depend on
this choice, we repeated clumping under three alternative settings (r²<0.01 within 1,000 kb;
r²<0.001 within 1,000 kb; and r²<0.10 within 250 kb).

At the main clumping parameters, the exploratory threshold produced 3,230 instrument–phenotype
associations spanning all 115 traits (median 25 per trait, minimum 10), and the strict threshold
produced 67 associations covering 39 traits. These counts refer to clumped instrument–phenotype
pairs and map to 3,009 distinct genomic positions.

## Harmonisation of exposure and outcome data

Clumped instruments were aligned to the BBJ summary file by chromosome and position. When the
exposure effect allele matched the outcome allele 2, effect estimates were used directly. When
the alleles were reversed relative to the outcome file, the outcome log-odds ratio was negated
and its effect-allele frequency converted to 1 minus that frequency. Palindromic (A/T or C/G)
variants for which both the exposure minor-allele frequency and the outcome effect-allele
frequency fell between 0.42 and 0.58 were excluded because the strand could not be resolved
reliably; variants without an outcome record or with incompatible alleles were also removed.

After harmonisation, 2,278 valid instruments remained for the exploratory analysis across all
115 traits (median 18 per trait; range 2–65). Six traits contained fewer than ten instruments
(indices 003, 011, 063 and 065 with nine each, 113 with eight, and 055 with two). For the strict
analysis, 43 valid instruments covered 32 traits; 25 of these traits had a single instrument,
seven had two or more (of which three had at least three, with a maximum of four).

## Instrument strength

The strength of each instrument was quantified by the per-variant F statistic, calculated as
(βₓ/SEₓ)², where βₓ and SEₓ are the exposure effect and its standard error. Across all 2,278
harmonised instruments the median F was 21.3 and no instrument had F<10 (the smallest observed
value was approximately 17.4); for the Cutaneotrichosporon continuous trait the mean F was 22.6
(minimum 19.8). Because an F above 10 is to a large extent a mechanical consequence of selecting
variants at P<10⁻⁵, where the squared z score is already close to 19, we treat this statistic as
a weak-instrument screen rather than as evidence of high statistical power, and we report the
variance explained and minimum detectable effect alongside it.

## Two-sample MR estimation

Estimates are expressed as the change in the log-odds of RA per one-standard-deviation increase
in genetically predicted CLR abundance for continuous traits, or per unit increase in the log-odds
of detection for binary traits. For traits with a single valid instrument, the causal estimate was
the Wald ratio (βᵧ/βₓ). For traits with multiple instruments, the primary estimate was the
fixed-effects inverse-variance weighted (IVW) method, which is a weighted regression of the
outcome effects on the exposure effects constrained through the origin, with weights equal to the
inverse of the squared outcome standard error. We implemented IVW directly as
β = Σw βₓβᵧ / Σw βₓ² with standard error 1/√(Σw βₓ²). Rather than inflating the fixed-effects
standard error through a multiplicative random-effects term, we treated over-dispersion as a
diagnostic signal in its own right and addressed it through Cochran's Q, MR-Egger, RAPS and
MR-PRESSO, and through leave-one-out and regional exclusion; the causal conclusion is established
from these sensitivity analyses rather than from the nominal fixed-effects precision.

Two complementary estimators were used. The weighted-median estimator returns the median of the
per-variant Wald ratios under inverse-variance weights and can provide a consistent estimate even
when up to half of the weight comes from invalid instruments, provided the valid instruments are
not themselves the minority; its standard error was obtained from 1,000 bootstrap resamples that
resample the variant set with replacement (a single shared index resamples each whole
instrument), with a fixed random seed of 42.
MR-Egger fits an unconstrained weighted regression whose intercept estimates directional
pleiotropy (the mean horizontal-pleiotropic effect) and whose slope estimates the causal effect under an Instrument Strength Independent
of Direct Effects assumption, although it is known to be less precise and can be sensitive when
a single instrument dominates the weights. The intercept standard error was computed on the
residual scale (the weighted residual sum of squares divided by the number of instruments minus
two), matching the default TwoSampleMR implementation. Agreement across estimators was interpreted as
supportive rather than definitive, because estimators that share the same dominant variant are
not independent. All P values used z statistics under a normal approximation (with unrounded
intercept and standard error), rather than a small-sample t distribution, which gives a slightly
larger intercept-test value.

## Pleiotropy and sensitivity analyses

Heterogeneity across instruments was assessed with Cochran's Q statistic, comparing each
variant against the IVW fit. We applied the MR Pleiotropy Residual Sum and Outlier
(MR-PRESSO) method with 1,000 simulated datasets: the global test evaluates excess pleiotropic
residuals, variants with studentised residuals of magnitude at least three were flagged as
outliers, and the distortion test compares the estimate before and after outlier removal. The
leave-one-out analysis re-estimated the IVW effect after dropping each instrument in turn to
identify estimates driven by a single variant. Directional testing of the single influential
variant was examined with the Steiger approach, but its result proved sensitive to the variance
scale and it is therefore treated as an auxiliary consistency check rather than as independent
evidence. Three scale pairings were evaluated transparently: the TwoSampleMR default, pairing the
exposure observed R² with the RA liability-scale R² from the Lee transformation, indicated
reverse causation; a consistent z/p formulation using the full RA sample size (N=178,616)
indicated forward causation; and the same z/p formulation using the RA effective sample size
(approximately 20,751) again indicated reverse causation. Because the conclusion flips according
to a near-boundary choice of scale rather than being stable, no directional inference is drawn
from it. MR-Egger intercepts were inspected for directional pleiotropy, and a robust adjusted profile
score (RAPS) was computed as an iteratively reweighted profile-score regression with Huber
weights (tuning constant 1.345), initialised at the IVW estimate, which down-weights rather than
deletes large-residual variants.

## Signal dissection and variance normalisation

Because the apparent strong signal arose from a single low-frequency variant in the HLA region,
we quantified the variance explained on a stable basis. The point-estimate formula
2 MAF(1−MAF) β² is unreliable for rare variants in sparse, zero-inflated CLR data because an
unstable large β is squared and amplified; we therefore used the association-strength formulation
R² = z²/(z²+N), an N-anchored estimate determined by the observed z score that is stable across
measurement scales, as
the primary measure for comparing the strength of the same variant across the two traits, and it
is not interpreted as a population R². For the exposure, N is the fungal sample size of 7,350;
for the outcome, the denominator is reported separately for the full sample (N=178,616) and the
effective sample size (approximately 20,751), because the two conventions place the variant on
opposite sides of the near-boundary Steiger comparison. Outcome variance on the liability scale
was estimated using
the Lee (2012) framework, correcting the case-enriched sample allele frequency to the population
frequency under an assumed RA population prevalence of K=0.01 (with K=0.005 as a sensitivity
value) and applying V_G/(V_G+π²/3). We additionally re-ran the entire analysis after removing
all instruments in the HLA interval (chromosome 6, 25–35 Mb), and annotated the driving variant
using Ensembl/VEP and its local haplotype structure through the Ensembl LD service in the
Japanese (JPT) reference.

## Positive and negative controls

To confirm that the pipeline could recover a known causal relationship and could expose the
specific bias affecting naive IVW analyses, we used body mass index as an exposure and type 2
diabetes as the outcome, clumping BMI instruments at P<5×10⁻⁸ and harmonising 98 independent
variants (median F 37.5, minimum 29.8). We report the naive IVW estimate, the IVW estimate after
a single round in which variants with studentised residuals above three were removed, and the
weighted-median estimate, expecting the latter two to recover the established positive direction
even when the naive estimate is distorted by diabetes-specific loci. BMI on RA was retained as a
power-calibration exercise using the same minimum-detectable-effect calculation.

As a negative control we used cataract, an age-related, non-immune phenotype with no expected
link to fungal-driven autoimmunity, and evaluated three layers. The first was the MR estimate for
the flagged fungal trait across its full instrument set; the second was the effect of the driving
variant alone, compared with its effect on RA; and the third was the genome-wide distribution of
all 2,109 matched fungal-instrument positions, including the observed versus expected counts of
small P values, the Benjamini–Hochberg results and the genomic-control inflation factor. A
single non-immune control cannot establish strict immunological specificity, but a null result
across all three layers argues against a generic technical association that would appear for any
outcome.

## Statistical power and minimum detectable effect

For each trait we derived the minimum detectable effect (MDE) at 80% power and a two-sided
α=0.05 as 2.80 times the IVW standard error (the corresponding nominal significance threshold is
1.96 standard errors). Across the 115 exploratory traits the median MDE in absolute standardised
effect was 0.083 (interquartile range 0.056–0.112); the median was 0.109 for continuous traits
and 0.062 for binary traits. These values frame the scope of the null findings: effects below
roughly 0.06 cannot be excluded.

## Multiple testing and software

Multiple-testing control used the Benjamini–Hochberg false-discovery rate (FDR) applied within
each threshold family, with FDR<0.05 considered evidence of an association. Analyses were
implemented in Python 3 using NumPy, pandas, SciPy, statsmodels and scikit-learn, with remote
genotype queries through pysam; the figure style used a colour-blind-safe palette. All analysis
scripts and the harmonised datasets are archived with the project for reproducibility.
