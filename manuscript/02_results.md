# Results

## Instruments and instrument strength

Across the 115 genus-level fungal traits, variant selection at the exploratory threshold
(P<1×10⁻⁵) followed by true r² clumping in the East Asian reference panel yielded 3,230
instrument–phenotype associations covering every trait (median 25 per trait, minimum 10), which
mapped to 3,009 distinct genomic positions. After alignment with the BBJ RA summary file, 2,278
valid instruments remained (median 18 per trait; range 2–65). Six traits retained fewer than ten
instruments (indices 003, 011, 063 and 065 with nine each, 113 with eight, and 055 with two).
At the conventional threshold (P<5×10⁻⁸), 67 clumped associations across 39 traits fell to 43
valid instruments across 32 traits after harmonisation; 25 of these traits had a single
instrument and seven had two or more, of which three had at least three (maximum four).

Instrument strength was broadly comparable across the full set. The per-variant F statistic had
a median of 21.3, and no instrument fell below 10 (the smallest observed value was 17.4). As
anticipated, these values reflect the selection threshold as much as genuine first-stage strength,
because an F near 19 is already implied by a P value of 10⁻⁵, so we interpret them as a
weak-instrument screen rather than as a power guarantee. On the variance scale, the median
per-instrument R² was 0.29%, and the median total R² across all instruments for a trait was
5.1% (range 0.6%–18.2%), indicating that the genetic instruments collectively captured a small
but non-negligible fraction of the variation in each fungal trait while leaving most of the
variance unexplained.

## Primary two-sample MR estimates

At the exploratory threshold, five of the 115 traits showed a nominal association with RA
(P<0.05), but only one survived Benjamini–Hochberg correction. That signal was the continuous
CLR abundance of Cutaneotrichosporon (trait 109), for which 22 instruments gave an IVW estimate
of β=−0.167 (SE 0.021; 95% CI −0.209 to −0.126; P=2.64×10⁻¹⁵; FDR=3.04×10⁻¹³), indicating that a
one-standard-deviation increase in genetically predicted abundance was associated with lower RA
log-odds. The remaining four nominal associations did not survive correction. Pichia (13
instruments) gave β=0.071 (CI 0.005 to 0.137; P=0.035; FDR=0.73) with no detectable heterogeneity
(P=0.31). The binary Cutaneotrichosporon detection trait (16 instruments) gave β=0.058
(CI 0.007 to 0.109; P=0.026; FDR=0.73). Mucor, based on only two instruments, gave β=−0.191
(CI −0.336 to −0.045; P=0.010; FDR=0.40). Nigrospora (31 instruments) gave β=−0.115
(CI −0.192 to −0.038; P=0.0034; FDR=0.20). Candida, the genus most often discussed in relation
to autoimmune disease, showed no association (β=−0.022, P=0.099), and the wider distribution of
the 115 IVW estimates was centred close to the null rather than showing a systematic directional
shift.

At the strict threshold, four traits were nominally associated, all through single-instrument
Wald ratios: Penicillium (β=0.493, P=0.023), Meyerozyma (β=−0.140, P=0.026), Acremonium
(β=−0.419, P=0.032) and Zygosaccharomyces (β=−0.405, P=0.044). A single-instrument ratio permits
no heterogeneity, pleiotropy or leave-one-out check, so each of these rests entirely on the
validity of one variant and is particularly vulnerable to the same HLA-driven confounding
documented below. After standard Benjamini–Hochberg correction across the 32 traits, their
adjusted values lay between 0.34 and 0.35 (three at 0.339 and Zygosaccharomyces at 0.351), so
none was significant.

## Sensitivity estimators

For the Cutaneotrichosporon continuous trait, the weighted-median point estimate was −0.244 with
a bootstrap standard error of 0.117 (P=0.037), and the MR-Egger slope was −0.205
(P=3.5×10⁻⁹ on the residual scale). The robust adjusted profile score, which down-weights variants with large
residuals through a Huber loss, gave −0.129 (P=9.8×10⁻¹⁰). All three estimators share the
direction of the IVW result, but they do not constitute independent confirmation. The
weighted-median point estimate is exactly the Wald ratio of the single dominant variant, because
that variant carries 70% of the total weight and the estimator collapses onto it; the MR-Egger
slope is likewise dominated by the same outlier and is uninformative once the intercept is
modelled from so few variants; and the robust profile score, while less extreme, still retains
the variant in its loss function. Agreement across methods is therefore expected given the
shared dominant instrument rather than evidence of a robust effect, and the relevant
discriminating evidence is what happens once that instrument is removed.

## Heterogeneity and pleiotropy

The 22 instruments for the Cutaneotrichosporon trait showed substantial heterogeneity
(Cochran's Q=58.19; P=2.4×10⁻⁵), inconsistent with a single homogeneous causal effect. The
MR-Egger intercept was 0.032 (SE 0.013; P=0.011 under the normal-z approximation using the
residual-scale standard error; a t-based value with 20 degrees of freedom would be about 0.019).
This is a nominally significant intercept that raises the possibility of directional pleiotropy,
but it is not corrected for the 115-trait screen and is driven by the same high-leverage variant.
Cochran's Q and MR-PRESSO, which detect horizontal pleiotropy and heterogeneity rather than its
direction, independently indicated pleiotropy, with MR-PRESSO global P=0.001. Because the dominant
variant carries most of the weight, the intercept is best interpreted alongside the outlier
analysis rather than as an independent test. Under the pre-specified studentised
residual threshold of three, three variants were flagged, and the dominant variant had the
strongest residual by a wide margin (−10.1), with the next strongest at 4.3. Removing the three
strongly outlying variants left 19 instruments and returned the IVW estimate to −0.015
(SE 0.056; P=0.79), and the distortion test indicated a significant difference between the raw
and corrected estimates (P<0.05 across the standard-error specifications examined). The
convergence of the heterogeneity, residual and corrected-estimate results on the same small set
of variants localises the pleiotropy rather than distributing it across the instrument set.

## Directionality (Steiger)

The Steiger direction test for the dominant variant was not stable. The TwoSampleMR default,
pairing the exposure observed R² (0.0027) with the RA liability-scale R² from the Lee
transformation (0.0031), indicated reverse causation. A consistent z/p formulation using the
full RA sample size (N=178,616) gave an outcome R² of 0.00052 and indicated forward causation.
The same formulation using the effective RA sample size (approximately 20,751) gave an outcome
R² of 0.0045 and again indicated reverse causation. Because the call flips with a near-boundary
choice of scale and sample-size denominator, and all these R² values are small and straddle the
boundary at which the direction reverses, we do not draw a directional conclusion from this
test and treat it only as a sensitivity indicator rather than as an independent evidence path.

## Signal dissection and variance normalisation

The leave-one-out analysis localised the entire Cutaneotrichosporon signal to one variant,
rs111394117 on chromosome 6 (position 32,187,213), located in the ninth intron of NOTCH4 within
the HLA complex. On the exposure side this is a rare variant (MAF 0.021) with an unstable point
estimate (β=1.70, SE 0.38; z=4.47; P≈8.8×10⁻⁶), well below the exposure study-wide threshold of
2.2×10⁻¹⁰. On the outcome side it is a strongly established RA locus (β=−0.414, SE 0.043;
effect-allele frequency 0.063; z=−9.64; P=5.27×10⁻²²). It is therefore a genome-wide-significant
RA variant that is only a sub-threshold fungal quantitative-trait locus. On the
association-strength (method B) estimate, it explained 0.27% of the fungal trait and 0.31% of
RA liability, so its association with RA is at least as strong as its association with the
proposed exposure.

The choice of variance formula was validated across all 2,278 instruments rather than only at
this variant. The point-estimate formula and the association-strength formula were applied to
every instrument: the point-estimate formula produced no direction flips, whereas the
association-strength formula produced exactly one, and that single flip was the dominant HLA
variant, for which the point estimate assigned an implausible 11.9% of variance. On ordinary
variants the point-estimate formula systematically overstated variance, with a median inflation
of 1.63-fold, confirming that the discrepancy is not confined to the outlier but reflects a
general bias in the point-estimate approach for these traits. The association-strength formula
is therefore used throughout, and the point-estimate value is retained only as a documented
sensitivity comparator.

When this single variant was removed, the remaining 21 instruments gave β=0.012 (SE 0.039;
P=0.763), and the signal disappeared. A global analysis excluding the entire HLA interval
(chromosome 6, 25–35 Mb) produced no FDR-significant association in either set: zero of 32
strict traits and zero of 115 exploratory traits, with the smallest exploratory FDR at 0.763.
Haplotype annotation in the Japanese reference showed the variant sits in a large HLA block (66
proxy relationships at r²≥0.8, four at r²=1.0), consistent with NOTCH4 acting as a marker rather
than the causal gene, most plausibly tagging class II DRB1/DQB1 variation.

The two measurement scales for the same genus gave opposite directions. The binary detection
trait (051), based on 16 instruments, yielded β=0.058 (SE 0.026; P=0.026; heterogeneity
P=0.602), and its weighted median was concordant at β=0.065 (bootstrap SE 0.028; P=0.020). The
dominant variant was absent from its instrument set altogether, because that variant reaches
candidate status only for the continuous CLR phenotype. A genuine causal effect should propagate
across defensible operationalisations of the exposure; instead the two scales point in opposite
directions and the only strong, corrected signal is confined to the continuous scale and to one
HLA variant that exists only on that scale. This is consistent with an artefactual signal arising
from a single pleiotropic HLA variant rather than a biological effect of the fungus itself.

## Positive and negative controls

The positive control behaved as expected and exposed the specific failure mode of uncorrected
IVW. For BMI on type 2 diabetes, 100 independent BMI loci were reduced to 98 harmonised
instruments (median F 37.5, minimum 29.8). The naive IVW estimate was β=−0.288 (SE 0.040;
P=8.56×10⁻¹³), pointing in the wrong direction. The largest standardised residual reached 27.7,
and the responsible variants were the well-established diabetes-specific loci TCF7L2, KCNQ1 and
CDKN2A/B, which carry direct effects on glycaemic traits and therefore violate the
exclusion-restriction assumption when treated as instruments for adiposity. A single round of
residual-based removal flagged 42 such variants and retained 56, after which the IVW estimate
reversed to β=0.184 (SE 0.064; P=0.0042). The weighted median independently gave β=0.589
(SE 0.157; P=1.81×10⁻⁴), and the MR-Egger intercept was near zero (−0.0039). Both corrected
estimators thus recover the established positive effect of BMI on diabetes, demonstrating that
the pipeline can recover a known effect in the presence of pleiotropic contamination and that
the residual-based correction behaves as expected. BMI on RA (50 instruments) gave β=−0.052
(SE 0.110; P=0.64) with an 80%-power minimum detectable effect of 0.31, indicating that the null
BMI–RA estimate is a power limitation rather than evidence of absence, and providing a benchmark
for the fungal-trait calibration.

The negative control, cataract, showed no signal across all three layers, and the layered design
was chosen so that no single statistic carried the inference. The first layer tested the same
instrument set that produced the RA signal: across the 22 Cutaneotrichosporon instruments the
IVW estimate was −0.010 (SE 0.009; P=0.291; 95% CI −0.027 to 0.008) with no heterogeneity
(P=0.786), so a relationship that was striking for RA was absent for the eye phenotype. The
second layer isolated the variant responsible for the RA signal, which had β=−0.011 (SE 0.018;
P=0.558) for cataract, an effect around 39-fold smaller than its RA effect and statistically
indistinguishable from zero. The third layer asked whether the fungal instrument positions as a
whole were enriched for small P values. Across all 2,109 matched positions the counts matched
their null expectation (109 observed versus 105 expected at P<0.05 and 17 versus 21 at P<0.01),
no association survived Benjamini–Hochberg, the minimum P was 5.18×10⁻⁴, and the
genomic-control inflation factor was 0.944, indicating slight deflation rather than any
systematic signal. Together these results argue against a generic technical association that
would appear for any outcome and against population stratification that would inflate all
phenotypes, although a single non-immune control cannot by itself establish strict immunological
specificity and we phrase the conclusion accordingly.

## Power and minimum detectable effect

At 80% power, the median minimum detectable effect across the 115 exploratory traits was 0.083
(interquartile range 0.056–0.112); it was 0.109 for continuous traits and 0.062 for binary
traits. The lower threshold for binary traits reflects their generally more precise instruments,
whereas the CLR continuous traits, with their sparse zero-inflated distributions, carry wider
standard errors. Seventy-six traits (66%) were powered to detect an effect below 0.10 and 39
below 0.06. The null findings therefore exclude moderate-to-large effects across the panel,
including effect sizes of the magnitude that would be considered biologically meaningful for a
log-odds change in disease risk. Standardised effects smaller than roughly 0.06, however, cannot
be ruled out and would require substantially larger exposure samples than the current N of 7,350;
this residual uncertainty is most relevant for the continuous traits, whose minimum detectable
effect is largest (0.109) and whose power is correspondingly lowest.

## Multiple-testing summary

After Benjamini–Hochberg correction within each threshold family, one association (the continuous
Cutaneotrichosporon trait) was significant at the exploratory threshold and none was significant
at the strict threshold. That single corrected association is fully attributable to one HLA
variant that is itself an established RA locus, is absent from the binary measurement scale, and
disappears when the HLA region is removed. Read together with the null cataract control, the
residual-based recovery of the positive BMI–diabetes control, and the power calibration that
excludes moderate-to-large effects, the study provides no robust evidence of a causal effect of
the gut mycobiome on RA, while documenting a worked example of how HLA pleiotropy can generate
a deceptively strong Mendelian-randomization signal.
