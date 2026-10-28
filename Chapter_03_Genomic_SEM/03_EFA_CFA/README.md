# 03 - Exploratory and Confirmatory Factor Analysis

## Purpose

This directory documents the exploratory factor analysis (EFA) and confirmatory factor analysis (CFA) undertaken as part of the Genomic Structural Equation Modelling (Genomic SEM) analyses in Chapter 3. The objective was to examine the genetic covariance structure of epilepsy subtypes and subcortical brain volumes, evaluate alternative latent-factor specifications, and identify interpretable models for subsequent multivariate GWAS. The analyses used LDSC-derived genetic covariance matrices (`S`) and their sampling covariance matrices (`V`). The provenance and verification of these inputs are documented separately in `02_GenomicSEM_LDSC`.

## Phenotype selection and model development

GGE and FE subtype groups were initially examined separately because joint modelling of all epilepsy and brain-volume phenotypes produced a complex covariance structure. The initial models included seven subcortical brain-volume phenotypes: accumbens, amygdala, brainstem, caudate, pallidum, putamen and thalamus. Following model evaluation, brainstem was excluded from the retained models.

| Model | Epilepsy indicators | Retained brain-volume indicators |
|---|---|---|
| GGE | CAE, GTCSA, JAE, JME | Accumbens, amygdala, caudate, pallidum, putamen, thalamus |
| FE | F-HS, F-Neg, F-Oth | Accumbens, amygdala, caudate, pallidum, putamen, thalamus |

The computational labels `GTCS`, `F_HS`, `F_Neg` and `F_Oth` correspond to GTCSA, F-HS, F-Neg and F-Oth in the thesis.

## Exploratory factor analysis

EFA was conducted in R using `factanal()` to investigate the latent genetic covariance structure and inform candidate CFA specifications. Two- to four-factor solutions were evaluated using promax rotation, with consideration given to factor-loading patterns, cross-loadings, interpretability, variance explained and statistical admissibility. Where the LDSC genetic covariance matrix was non-positive-definite, `Matrix::nearPD()` was used to obtain a positive-definite covariance matrix for EFA. EFA findings were used to inform subsequent CFA model development rather than interpreted as independent validation of the factor structures.

## Confirmatory factor analysis

Candidate CFA models were fitted using `GenomicSEM::usermodel()` with diagonally weighted least-squares (DWLS) estimation and latent variances fixed to one. Models were evaluated using the chi-square statistic, comparative fit index (CFI), standardised root mean square residual (SRMR) and Akaike information criterion (AIC).

CFI ≥ 0.90 and SRMR ≤ 0.08 were considered acceptable fit criteria. Model selection also considered interpretability, statistical admissibility and parsimony. AIC comparisons were restricted to models fitted to the same phenotype sets. Non-positive-definite covariance matrices were handled using the internal smoothing procedures provided by Genomic SEM.

### Retained GGE model

The retained GGE specification comprised two correlated latent factors:

- **GGE factor:** CAE, GTCSA, JAE and JME.
- **Brain-volume factor:** Accumbens, amygdala, caudate, pallidum, putamen and thalamus.

The residual variance of CAE was constrained to zero to address a boundary estimate. The retained model, evaluated using the corrected LDSC population-prevalence specification, produced the following fit statistics:

| Statistic | Estimate |
|---|---|
| χ² | 82.47 |
| Degrees of freedom | 35 |
| CFI | 0.916 |
| SRMR | 0.075 |
| AIC | 122.47 |
| Latent-factor correlation | Approximately −0.17 |

Both CFI and SRMR met the prespecified fit criteria.

The corrected GGE CFA was verified against the prevalence-corrected August LDSC covariance object, reproducing the fit statistics reported in Appendix Table A8. Standardised factor loadings and their standard errors agreed with Appendix Table A9 to the reported precision.

An earlier reproduction using the original July six-region LDSC object produced χ² = 82.54 and AIC = 122.54. Those values reflect the earlier covariance specification and should not be substituted for the corrected results reported above.

### Retained FE model

The retained FE specification comprised two correlated latent factors:

- **FE factor:** F-HS, F-Neg and F-Oth.
- **Brain-volume factor:** Accumbens, amygdala, caudate, pallidum, putamen and thalamus.

The retained model produced the following fit statistics:

| Statistic | Estimate |
|---|---|
| χ² | 70.30 |
| Degrees of freedom | 26 |
| CFI | 0.863 |
| SRMR | 0.077 |
| AIC | 108.30 |
| Latent-factor correlation | −0.095 |

The FE model met the SRMR criterion but not the prespecified CFI threshold, indicating mixed evidence of model fit.

The retained FE CFA model was reproduced using `ldsc_FE_brain_6region.rds`. Its standardised factor loadings and latent-factor correlation agreed with Appendix Table A15 to the reported precision.

## Sensitivity analyses

Additional CFA sensitivity analyses were conducted to assess model stability and the influence of individual indicators.

For GGE, these included excluding CAE and evaluating an alternative specification that separated the six retained brain volumes into correlated striatal/thalamic and limbic factors.

For FE, sensitivity analyses separately excluded F-Oth and F-HS. Where only two FE indicators remained, their factor loadings were constrained to equality for model identification.

The full model comparisons and sensitivity findings are reported in the Chapter 3 appendix.

## Reproducibility and archived inputs

The covariance objects used during model development include both original and subsequently corrected LDSC specifications. These must be distinguished when reproducing the reported analyses.

**Original July covariance objects:**

- `ldsc_GGE_brain_6region.rds`
- `ldsc_FE_brain_6region.rds`

The original GGE object reflects the earlier population-prevalence specification in which CAE and JME prevalences were interchanged. It reproduces the earlier GGE CFA fit statistics but is not the corrected reference for the final GGE model.

**Corrected GGE LDSC specification:**

The corrected population prevalences were:

- CAE: 0.00035
- GTCSA: 0.00020
- JAE: 0.00015
- JME: 0.00015

The August corrected LDSC output was archived as `GGE_ldsc_July2026.rds`, despite the filename referring to July. Its provenance was established by verifying the underlying script, covariance estimates and reproduced CFA results rather than relying on the filename alone.

This corrected object reproduced the final GGE CFA fit statistics reported in Appendix Table A8. It should be distinguished from the separate corrected LDSC object subsequently prepared for the multivariate GWAS.

The original model-development scripts include `GGE_modelling_July2026.R` and `FE_modelling_July2026.R`. These contain exploratory specifications, including models that were not retained.

The initial standalone retained GGE CFA reproduction was conducted using the original July covariance specification. The GitHub reproduction script must therefore be checked and updated to use the corrected LDSC inputs before it can be treated as reproducing the final GGE CFA results.

Additional reproducibility scripts for the retained models and sensitivity analyses should be included following verification against the reported thesis results.

Original covariance objects and large intermediate data files are retained in research archives rather than distributed in this repository.
