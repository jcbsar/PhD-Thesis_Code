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
| Latent-factor correlation | −0.166 |

Both CFI and SRMR met the prespecified fit criteria.

The retained GGE CFA model was reproduced using the corrected August LDSC covariance object (`GGE_ldsc_corrected_prevalence_Aug2026.rds`). The resulting model-fit statistics, standardised factor loadings, standard errors and latent-factor correlation agreed with Appendix Tables A8 and A9.

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

Additional CFA sensitivity analyses were conducted to assess model stability and the influence of
