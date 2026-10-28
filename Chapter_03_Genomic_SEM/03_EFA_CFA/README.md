# 03 - Exploratory and Confirmatory Factor Analysis

## Purpose

This directory documents the exploratory factor analysis (EFA) and confirmatory factor analysis (CFA) undertaken as part of the Genomic Structural Equation Modelling (Genomic SEM) analyses in Chapter 3. The objective was to examine the genetic covariance structure of epilepsy subtypes and subcortical brain volumes, evaluate alternative latent-factor specifications, and identify interpretable models for subsequent multivariate GWAS.

The analyses used LDSC-derived genetic covariance matrices (`S`) and their sampling covariance matrices (`V`). The provenance and verification of these inputs are documented separately in `02_GenomicSEM_LDSC`.

## Phenotype selection and model development

GGE and FE subtype groups were initially examined separately because joint modelling of all epilepsy and brain-volume phenotypes produced a complex covariance structure. The initial models included seven subcortical brain-volume phenotypes: accumbens, amygdala, brainstem, caudate, pallidum, putamen and thalamus. Following model evaluation, brainstem was excluded from the retained models.

| Model | Epilepsy indicators | Retained brain-volume indicators |
|---|---|---|
| GGE | CAE, GTCSA, JAE, JME | Accumbens, amygdala, caudate, pallidum, putamen, thalamus |
| FE | F-HS, F-Neg, F-Oth | Accumbens, amygdala, caudate, pallidum, putamen, thalamus |

The computational labels `GTCS`, `F_HS`, `F_Neg` and `F_Oth` correspond to GTCSA, F-HS, F-Neg and F-Oth in the thesis.

## Exploratory factor analysis

EFA was conducted in R using `factanal()` to investigate the latent genetic covariance structure and inform candidate CFA specifications.

For GGE, one- to four-factor solutions were evaluated. The one-factor solution was unrotated, while the two- to four-factor solutions used promax rotation. Model evaluation considered factor-loading patterns, cross-loadings, interpretability, variance explained and statistical admissibility.

Where the LDSC genetic covariance matrix was non-positive-definite, `Matrix::nearPD()` was used to obtain a positive-definite covariance matrix for EFA.

EFA findings were used to inform subsequent CFA model development rather than interpreted as independent validation of the factor structures.

The GGE EFA analyses were reproduced using the corrected August LDSC covariance object (`GGE_ldsc_corrected_prevalence_Aug2026.rds`). The reproduced two-, three- and four-factor loading patterns and variance summaries agreed with Appendix Tables A5–A7.

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

### GGE CFA model comparison

One-, two- and three-factor CFA specifications were evaluated using the same corrected LDSC input.

| Model | χ² | df | AIC | CFI | SRMR |
|---|---|---|---|---|---|
| One-factor | 240.18 | 36 | 278.18 | 0.640 | 0.254 |
| Two-factor | 82.47 | 35 | 122.47 | 0.916 | 0.075 |
| Three-factor | 78.72 | 33 | 122.72 | 0.919 | 0.072 |

The one-factor specification modelled epilepsy subtypes and subcortical brain volumes using a single common factor. The two-factor specification separated GGE liability from shared subcortical brain-volume variation. The three-factor specification further separated the brain-volume indicators into striatal/thalamic and limbic factors.

The residual variance of CAE was constrained to zero in all three specifications.

Although the three-factor model showed marginally improved CFI and SRMR, the two-factor model was retained based on its parsimony, interpretability and slightly lower AIC.

All three CFA specifications were reproduced using the corrected August LDSC object, with model-fit statistics agreeing with Appendix Table A8.

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

The retained FE CFA model was reproduced using the archived LDSC covariance object (`FE_subtypes_brain_ldsc_July2026.rds`). The resulting model-fit statistics, standardised factor loadings, standard errors and latent-factor correlation agreed with Appendix Tables A14 and A15.

## Sensitivity analyses

Additional CFA sensitivity analyses were conducted to assess model stability and the influence of individual indicators.

For GGE, these included excluding CAE and evaluating an alternative specification that separated the six retained brain volumes into correlated striatal/thalamic and limbic factors.

For FE, sensitivity analyses separately excluded F-Oth and F-HS. Where only two FE indicators remained, their factor loadings were constrained to equality for model identification.

The full model comparisons and sensitivity findings are reported in the Chapter 3 appendix. The standalone reproduction of the remaining sensitivity analyses is not yet documented in this directory.

## Reproducibility and archived inputs

The original model-development scripts, including `GGE_modelling_July2026.R` and `FE_modelling_July2026.R`, contain exploratory and retained model specifications.

The verified reproducibility scripts in this directory use archived LDSC covariance objects and do not rerun the original summary-statistic preparation or LDSC analyses.

### LDSC inputs

**GGE analyses**

The final GGE EFA and CFA analyses used the corrected August LDSC covariance object:

`GGE_ldsc_corrected_prevalence_Aug2026.rds`

The corrected population prevalences were:

| Phenotype | Population prevalence |
|---|---|
| CAE | 0.00035 |
| GTCSA | 0.00020 |
| JAE | 0.00015 |
| JME | 0.00015 |

The corrected LDSC object contains 11 phenotypes, including brainstem, which is excluded from the retained six-region CFA specifications.

**FE analyses**

The retained FE CFA used:

`FE_subtypes_brain_ldsc_July2026.rds`

This object contains three FE subtypes and seven brain-volume phenotypes, including brainstem.

### Verified scripts

| Script | Purpose | Thesis reference |
|---|---|---|
| `01_fit_GGE_retained_CFA.R` | Reproduce the retained two-factor GGE CFA | Tables A8–A9 |
| `02_fit_FE_retained_CFA.R` | Reproduce the retained two-factor FE CFA | Tables A14–A15 |
| `03_run_GGE_EFA.R` | Reproduce the one- to four-factor GGE EFA solutions | Tables A5–A7 |
| `04_compare_GGE_CFA_models.R` | Compare the one-, two- and three-factor GGE CFA specifications | Table A8 |

### Running the analyses

The scripts require R and the relevant packages, including `GenomicSEM` for CFA and `Matrix` for covariance-matrix smoothing in EFA.

From the `03_EFA_CFA` directory, execute the scripts using the appropriate archived LDSC objects.

**Retained GGE CFA**

```bash
Rscript scripts/01_fit_GGE_retained_CFA.R /path/to/GGE_ldsc_corrected_prevalence_Aug2026.rds
```

**Retained FE CFA**

```bash
Rscript scripts/02_fit_FE_retained_CFA.R /path/to/FE_subtypes_brain_ldsc_July2026.rds
```

**GGE exploratory factor analysis**

```bash
Rscript scripts/03_run_GGE_EFA.R /path/to/GGE_ldsc_corrected_prevalence_Aug2026.rds
```

**GGE CFA model comparison**

```bash
Rscript scripts/04_compare_GGE_CFA_models.R /path/to/GGE_ldsc_corrected_prevalence_Aug2026.rds
```

### Outputs

Each script generates its own output directory relative to the location of the script.

| Analysis | Output directory |
|---|---|
| Retained GGE CFA | `outputs/GGE/` |
| Retained FE CFA | `outputs/FE/` |
| GGE EFA | `outputs/GGE_EFA/` |
| GGE CFA comparison | `outputs/GGE_CFA_comparison/` |

Depending on the analysis, generated outputs include model-fit statistics, parameter estimates, factor loadings, variance summaries, fitted model objects, model syntax, input provenance and R session information.

### Verification

All four scripts were tested independently using their respective archived LDSC objects.

The retained GGE and FE CFA scripts reproduced the reported model-fit statistics and parameter estimates. The GGE EFA script reproduced the factor-loading patterns and variance summaries in Appendix Tables A5–A7. The GGE CFA comparison script reproduced all three model specifications and the fit statistics reported in Appendix Table A8.

Input preparation and LDSC verification are documented separately in `02_GenomicSEM_LDSC`.

Original covariance objects and large intermediate data files are retained in research archives rather than distributed in this repository.
