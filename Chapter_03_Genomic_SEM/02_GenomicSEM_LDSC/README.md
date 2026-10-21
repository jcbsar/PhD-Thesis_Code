# 02 - Multivariable LDSC for Genomic SEM

## Purpose

This directory documents the multivariable linkage disequilibrium score regression (LDSC) covariance inputs used for the Genomic Structural Equation Modelling (Genomic SEM) analyses in Chapter 3. It is distinct from `01_standalone_LDSC`, which documents the individual SNP-based heritability estimates, pairwise genetic correlations and Figure 3.1.

Multivariable LDSC was performed using `GenomicSEM::ldsc()` to estimate the genetic covariance matrix (`S`), the sampling covariance matrix of its unique elements (`V`), and associated intercept and sample-size information (`I`, `N`, `m`). These outputs provide the covariance structure for subsequent exploratory factor analysis (EFA), confirmatory factor analysis (CFA) and multivariate GWAS.

## Phenotype sets and archived inputs

The initial LDSC analyses included seven subcortical brain-volume phenotypes: accumbens, amygdala, brainstem, caudate, pallidum, putamen and thalamus. Following model development, brainstem was excluded from the retained GGE and FE models, leaving six brain-volume phenotypes.

| Group | Original LDSC object | Original dimensions (`S`; `V`) | Retained model input | Retained dimensions (`S`; `V`) |
|---|---|---|---|---|
| GGE: CAE, GTCSA, JAE, JME | `GGE_ldsc_July2026.rds` | 11 × 11; 66 × 66 | `ldsc_GGE_brain_6region.rds` | 10 × 10; 55 × 55 |
| FE: F_HS, F_Neg, F_Oth | `FE_subtypes_brain_ldsc_July2026.rds` | 10 × 10; 55 × 55 | `ldsc_FE_brain_6region.rds` | 9 × 9; 45 × 45 |

`GTCS` is the computational label corresponding to GTCSA in the thesis. Similarly, `F_HS`, `F_Neg` and `F_Oth` correspond to F-HS, F-Neg and F-Oth. The original computational names are preserved to maintain consistency with the archived covariance matrices and model syntax.

Numerical comparisons confirmed that the retained six-region objects match the corresponding elements of the July LDSC objects for `S`, `V`, `I` and `N`, after accounting for phenotype order and the lower-triangular indexing of covariance elements.

The retained GGE and FE objects were subsequently used to reproduce the standardised CFA parameter estimates reported in Appendix Tables A9 and A15, respectively, to the reported precision. The CFA specifications and results are addressed separately in `03_EFA_CFA`.

## LDSC inputs and parameters

The original analyses were performed using European LD scores and regression weights from `eur_w_ld_chr/`.

GWAS summary statistics were harmonised using `GenomicSEM::munge()` with the HapMap3 reference SNP list and a minor allele frequency threshold of MAF ≥ 0.01. SNP-specific sample sizes were retained from the source GWAS summary statistics.

For the epilepsy phenotypes, a sample case fraction of 0.50 was specified, reflecting the effective-sample-size convention used in the source GWAS. Population prevalences were specified for liability-scale transformation. Both sample and population prevalence parameters were set to `NA` for continuous brain-volume phenotypes.

### Population prevalences used in the July LDSC analyses

| Epilepsy phenotype | Population prevalence |
|---|---:|
| CAE | 0.00015 |
| GTCS (GTCSA) | 0.00020 |
| JAE | 0.00015 |
| JME | 0.00035 |
| F_HS | 0.00075 |
| F_Neg | 0.00090 |
| F_Oth | 0.00135 |

An archived munging log, `GGE_subtypes_brainvolumes_July2026_munge.log`, records execution on 12 August 2026. The log documents SNP retention following harmonisation and notes that INFO filtering could not be applied to source files without an INFO field. Reference-allele matching and MAF filtering were performed as recorded in the log.

Original GWAS summary statistics and LD reference files are not distributed in this repository.

## July and corrected August LDSC versions

Subsequent corrected August LDSC analyses used updated population prevalence specifications for two GGE subtypes:

| Phenotype | July specification | Corrected August specification |
|---|---:|---:|
| CAE | 0.00015 | 0.00035 |
| JME | 0.00035 | 0.00015 |

The corrected August LDSC analyses produced different covariance estimates for affected GGE phenotypes.

The retained GGE CFA estimates reported in Appendix Table A9 were reproduced using `ldsc_GGE_brain_6region.rds`, which matches the July-derived covariance estimates rather than the corrected August object.

The corrected August covariance matrices are therefore treated as a distinct analytical version and are not substituted retrospectively into the retained models. Evaluating the consequences of the updated prevalence specification for factor-model parameters or multivariate GWAS results would require separate downstream analyses.

Both LDSC versions should be preserved in research archives with their original filenames and associated execution records.

## Reproducibility and verification

The script `scripts/01_verify_six_region_objects.R` compares the original July LDSC covariance objects with the six-region objects retained for modelling.

It verifies numerical agreement for:

- The genetic covariance matrix (`S`).
- The LDSC intercept matrix (`I`).
- The sampling covariance matrix (`V`).
- The sample-size matrix (`N`).

For `V` and `N`, the comparison accounts for the ordering of unique covariance elements in the lower triangle of `S`.

### Running the verification script

The script requires R and access to the four archived `.rds` files.

Run it from the repository root using the following positional arguments, in the order shown:

```bash
Rscript Chapter_03_Genomic_SEM/02_GenomicSEM_LDSC/scripts/01_verify_six_region_objects.R "/path/to/GGE_ldsc_July2026.rds" "/path/to/ldsc_GGE_brain_6region.rds" "/path/to/FE_subtypes_brain_ldsc_July2026.rds" "/path/to/ldsc_FE_brain_6region.rds"
```

Replace the example paths with the locations of the archived files. The command can also be run from a Windows terminal with Rscript available.

The script uses base R and writes its verification report to:

`outputs/six_region_object_verification.csv`

The verification script does not rerun LDSC, perform matrix smoothing, modify the original `.rds` files or fit CFA models.

Agreement between archived covariance objects establishes consistency of the retained modelling inputs with their identified source objects. It does not independently verify all upstream GWAS processing or alternative LDSC specifications.

## Archived-data handling

Original covariance objects, harmonised GWAS summary statistics, LD reference files and execution logs should remain in appropriately controlled research storage.

Large intermediate files and restricted GWAS data are not included in this repository. The verification script accepts local paths to these archived inputs.

Generated verification reports can be retained locally and do not need to be committed to the repository.
