# 02 — Multivariable LDSC for Genomic SEM

## Purpose

This directory documents the multivariable linkage disequilibrium score regression (LDSC) covariance inputs used for the Chapter 3 Genomic SEM factor models. It is distinct from `01_standalone_LDSC`, which documents the individual heritability estimates, pairwise genetic correlations and Figure 3.1.

Multivariable `GenomicSEM::ldsc()` provides the genetic covariance matrix (`S`), the sampling covariance matrix of its unique elements (`V`), and associated intercept and sample-size information (`I`, `N`, `m`). These quantities are inputs to the exploratory and confirmatory factor analyses and subsequent multivariate GWAS.

## Phenotype sets and archived inputs

Initial July 2026 LDSC analyses included seven subcortical brain-volume phenotypes: accumbens, amygdala, brainstem, caudate, pallidum, putamen and thalamus. The retained models used the six regions other than brainstem.

| Group | Original LDSC input object | Original dimensions (`S`; `V`) | Retained model input | Retained dimensions (`S`; `V`) |
| --- | --- | --- | --- | --- |
| GGE: CAE, GTCS, JAE, JME | `GGE_ldsc_July2026.rds` | 11 × 11; 66 × 66 | `ldsc_GGE_brain_6region.rds` | 10 × 10; 55 × 55 |
| FE: F_HS, F_Neg, F_Oth | `FE_subtypes_brain_ldsc_July2026.rds` | 10 × 10; 55 × 55 | `ldsc_FE_brain_6region.rds` | 9 × 9; 45 × 45 |

Here, `GTCS` is the computational phenotype label for GTCSA in the thesis. The computational labels `F_HS`, `F_Neg` and `F_Oth` correspond to F-HS, F-Neg and F-Oth in the thesis. Phenotype names are retained as stored in the original objects to preserve matrix indexing and model syntax.

The retained six-region objects have been checked against their July source objects: the corresponding numerical entries in `S`, `I`, `V` and `N` agree exactly after retaining the appropriate phenotypes and matching lower-triangular covariance-element indices. The retained GGE covariance object was also used to reproduce the standardized CFA estimates in Appendix Table A9; the FE object reproduced the CFA estimates in Appendix Table A15 to the reported precision. These CFA results are addressed separately in `03_EFA_CFA`.

## LDSC inputs and parameters

The original July R scripts use European LD scores and regression weights from `eur_w_ld_chr/`. GWAS summary statistics were harmonised using `GenomicSEM::munge()`, aligned to a HapMap3 reference SNP list and filtered at MAF ≥ 0.01. SNP-specific sample sizes from the source summary statistics were retained where present. A sample case fraction of 0.50 was specified for epilepsy phenotypes, reflecting the effective-sample-size convention; prevalence arguments were set to `NA` for continuous brain-volume phenotypes.

The July scripts specified the following population prevalences:

| Epilepsy trait | July LDSC population prevalence |
| --- | ---: |
| CAE | 0.00015 |
| GTCS (GTCSA) | 0.00020 |
| JAE | 0.00015 |
| JME | 0.00035 |
| F_HS | 0.00075 |
| F_Neg | 0.00090 |
| F_Oth | 0.00135 |

The original munging log records that an INFO filter could not be applied to input files lacking an INFO field; the requested MAF threshold and reference-allele checks are documented in that log. Source GWAS files and reference LD files are not bundled in this directory.

## July and corrected August covariance versions

The corrected August LDSC analyses used **CAE = 0.00035** and **JME = 0.00015**, rather than the July settings shown above. The August GGE covariance estimates are not interchangeable with the July-derived covariance object used to reproduce the retained GGE CFA in Appendix Table A9. The August-corrected objects are separate analytical versions; replacing the retained model's input would constitute a new sensitivity or reanalysis and would require checking the downstream estimates and model fit.

The purpose of the scripts here is to verify the provenance of the **reported** models, not to imply that the corrected August specification was used to generate them. Both versions should be retained in archived research storage with their original filenames and execution logs.

## Reproducibility check

The script `scripts/01_verify_six_region_objects.R` checks the full July GGE and FE LDSC objects against the retained six-region inputs. It tests all entries of `S`, `I`, `V` and `N` after matching phenotype names and the ordering of the lower-triangular covariance elements in `V` and `N`.

From the repository root, with R installed:

```bash
Rscript Chapter_03_Genomic_SEM/02_GenomicSEM_LDSC/scripts/01_verify_six_region_objects.R \
  "/path/to/GGE_ldsc_July2026.rds" \
  "/path/to/ldsc_GGE_brain_6region.rds" \
  "/path/to/FE_subtypes_brain_ldsc_July2026.rds" \
  "/path/to/ldsc_FE_brain_6region.rds"
```

The script uses **base R only** and writes a verification table to `outputs/six_region_object_verification.csv` within this directory. It does **not** rerun LDSC, perform matrix smoothing, overwrite the saved `.rds` files or fit CFA models. Numerical agreement documents consistency between archived objects; it does not independently establish every upstream GWAS processing step.

## Archived-data handling

Keep the original `.rds` files, harmonised GWAS statistics, LD reference files and execution logs in appropriately controlled research storage. Do not commit restricted summary statistics or large intermediate data to the public repository. The scripts can be run by providing paths to the archived files. Generated verification outputs are derived reports and need not be committed.
