# 02 - Multivariable LDSC for Genomic SEM

## Purpose

This directory documents the multivariable linkage disequilibrium score regression (LDSC) covariance inputs used for the Genomic Structural Equation Modelling (Genomic SEM) analyses in Chapter 3. It is distinct from `01_standalone_LDSC`, which documents individual SNP-based heritability estimates, pairwise genetic correlations and Figure 3.1.

Multivariable LDSC was performed using `GenomicSEM::ldsc()` to estimate the genetic covariance matrix (`S`), the sampling covariance matrix of its unique elements (`V`), and associated intercept and sample-size information (`I`, `N`, `m`). These outputs provided the genetic covariance structure for subsequent exploratory factor analysis (EFA), confirmatory factor analysis (CFA) and multivariate GWAS.

This directory documents the original and corrected GGE covariance specifications, the retained FE covariance input, and the independently estimated nine-phenotype LDSC analysis used for the GGE sensitivity analysis excluding childhood absence epilepsy (CAE).

## Phenotype sets and archived inputs

The initial LDSC analyses included seven subcortical brain-volume phenotypes: accumbens, amygdala, brainstem, caudate, pallidum, putamen and thalamus. Following model development, brainstem was excluded from the retained GGE and FE models, leaving six brain-volume phenotypes.

| Analysis | Archived LDSC object | `S` dimensions | `V` dimensions |
|---|---|---|---|
| Original GGE, seven brain regions | `GGE_ldsc_July2026.rds` | 11 × 11 | 66 × 66 |
| Original FE, seven brain regions | `FE_subtypes_brain_ldsc_July2026.rds` | 10 × 10 | 55 × 55 |
| GGE, six brain regions (July-derived) | `ldsc_GGE_brain_6region.rds` | 10 × 10 | 55 × 55 |
| FE, six brain regions | `ldsc_FE_brain_6region.rds` | 9 × 9 | 45 × 45 |
| Corrected GGE, six brain regions | `ldsc_GGE_brain_6region_corrected.rds` | 10 × 10 | 55 × 55 |
| GGE excluding CAE, six brain regions | `ldsc_GGE_no_CAE_brain_6region.rds` | 9 × 9 | 45 × 45 |

The GGE phenotype set comprised CAE, GTCSA, JAE and JME, while the FE phenotype set comprised F-HS, F-Neg and F-Oth. The computational label `GTCS` corresponds to GTCSA in the thesis. Similarly, `F_HS`, `F_Neg` and `F_Oth` correspond to F-HS, F-Neg and F-Oth. These computational labels are preserved for consistency with the archived covariance matrices and model syntax.

Numerical comparisons confirmed that the July-derived six-region GGE and FE objects match the corresponding elements of their original seven-region LDSC objects for `S`, `V`, `I` and `N`, accounting for phenotype order and the lower-triangular indexing of covariance elements.

The six-region July-derived GGE and FE objects were also used to reproduce the standardised CFA parameter estimates reported in Appendix Tables A9 and A15, respectively, to the reported precision. The CFA specifications and results are documented separately in `03_EFA_CFA`.

## LDSC inputs and parameters

The analyses used European LD scores and regression weights from `eur_w_ld_chr/`.

GWAS summary statistics were prepared using `GenomicSEM::munge()` with the HapMap3 reference SNP list and a minor allele frequency threshold of MAF ≥ 0.01. SNP-specific sample sizes were retained from the source GWAS summary statistics.

For the epilepsy phenotypes, a sample case fraction of 0.50 was specified, consistent with the effective-sample-size convention used in the source GWAS. Population prevalences were supplied for liability-scale transformation. Both sample and population prevalence parameters were set to `NA` for the continuous brain-volume phenotypes.

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

Original GWAS summary statistics, munged datasets and LD reference files are not distributed in this repository.

## Corrected GGE LDSC specification and downstream analyses

A subsequent correction updated the population prevalence specifications for two GGE subtypes:

| Phenotype | Original July specification | Corrected specification |
|---|---:|---:|
| CAE | 0.00015 | 0.00035 |
| JME | 0.00035 | 0.00015 |

The revised prevalence assumptions produced a distinct GGE LDSC covariance object, archived as `ldsc_GGE_brain_6region_corrected.rds`.

The original and corrected covariance objects were used in different downstream analyses and should not be treated as interchangeable.

- **Retained CFA:** The GGE CFA estimates reproduced in Appendix Table A9 used `ldsc_GGE_brain_6region.rds`, matching the July-derived covariance estimates.
- **Definitive primary GGE multivariate GWAS:** The corrected GGE covariance object, `ldsc_GGE_brain_6region_corrected.rds`, was used for the final genome-wide analysis.
- **FE analyses:** The FE models used `ldsc_FE_brain_6region.rds`; the GGE prevalence correction did not alter this object.
- **GGE excluding CAE:** The sensitivity analysis used a separately estimated nine-phenotype LDSC object, described below.

These distinctions are maintained to preserve the provenance of each reported analysis. Reproducing the original CFA estimates using the July-derived covariance object should not be interpreted as validation of CFA estimates under the corrected GGE prevalence specification.

Both GGE LDSC versions are retained in the research archives under their original filenames.

## GGE sensitivity analysis excluding CAE

To assess the influence of CAE on the multivariate GGE findings, a sensitivity analysis excluded CAE from the epilepsy factor.

Rather than extracting a reduced covariance matrix from the ten-phenotype GGE object, multivariable LDSC was independently re-estimated using nine phenotypes:

- **GGE indicators:** GTCS, JAE and JME.
- **Subcortical brain-volume indicators:** accumbens, amygdala, caudate, pallidum, putamen and thalamus.

The analysis used the original munged GWAS summary-statistic files and European LD scores and regression weights from `eur_w_ld_chr/`.

### Prevalence assumptions

| Epilepsy phenotype | Sample prevalence | Population prevalence |
|---|---:|---:|
| GTCS | 0.50 | 0.00020 |
| JAE | 0.50 | 0.00015 |
| JME | 0.50 | 0.00015 |

For all six continuous brain-volume phenotypes, both prevalence arguments were `NA`.

The independent LDSC estimation generated `ldsc_GGE_no_CAE_brain_6region.rds`, containing a 9 × 9 genetic covariance matrix (`S`) and a 45 × 45 sampling covariance matrix (`V`), together with the associated LDSC outputs.

The preparation procedure is preserved in:

`scripts/02_prepare_GGE_no_CAE_LDSC.R`

The script retains the original phenotype order, munged input filenames, prevalence assumptions and European LD-score reference configuration.

The CAE-exclusion multivariate GWAS subsequently uses this nine-phenotype covariance object. At the SNP level, `beta.CAE` and `se.CAE` are removed from the harmonised GGE input data without re-harmonising the remaining variants. This operation is implemented in `04_Multivariate_GWAS/02_run_userGWAS.R`; a separate CAE-excluded SNP input file is not required for the consolidated workflow.

## Reproducibility and verification

Two scripts are provided in the `scripts/` directory.

### 1. Verification of retained six-region covariance objects

`scripts/01_verify_six_region_objects.R`

This script compares the original July LDSC covariance objects with the six-region objects retained for modelling.

It verifies numerical agreement for:

- The genetic covariance matrix (`S`).
- The LDSC intercept matrix (`I`).
- The sampling covariance matrix (`V`).
- The sample-size matrix (`N`).

For `V` and `N`, the comparison accounts for the ordering of unique covariance elements in the lower triangle of `S`.

Run the script from the repository root:

```bash
Rscript Chapter_03_Genomic_SEM/02_GenomicSEM_LDSC/scripts/01_verify_six_region_objects.R \
  "/path/to/GGE_ldsc_July2026.rds" \
  "/path/to/ldsc_GGE_brain_6region.rds" \
  "/path/to/FE_subtypes_brain_ldsc_July2026.rds" \
  "/path/to/ldsc_FE_brain_6region.rds"
```

Replace the example paths with the corresponding archived input locations.

The script uses base R and writes its verification report to:

`outputs/six_region_object_verification.csv`

It does not rerun LDSC, perform matrix smoothing, modify the archived covariance objects or fit CFA models.

### 2. Reproduction of the CAE-exclusion LDSC preparation

`scripts/02_prepare_GGE_no_CAE_LDSC.R`

This script documents the independently estimated nine-phenotype LDSC analysis underlying the GGE sensitivity model.

It requires R, Genomic SEM, the nine munged GWAS datasets and the European LD-score reference files.

To verify the presence of the required inputs without rerunning LDSC:

```bash
Rscript Chapter_03_Genomic_SEM/02_GenomicSEM_LDSC/scripts/02_prepare_GGE_no_CAE_LDSC.R \
  "/path/to/gsem" \
  "/path/to/enigma_ilae_sumstats/gsem" \
  --check-only
```

To execute the nine-phenotype LDSC analysis:

```bash
Rscript Chapter_03_Genomic_SEM/02_GenomicSEM_LDSC/scripts/02_prepare_GGE_no_CAE_LDSC.R \
  "/path/to/gsem" \
  "/path/to/enigma_ilae_sumstats/gsem"
```

The first path identifies the analysis root containing `input/`. The second identifies the source directory containing `munged/` and `ldsc_refdata/eur_w_ld_chr/`.

The script checks the required inputs, estimates the nine-phenotype LDSC covariance structure, validates the resulting matrix dimensions and saves the covariance object and session information. It refuses to overwrite an existing output object.

The refactored script has passed a read-only HPC input check confirming the presence of all nine munged summary-statistic files and the European LD-score reference directory. The original nine-phenotype LDSC output is archived; the refactored script has not been used to rerun the full LDSC estimation.

## Interpretation of verification

The archived-object comparisons establish numerical consistency between the original July LDSC analyses and their six-region derivatives. The CAE-exclusion input check establishes that the refactored LDSC procedure can locate its required source data.

These checks do not independently reproduce all upstream GWAS processing or establish equivalence between the original and corrected GGE covariance specifications.

The multivariate GWAS execution, chromosome-level combination and downstream result verification are documented separately in `04_Multivariate_GWAS`.

## Archived-data handling

Original and corrected covariance objects, harmonised GWAS summary statistics, munged datasets, LD reference files and execution logs are maintained in controlled research storage.

Large intermediate files and restricted GWAS data are not included in this repository. The scripts accept locally supplied paths to the archived inputs, allowing the analysis and verification procedures to be executed without embedding institution-specific storage locations in the repository.

Generated verification reports and large analytical outputs may be retained locally rather than committed to GitHub.
