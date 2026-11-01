# Multivariable Mendelian Randomisation

## Overview

This directory contains the R scripts supporting the exploratory multivariable Mendelian randomisation (MVMR) analyses presented in Chapter 5 of the thesis *From genetic risk to brain structure: An integrative investigation of epilepsy and functional dissociative seizures*.

The analyses investigated the associations of genetically proxied educational attainment and type 2 diabetes (T2D) with six epilepsy outcomes: all epilepsy, focal epilepsy (FE), genetic generalised epilepsy (GGE), focal epilepsy with hippocampal sclerosis (F-HS), focal lesion-negative epilepsy (F-Neg), and juvenile myoclonic epilepsy (JME).

The analyses used jointly LD-clumped genetic instruments, followed by cross-exposure and outcome-specific harmonisation. MVMR estimates were obtained using inverse-variance-weighted regression without an intercept.

## Scripts

| Script | Description |
|---|---|
| `01_validate_instrument_flow.R` | Validates the archived instrument-selection and harmonisation records against the final instrument counts reported in Appendix Table A6. |
| `02_run_mvmr.R` | Reproduces MVMR estimates for all six epilepsy outcomes using the final archived harmonised datasets and validates these estimates against the original results. |
| `03_conditional_F_sensitivity.R` | Recalculates conditional F-statistics under zero covariance and across nine assumed correlations between SNP–exposure effect estimates, validating the results against Appendix Tables A9 and A10. |

## Instrument selection

Educational attainment and T2D instruments were initially selected separately using a genome-wide significance threshold of p < 5 × 10⁻⁸ and LD clumping with a 10,000 kb window, r² < 0.001, and a European ancestry reference population.

The final MVMR instrument-selection workflow was:

| Stage | Number of SNPs |
|---|---:|
| Educational attainment instruments | 418 |
| T2D instruments | 460 |
| Unique variants in the combined candidate set | 877 |
| Variants retained after joint LD clumping | 522 |
| Variants with valid associations in both exposure GWAS datasets | 519 |
| Variants retained after cross-exposure harmonisation | 503 |

Following harmonisation with the epilepsy outcome datasets, the final MVMR models included 351–364 SNPs, depending on the outcome.

## Reproducibility

The archived instrument-selection audits were checked for consistency in SNP identities, instrument counts, exposure associations and harmonisation exclusions.

The MVMR estimates for all six epilepsy outcomes were independently recalculated from the archived harmonised datasets and agreed with the original results, including the effect estimates, standard errors and p-values.

The conditional F-statistics and covariance sensitivity analyses were also recalculated and validated against the archived results.

**Scope of validation:** The instrument-flow script validates the recorded joint LD-clumping and harmonisation outputs. It does not independently rerun the original joint LD-clumping procedure from the source GWAS summary statistics.

## Input data

The scripts use files from the final archived MVMR analysis, including:

- Instrument-selection and harmonisation audit CSV files.
- Six outcome-specific `harmonised_MVMR_data_*.rds` datasets.
- The original MVMR results and conditional F-statistic reference tables.

The scripts require these archived inputs to be available locally. The complete source GWAS summary statistics and archived analysis datasets are not included in this directory.

## Running the analyses

The scripts were developed and validated using R 4.4.2, with `TwoSampleMR` 0.7.5 and `MVMR` 0.4.6.

Set the environment variable `MVMR_DATA_DIR` to the directory containing the final archived MVMR input files:

```r
Sys.setenv(
  MVMR_DATA_DIR = "path/to/MVMR_Education_T2D_JOINT_CLUMPED"
)
```

Run the scripts in numerical order:

```r
source("01_validate_instrument_flow.R")
source("02_run_mvmr.R")
source("03_conditional_F_sensitivity.R")
```

Run these commands from the directory containing the scripts, or provide their full file paths. The second and third scripts write reproduced results into an `outputs` directory relative to the current R working directory.

The scripts contain validation checks that stop execution if the recalculated results differ from the archived reference outputs.

## Interpretation

In the exploratory MVMR analyses, higher genetically proxied educational attainment was associated with lower odds of all epilepsy, FE and F-Neg after accounting for T2D liability, with these associations surviving correction for multiple testing. No T2D association survived correction.

Conditional F-statistics remained above the conventional threshold of 10 across the examined covariance assumptions.

These results should be interpreted in the context of the assumptions and limitations of Mendelian randomisation, as discussed in Chapter 5 of the thesis.
