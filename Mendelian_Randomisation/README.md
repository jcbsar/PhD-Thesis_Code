# Mendelian Randomisation

The study investigated the potential causal relevance of selected exposures to epilepsy susceptibility using complementary univariable, reverse-direction and multivariable MR approaches. This directory contains R scripts and reproducibility documentation for the Mendelian randomisation (MR) analyses.

## Repository structure

| Directory | Description |
|---|---|
| [`01_Univariable_MR`](01_Univariable_MR/) | Forward two-sample MR analyses of educational attainment, type 2 diabetes, stroke and major depressive disorder across ten epilepsy outcomes, including sensitivity analyses and pleiotropy diagnostics |
| [`02_Reverse_MR`](02_Reverse_MR/) | Reverse-direction MR analyses examining genetically proxied epilepsy liability in relation to educational attainment |
| [`03_Multivariable_MR`](03_Multivariable_MR/) | Exploratory multivariable MR models evaluating educational attainment and type 2 diabetes jointly, including instrument-flow validation and conditional instrument-strength sensitivity analyses |

Each directory contains its own README describing the analytical methods, input requirements, scripts, outputs and reproducibility validation.

## Analytical overview

### Forward univariable MR

Four exposures were investigated across ten epilepsy outcomes, comprising 40 primary comparisons. Inverse variance weighted (IVW) MR was the primary method, with MR-Egger regression, weighted median and MR-RAPS as sensitivity analyses. Additional diagnostics included Cochran's Q heterogeneity tests, MR-Egger intercept tests and MR-PRESSO.

### Reverse-direction MR

Reverse MR assessed associations of genetically proxied liability to all epilepsy, genetic generalised epilepsy (GGE) and juvenile myoclonic epilepsy (JME) with educational attainment. Other epilepsy phenotypes did not have sufficient independent instruments for the planned analyses.

### Multivariable MR

Exploratory multivariable MR evaluated associations of educational attainment and type 2 diabetes with six epilepsy outcomes, conditioning each exposure on the other. Reproducibility analyses included verification of archived instrument selection and harmonisation records, reproduction of MR estimates, and assessment of conditional F-statistics across alternative covariance assumptions.

## Reproducibility

The scripts were reconstructed and validated using archived analysis inputs and outputs from the original thesis analyses.

Primary IVW estimates were reproduced exactly for all 40 forward univariable comparisons and all three reverse-direction comparisons. Multivariable MR estimates and instrument-strength calculations were also validated against archived outputs.

Simulation-based methods, including weighted median and MR-PRESSO, may produce slightly different p-values between runs. These differences and the scope of validation are documented in the relevant subdirectories.

## Input data and availability

The analyses use previously prepared GWAS summary statistics, clumped genetic instruments and harmonised datasets. These input datasets are not distributed in this repository.

Scripts requiring local MR data use the `MR_DATA_DIR` environment variable. Additional file requirements and execution instructions are documented within each analytical directory.

The repository provides reproducible analysis code and verification procedures but does not independently reconstruct every upstream GWAS processing or LD-clumping step from raw summary statistics.

## Interpretation

The analyses should be interpreted in the context of the assumptions and limitations of MR, including genetic instrument validity, horizontal pleiotropy, population ancestry, statistical power and potential sample overlap.

Detailed study findings, interpretation and limitations are presented in Chapter 5 of the thesis.
