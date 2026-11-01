# Reverse Mendelian Randomisation

## Overview

This directory contains R scripts for reproducing the reverse-direction two-sample Mendelian randomisation (MR) analyses presented in Chapter 5 of the PhD thesis, *From genetic risk to brain structure: An integrative investigation of epilepsy and functional dissociative seizures*.

The analyses investigated whether genetically proxied liability to epilepsy was associated with educational attainment, complementing the forward MR analyses examining educational attainment as an exposure.

Reverse MR was feasible for three epilepsy phenotypes:

- All epilepsy
- Genetic generalised epilepsy (GGE)
- Juvenile myoclonic epilepsy (JME)

The remaining epilepsy phenotypes lacked sufficient independent genome-wide significant instruments for the planned multi-instrument analyses. Childhood absence epilepsy (CAE) had one independent instrument and was not included.

## Scripts

| Script | Description |
|---|---|
| `01_run_reverse_mr.R` | Performs SNP matching, harmonisation, MR estimation, Cochran's Q heterogeneity tests and MR-Egger intercept tests |
| `02_run_reverse_mr_presso.R` | Reproduces MR-PRESSO global, outlier and distortion analyses for the two eligible reverse MR comparisons |

## Input datasets

Epilepsy instruments were derived from the ILAE 2023 epilepsy GWAS and previously clumped using a genome-wide significance threshold of p < 5 × 10⁻⁸.

Educational attainment was treated as the outcome, using the full SSGAC 2022 European-ancestry GWAS dataset, excluding 23andMe.

| Epilepsy exposure | Clumped SNPs | SNPs retained after harmonisation |
|---|---:|---:|
| All epilepsy | 4 | 4 |
| GGE | 23 | 19 |
| JME | 3 | 3 |

Four GGE variants were removed during harmonisation because they were palindromic with intermediate allele frequencies.

The scripts use archived epilepsy instrument RDS files and the full educational attainment GWAS summary statistics. The original data are not distributed in this repository.

## MR methods

The following methods were implemented using `TwoSampleMR`:

- Inverse variance weighted (IVW) MR as the primary analysis
- MR-Egger regression
- Weighted median
- Robust adjusted profile score (MR-RAPS)

Exposure and outcome alleles were harmonised using `TwoSampleMR::harmonise_data(action = 2)`.

Cochran's Q statistics were calculated for IVW and MR-Egger regression, and MR-Egger intercept tests were performed to assess evidence of directional horizontal pleiotropy.

Reverse MR estimates are reported as beta coefficients on the educational attainment scale, rather than odds ratios.

The Bonferroni-corrected significance threshold for the three primary reverse MR comparisons was p < 0.0167.

## MR-PRESSO

MR-PRESSO was performed for all epilepsy and GGE using the original settings:

| Parameter | Setting |
|---|---|
| `OUTLIERtest` | `TRUE` |
| `DISTORTIONtest` | `TRUE` |
| `NbDistribution` | `10000` |
| `SignifThreshold` | `0.05` |

MR-PRESSO was not performed for JME because only three instruments were available after harmonisation.

The reconstructed script uses fixed random seeds. The original historical simulation seeds were not documented, so simulated p-values may differ slightly between the archived and reconstructed analyses.

## Running the scripts

Set the `MR_DATA_DIR` environment variable to the directory containing the archived MR datasets:

`Sys.setenv(MR_DATA_DIR = "/path/to/MR_data")`

Execute `01_run_reverse_mr.R` followed by `02_run_reverse_mr_presso.R`.

The first script requires the archived reverse MR results CSV to perform its built-in IVW validation. The MR-PRESSO script additionally checks its observed RSS statistics against the original archived values.

## Outputs

The scripts save their outputs under `Results/Reverse_MR_Reproduction/`.

| Output file | Description |
|---|---|
| `01_reverse_MR_results.csv` | Reverse MR estimates and beta-scale confidence intervals |
| `02_reverse_heterogeneity_results.csv` | IVW and MR-Egger Cochran's Q results |
| `03_reverse_egger_intercept_results.csv` | MR-Egger intercept tests |
| `04_reverse_instrument_flow.csv` | SNP counts before and after harmonisation |
| `05_reverse_sessionInfo.txt` | Software environment information |

The `MR_PRESSO/` subdirectory contains the reproduced global-test CSV, session information and complete saved MR-PRESSO result objects.

## Reproducibility validation

The reconstructed analyses were compared against the original archived results.

### Primary IVW results

| Exposure | SNPs | β | SE | p-value |
|---|---:|---:|---:|---:|
| All epilepsy | 4 | −0.08109 | 0.17262 | 0.6385 |
| GGE | 19 | −0.03181 | 0.02270 | 0.1611 |
| JME | 3 | 0.02120 | 0.01897 | 0.2638 |

All three IVW estimates, standard errors and p-values reproduced exactly.

MR-Egger and MR-RAPS estimates also agreed with the archived analyses. Weighted-median effect estimates reproduced exactly, while bootstrap-derived standard errors and p-values showed minor variation.

### Diagnostic validation

All six Cochran's Q results and all three MR-Egger intercept tests reproduced exactly, with zero numerical differences from the archived results.

Both eligible reverse MR-PRESSO global tests reproduced the original observed residual sums of squares:

| Exposure | RSSobs | Reproduced global p-value |
|---|---:|---|
| All epilepsy | 252.1128 | < 0.0001 |
| GGE | 220.4797 | < 0.0001 |

The complete MR-PRESSO result objects were preserved. Individual outlier and distortion estimates were not independently validated against every historical estimate.

## Interpretation and limitations

None of the three primary IVW reverse MR associations met the corrected significance threshold.

Significant heterogeneity was observed for all three epilepsy exposures, while MR-Egger intercept tests provided no statistically significant evidence of directional horizontal pleiotropy.

The small number of instruments, particularly for all epilepsy and JME, limits precision and the ability to assess certain MR assumptions. Consequently, the absence of statistically significant reverse associations should not be interpreted as proof of no causal effect.

These scripts reproduce the analyses from archived instruments and GWAS data. They do not independently repeat the original GWAS processing or LD clumping from raw summary statistics.
