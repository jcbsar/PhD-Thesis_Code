# Univariable Mendelian Randomisation

This directory contains R scripts for reproducing the forward two-sample univariable Mendelian randomisation (MR) analyses presented in Chapter 5 of the PhD thesis, *From genetic risk to brain structure: An integrative investigation of epilepsy and functional dissociative seizures*.

## Overview

The analyses investigated associations between genetically proxied exposures and epilepsy outcomes using GWAS summary statistics.

Four exposures were investigated: educational attainment, type 2 diabetes (T2D), any ischaemic stroke, and major depressive disorder (MDD).

Ten epilepsy outcomes were investigated: all epilepsy, focal epilepsy (FE), genetic generalised epilepsy (GGE), focal hippocampal sclerosis epilepsy (F-HS), focal lesion-negative epilepsy (F-Neg), focal other lesion epilepsy (F-Oth), generalised tonic-clonic seizures alone (GTCSA), juvenile absence epilepsy (JAE), juvenile myoclonic epilepsy (JME), and childhood absence epilepsy (CAE).

The analysis comprised **40 exposure–outcome comparisons**.

## Scripts

| Script | Description |
|---|---|
| `01_run_univariable_mr.R` | Performs SNP matching, allele harmonisation and four MR estimation methods |
| `02_run_heterogeneity_egger.R` | Calculates IVW and MR-Egger Cochran's Q heterogeneity statistics and MR-Egger intercept tests |
| `03_run_mr_presso.R` | Performs MR-PRESSO global, outlier and distortion tests and preserves the complete results |

### MR estimation methods

Inverse variance weighted (IVW) MR was the primary method. MR-Egger regression, weighted median and robust adjusted profile score (MR-RAPS) were used as sensitivity analyses.

The primary IVW estimates were evaluated using a Bonferroni-corrected threshold of **p < 0.00125**, accounting for 40 exposure–outcome comparisons.

### Diagnostic analyses

Cochran's Q was calculated for IVW and MR-Egger estimates. MR-Egger intercept tests were used to investigate evidence of directional horizontal pleiotropy.

MR-PRESSO was implemented using the original analytical settings:

- `OUTLIERtest = TRUE`
- `DISTORTIONtest = TRUE`
- `NbDistribution = 10000`
- `SignifThreshold = 0.05`

MR-PRESSO analyses required at least four retained instruments. The reconstructed script introduces fixed, comparison-specific random seeds to support reproducible reruns. The original historical MR-PRESSO simulations did not have a documented fixed seed.

MR-PRESSO output components are retained when available. The absence of an outlier or distortion component does not necessarily indicate an unsuccessful analysis.

## Input datasets

The scripts use previously prepared, LD-clumped exposure instruments and cleaned ILAE 2023 epilepsy GWAS outcome datasets, stored locally as RDS files.

| Exposure | Archived clumped instruments |
|---|---:|
| Educational attainment | 418 |
| T2D | 460 |
| Stroke | 23 |
| MDD | 147 |

Outcome-specific instrument counts were determined following SNP matching and allele harmonisation using `TwoSampleMR::harmonise_data(action = 2)`.

The original input data are not distributed in this repository. The scripts retain the archived directory structure and input filenames to support verification against the original analyses.

**Reproducibility scope:** These scripts reproduce MR estimation and diagnostic analyses using the archived clumped instruments and cleaned outcome datasets. They do not independently reproduce raw GWAS processing, instrument selection or LD clumping.

## Software requirements

The analyses were performed in R using `TwoSampleMR` and `MRPRESSO`.

Session information is saved by the scripts to record the software environment used for each run.

## Running the analyses

Before executing the scripts, configure the `MR_DATA_DIR` environment variable in R:

```r
Sys.setenv(MR_DATA_DIR = "/path/to/MR_data")
```

Execute the scripts in numerical order.

The MR-PRESSO script supports an optional exposure selection:

```r
Sys.setenv(MR_PRESSO_EXPOSURES = "MDD")
```

To run all four exposures:

```r
Sys.unsetenv("MR_PRESSO_EXPOSURES")
```

The MR-PRESSO script saves individual completed analyses as RDS files. Compatible saved results are reused on subsequent runs, allowing interrupted analyses to resume without unnecessarily repeating completed computations.

## Outputs

The first two scripts save results to:

`Results/Univariable_MR_Reproduction/`

| Output | Description |
|---|---|
| `01_univariable_MR_results.csv` | MR effect estimates, SEs, p-values, odds ratios and confidence intervals |
| `02_univariable_MR_run_log.csv` | Instrument-flow records |
| `03_sessionInfo.txt` | MR estimation software environment |
| `04_heterogeneity_results.csv` | IVW and MR-Egger Cochran's Q results |
| `05_egger_intercept_results.csv` | MR-Egger intercept estimates |
| `06_diagnostics_run_log.csv` | Diagnostic analysis completion records |
| `07_diagnostics_sessionInfo.txt` | Diagnostic analysis software environment |

MR-PRESSO outputs are saved under `Results/Univariable_MR_Reproduction/MR_PRESSO/`, including `08_MR_PRESSO_global_tests.csv`, `09_MR_PRESSO_sessionInfo.txt`, and individual full MR-PRESSO RDS results in `raw_results/`.

These generated datasets and analysis objects are retained locally rather than committed to the repository.

## Reproducibility validation

### MR estimates

The reconstructed workflow was executed using the original archived inputs and compared with the historical MR results.

| Method | Validation outcome |
|---|---|
| IVW | All 40 results reproduced exactly |
| MR-Egger | All 40 results reproduced exactly |
| MR-RAPS | All 40 results reproduced to near-numerical precision |
| Weighted median | All 40 effect estimates reproduced exactly; bootstrap SEs and p-values varied |

Weighted-median bootstrap variability caused two secondary significance classifications to change between the archived and reconstructed runs: educational attainment–GGE at p = 0.00125 and T2D–GTCSA at p = 0.05. The primary IVW results were unaffected.

Five educational attainment associations met the primary IVW Bonferroni threshold. No T2D, stroke or MDD association met this threshold.

### Heterogeneity and MR-Egger intercept tests

The reconstructed analyses generated 80 Cochran's Q estimates and 40 MR-Egger intercept results.

These were compared against the thesis supplementary results. The analyses agreed with the archived outputs, with minor differences identified in the rounded values reported in Appendix Table A4.

### MR-PRESSO

All 40 MR-PRESSO analyses completed successfully.

The reconstructed global tests reproduced the archived observed residual sums of squares (RSSobs) at the precision examined. Simulation-derived p-values showed expected variation between runs.

One comparison, stroke–focal other lesion epilepsy, crossed the nominal p = 0.05 threshold between the archived and reconstructed simulations. The original archived values remain the reference for thesis reporting.

The saved MR-PRESSO objects were audited:

| Component | Available results |
|---|---:|
| Global tests | 40/40 |
| Main MR estimates | 40/40 |
| Outlier-test components | 22/40 |
| Distortion-test components | 22/40 |

The presence of an outlier-test component does not imply that a statistically significant individual outlier was identified. Detailed outlier-corrected and distortion estimates were preserved in the full analysis objects but were not independently validated against every historical estimate.

## Interpretation and limitations

MR findings should be interpreted in the context of instrument validity, horizontal pleiotropy, statistical power, population ancestry and the methodological limitations discussed in Chapter 5.

The reconstructed analyses provide a reproducibility record of the original forward MR workflow. They should not be interpreted as independent replication using new GWAS datasets.
