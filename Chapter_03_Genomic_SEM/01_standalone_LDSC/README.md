# Standalone Linkage Disequilibrium Score Regression (LDSC)

## Overview

Standalone linkage disequilibrium score regression (LDSC) was performed to estimate SNP-based heritability and pairwise genetic correlations for epilepsy phenotypes and subcortical brain volumes, as presented in **Figure 3.1 and Appendix Tables A3-A4** of Chapter 3. The analyses were conducted using the standalone Python implementation of LDSC (v1.0.1), with European-ancestry reference LD scores.

Two sets of analyses were performed:

1. **Univariate LDSC:** Estimation of SNP-based heritability (h²) for nine epilepsy phenotypes and seven subcortical brain volumes.
2. **Bivariate LDSC:** Estimation of pairwise genetic correlations (rg) across seven epilepsy subtypes and seven subcortical brain volumes.

These analyses were conducted separately from the LDSC procedure implemented within Genomic SEM, which was subsequently used to estimate genetic covariance and sampling covariance matrices for the multivariate analyses.

## Phenotypes

**Epilepsy phenotypes**

- Genetic generalised epilepsy (GGE)
- Focal epilepsy (FE)
- Generalised tonic-clonic seizures alone (GTCSA)
- Childhood absence epilepsy (CAE)
- Juvenile absence epilepsy (JAE)
- Juvenile myoclonic epilepsy (JME)
- Focal epilepsy with hippocampal sclerosis (F-HS)
- Focal lesion-negative epilepsy (F-Neg)
- Focal epilepsy with other lesions (F-Oth)

**Subcortical brain volumes**

- Nucleus accumbens (ACC)
- Amygdala (AMY)
- Brainstem (BRS)
- Caudate (CAU)
- Pallidum (PAL)
- Putamen (PUT)
- Thalamus (THA)

GGE and FE were included in the univariate heritability analyses but excluded from the pairwise genetic correlation analyses.

## Analytical Methods

### SNP-based Heritability

Univariate LDSC was used to estimate SNP-based heritability for each phenotype.

Heritability estimates for binary epilepsy phenotypes were expressed on the liability scale using the sample and population prevalence specifications provided in `config/trait_prevalence_specification.tsv`.

A sample prevalence of 0.50 was specified for the epilepsy phenotypes, with phenotype-specific population prevalence estimates. Brain-volume phenotypes were analysed as quantitative traits without prevalence transformations. Estimates and standard errors were used to calculate 95% confidence intervals. For binary epilepsy phenotypes, effective sample sizes were used with a sample prevalence of 0.5 for liability-scale transformation, following recommendations for case–control GWAS meta-analyses. Phenotype-specific population prevalence values were specified based on the ILAE Consortium's 2018 GWAS (Supplementary Table 11). Quantitative brain-volume phenotypes did not require prevalence correction.

### Genetic Correlations

Bivariate LDSC was used to estimate pairwise genetic correlations between seven epilepsy subtypes and seven subcortical brain volumes. A total of 91 unique pairwise correlations were estimated across the 14 phenotypes, including 49 epilepsy–brain-volume comparisons. The Benjamini–Hochberg false discovery rate (BH-FDR) procedure was applied to the 49 epilepsy–brain-volume comparisons to account for multiple testing.

## Repository Contents

| File | Description |
|---|---|
| `config/trait_prevalence_specification.tsv` | GWAS input filenames and prevalence specifications |
| `data/figure_3_1A_heritability.csv` | Heritability estimates, standard errors and 95% confidence intervals |
| `data/figure_3_1B_rg_summary.tsv` | Pairwise genetic correlation estimates, standard errors, Z statistics and p-values |
| `data/figure_3_1B_rg_matrix.csv` | Genetic correlation matrix used for visualisation |
| `results/figure_3_1_cross_domain_fdr.csv` | BH-FDR-adjusted results for the 49 epilepsy–brain-volume comparisons |
| `scripts/01_run_h2.sh` | Standalone LDSC heritability analysis |
| `scripts/02_run_rg.sh` | Standalone LDSC genetic correlation analysis |
| `scripts/03_verify_results.py` | Validation of summary-level outputs and BH-FDR calculations |
| `scripts/04_plot_figure_3_1.R` | Generation of heritability and genetic correlation plots |

## Software Requirements

The LDSC analyses require:

- Standalone Python LDSC v1.0.1
- Python and the dependencies required by LDSC
- European-ancestry reference LD scores and regression weights (`eur_w_ld_chr`)
- Harmonised and munged GWAS summary statistics

The accompanying visualisation script requires R with the following packages:

- `ggplot2`
- `pheatmap`

The summary-level verification script uses the Python standard library and does not require GWAS summary statistics or LDSC.

## Running the Analyses

### 1. Input Data

The LDSC scripts require harmonised GWAS summary statistics in `.sumstats.gz` format.

The expected input filenames and prevalence parameters are defined in:

`config/trait_prevalence_specification.tsv`

GWAS summary statistics and external LD reference files are not included in this repository and must be obtained separately.

### 2. Heritability and Genetic Correlation Analyses

The shell scripts use environment variables to specify the locations of the LDSC installation, GWAS summary statistics, LD reference files and output directory.

From the `01_standalone_LDSC` directory:

```bash
export LDSC_PY=/path/to/ldsc.py
export MUNGED_DIR=/path/to/munged_sumstats
export EUR_LD_DIR=/path/to/eur_w_ld_chr
export ANALYSIS_OUT_DIR=/path/to/analysis_outputs
export PYTHON_BIN=python
```

Run the heritability analyses:

```bash
bash scripts/01_run_h2.sh
```

Run the genetic correlation analyses:

```bash
bash scripts/02_run_rg.sh
```

The scripts use the recorded LDSC command configurations and input specifications. They are portable shell implementations rather than copies of the original HPC job-submission scripts.

**Note:** The heritability execution manifest covers 14 phenotypes. Although broad GGE and FE estimates are included in the final Figure 3.1A data, their original standalone execution commands and logs are not available in the accompanying analysis records. These two estimates therefore cannot be independently regenerated using the supplied shell script alone.

The shell scripts execute LDSC commands directly and do not submit jobs through SLURM.

### 3. Verification of Summary-Level Results

The supplied Python script checks the consistency of the stored heritability and genetic correlation results and recalculates BH-FDR-adjusted p-values.

```bash
python3 scripts/03_verify_results.py
```

Where the original standalone LDSC logs are available, the exported estimates can additionally be compared against them:

```bash
python3 scripts/03_verify_results.py \
  --h2-logs /path/to/h2_logs \
  --rg-logs /path/to/rg_logs
```

The verification procedure evaluates:

- The 16 heritability entries and corresponding confidence intervals.
- The completeness and consistency of the 91 pairwise genetic correlation estimates.
- The 49 epilepsy–brain-volume comparisons.
- The BH-FDR correction for multiple testing.

Comparison with available original execution logs supports the 14 subtype/brain-volume heritability estimates and all 91 pairwise genetic correlation estimates. The broad GGE and FE heritability estimates are retained from the final figure-data records.

### 4. Generating Figure 3.1

The accompanying R script generates the heritability plot and genetic correlation heatmap from the supplied summary-level results.

Install the required R packages:

```r
install.packages(c("ggplot2", "pheatmap"))
```

Run:

```bash
Rscript scripts/04_plot_figure_3_1.R
```

The script generates PNG and PDF outputs in:

`outputs/figures/`

The figure-generation script uses the stored LDSC estimates and does not independently rerun the GWAS analyses.

## Summary of Findings

The standalone LDSC analyses provided estimates of SNP-based heritability across the epilepsy and subcortical brain-volume phenotypes.

Four epilepsy–brain-volume genetic correlations had unadjusted p-values below 0.05:

- GTCSA-amygdala
- JAE-caudate
- JME-brainstem
- JME-pallidum

None of the 49 epilepsy-brain-volume genetic correlations remained statistically significant following BH-FDR correction.

Detailed estimates are provided in the accompanying data files and in Figure 3.1 and Appendix Tables A3-A4 of the thesis.

## Data Availability

The repository contains analysis scripts, input specifications and selected summary-level outputs. Original GWAS summary statistics and external LD reference datasets are not redistributed.

The standalone LDSC estimates presented here are distinct from the genetic covariance and sampling covariance matrices generated using `GenomicSEM::ldsc()` for the subsequent factor analyses.
