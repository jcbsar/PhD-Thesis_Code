# Standalone LD score regression — Chapter 3

This folder supports thesis Figure 3.1 and Appendix Tables A3–A4. Original standalone Python LDSC v1.0.1 logs date from **11 August 2026**. They were checked against archived, summary-level figure exports. The local original logs were *not* committed because they contain user-specific HPC paths.

## Contents

| File | Purpose | Origin/status |
| --- | --- | --- |
| `config/trait_prevalence_specification.tsv` | Ordered 14-trait file names; sample and population prevalences | Archived input specification, unchanged |
| `data/figure_3_1A_heritability.csv` | 16 phenotypes: h², SE, 95% CI | Archived final figure data, unchanged |
| `data/figure_3_1B_rg_summary.tsv` | 91 pairwise LDSC rg, SE, Z and p | Archived results, unchanged |
| `data/figure_3_1B_rg_matrix.csv` | 14×14 rg matrix for heatmap | Archived results, unchanged |
| `results/figure_3_1_cross_domain_fdr.csv` | 49 cross-domain BH-FDR q-values | Regenerated from original archived p-values |
| `scripts/01_run_h2.sh` | Run each phenotype's `--h2` on HPC/Linux | **Reconstructed** from saved standalone command logs |
| `scripts/02_run_rg.sh` | Run 13 ordered `--rg` jobs for 91 unique pairs | **Reconstructed** from saved standalone command logs |
| `scripts/03_verify_results.py` | Validate saved figure outputs, recompute BH-FDR; optional comparison to local original logs | New independent validation code |
| `scripts/04_plot_figure_3_1.R` | Plot Figure 3.1A/B and combined Figure 3.1 using archived results | **Adapted/reconstructed**, not a verbatim final executed script |

## Reproduce the verification (no GWAS data required)

From this folder:

```bash
python3 scripts/03_verify_results.py
```

Expected: 16 heritability entries, 91 correlations, 49 cross-domain tests, **four unadjusted p < 0.05, zero BH-FDR q < 0.05**. The script writes its BH results to `results/figure_3_1_cross_domain_fdr.csv`.

If original log files are available **privately** on your computer, use:

```bash
python3 scripts/03_verify_results.py \
  --h2-logs /path/to/private/h2-logs \
  --rg-logs /path/to/private/rg-logs
```

The 14 archived univariate logs and 13 correlation logs available during preparation pass these checks. **The original broad GGE and FE univariate logs were not recovered.** Their h² estimates are preserved from the archived final figure-data CSV, and must not be represented as independently checked against executed original log files.

## Reproduce Figure 3.1 from archived numbers

Requires **R**, `ggplot2` and `pheatmap`:

```r
install.packages(c("ggplot2", "pheatmap"))
```

Then from this folder run:

```bash
Rscript scripts/04_plot_figure_3_1.R
```

The script creates PNG and PDF versions of Panels A and B and the combined figure under the ignored `outputs/figures/` directory. The numerical values and 95% CIs are read from the recovered export (not recalculated from GWAS). Heatmap colours are capped at −1 and 1, but cell text retains unconstrained rg values above 1. GTCS in source data is labelled **GTCSA** in plots.

**Environment limit:** The R plot script was reviewed but was **not executed** while assembling this repository because R was not available in the preparation environment. Verify its appearance on a computer with R installed before calling the figure regenerated.

## Rerun standalone LDSC from harmonised summary statistics (advanced)

This step requires the original HapMap3-aligned, munged `.hm3.sumstats.gz` files, an installation of standalone Python LDSC v1.0.1 (or a compatible pinned build), and the European reference and regression LD scores (`eur_w_ld_chr`). Those reference/data files are **not** redistributed here.

The original analyses used effective per-SNP sample sizes for binary epilepsy traits; the sample prevalence is **0.50** and population prevalences are given in `config/trait_prevalence_specification.tsv`. Brain volumes use `nan` for both prevalences. The 14-trait manifest excludes broad GGE and FE because the archived original `--h2` input commands for them have not been recovered.

Set the paths in your environment **outside the Git repository**, for example:

```bash
export LDSC_PY=/path/to/ldsc.py
export MUNGED_DIR=/path/to/private/munged
export EUR_LD_DIR=/path/to/reference/eur_w_ld_chr
export ANALYSIS_OUT_DIR=/path/to/private/analysis_outputs
export PYTHON_BIN=python

bash scripts/01_run_h2.sh
bash scripts/02_run_rg.sh
```

Run within your institution's permitted computing environment, with suitable SLURM allocations if required. The scripts do not automatically submit `sbatch` jobs and make **no** promises about memory or walltime. They reproduce the observed command configurations but were not rerun with the full GWAS input files during repository preparation.

## Findings checked

- Broad GGE: h² = **0.2317**, SE = **0.0178**, taken from recovered final Figure 3.1A data; broad FE: h² = **0.0287**, SE = **0.0064**, same provenance.
- Original 14 available `--h2` logs match their figure-data values to reported precision.
- All 91 archived correlations matched the 13 August 2026 `--rg` logs (rg, SE, Z and p).
- Four nominal epilepsy/brain correlations: GTCSA–amygdala, JAE–caudate, JME–brainstem and JME–pallidum. All have BH-FDR q ≈ 0.52185; none is significant after correcting across 49 cross-domain tests.

## Scope and provenance

Standalone LDSC outputs should **not** be substituted for `GenomicSEM::ldsc()` S/V matrices; those were separate analyses and will be documented under `../02_GenomicSEM_LDSC/`. The preliminary "FINAL CORRECTED STANDALONE LDSC + FIGURE 1A + FIGURE 1B.R" contained older **GGE = 0.584** and **FE = 0.067** plot values, so that original working script is deliberately **not** presented as the final executed Figure 3.1 workflow. The reported final Figure 3.1A CSV is used instead.

For exact original processing logs, consult the author's private research archive. Do not place logs exposing original HPC user account paths or restricted summary-statistic files into a public repository without review.
