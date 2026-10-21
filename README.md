# Epilepsy PhD — analysis code

Code and supporting **summary-level** outputs for a PhD thesis examining epilepsy, genetic architecture, and subcortical brain volumes. The repository is organised by thesis chapter and will be populated incrementally after checking the research records.

## Status

| Thesis component | Status |
| --- | --- |
| Chapter 3: standalone LDSC; Figure 3.1A/B; Tables A3/A4 | Numerical outputs verified; scripts for rerunning LDSC reconstructed from command logs; figure code adapted from recovered source |
| Chapter 3: LDSC within Genomic SEM | Under review — not yet released here |
| Chapter 3: EFA/CFA, multivariate GWAS, sensitivity analyses, FUMA | Under review — not yet released here |
| Other thesis chapters | To be added after verification |

**Important:** Matching an archived result to an original log does not mean the complete analysis was rerun while preparing this repository. The original GWAS summary statistics and LD-reference files are not included. The reconstructed scripts are explicitly marked as such.

## Start here

See [Chapter 3 — standalone LDSC](Chapter_03_Genomic_SEM/01_standalone_LDSC/README.md). That folder includes verified published figure data, an FDR verification script, and code to reproduce Figure 3.1 from those data.

## Responsible sharing

- Do not commit participant-level data, restricted GWAS files, data-use agreements, account credentials, home-directory paths, or private HPC logs.
- Confirm permission to redistribute any source data before adding them, even if they are GWAS summary statistics.
- Run scripts using locally available reference datasets and working directories outside this Git repository.
- If you use or adapt this code in research, cite the related thesis/manuscript and methodological software. Formal citation metadata will be added after thesis details are confirmed.

## License

No reuse license has been selected yet. The author should choose one before broad code reuse or public release.
