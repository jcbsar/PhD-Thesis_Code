
# 04 - Multivariate Genome-Wide Association Studies

## Purpose

This directory documents the multivariate genome-wide association studies (GWAS) conducted using Genomic Structural Equation Modelling (Genomic SEM) for Chapter 3 of the thesis.

The analyses evaluated SNP associations with latent factors representing genetic liability to generalised genetic epilepsy (GGE), focal epilepsy (FE), and subcortical brain volumes. They were conducted using chromosome-partitioned summary statistics on the University of Cape Town high-performance computing (HPC) cluster.

This directory currently preserves the original execution scripts for the corrected GGE–brain GWAS undertaken in October 2026. Documentation of the FE analysis, CAE-exclusion sensitivity analysis, and upstream summary-statistic preparation is still being completed.

## Analysis overview

The retained GGE model included:

- GGE indicators: CAE, GTCS, JAE, JME
- Brain indicators: accumbens, amygdala, caudate, pallidum, putamen, thalamus
- Correlated GGE and brain-volume factors
- CAE residual variance constrained to zero

The multivariate GWAS estimated SNP associations with both latent factors simultaneously.

The corrected October analysis used:

- `GenomicSEM::userGWAS()`
- Diagonally weighted least squares (DWLS)
- Fixed measurement parameters (`fix_measurement = TRUE`)
- SNP heterogeneity testing (`Q_SNP = TRUE`)
- Smoothing diagnostics (`smooth_check = TRUE`)
- Genome-wide significance threshold: p < 5 × 10⁻⁸

The corrected LDSC covariance input was:

`ldsc_GGE_brain_6region_corrected.rds`

The chromosome-level SNP input files were reused from the harmonised GGE–brain GWAS dataset.

## Corrected GGE–brain GWAS workflow

The following scripts preserve the original HPC filenames and settings.

| Script | Function |
|---|---|
| `00_preflight_GGE_corrected.sh` | Validate corrected LDSC input, chromosome files, and output directory |
| `02_run_userGWAS_chr_corrected.R` | Run multivariate GWAS for one chromosome |
| `02_GGE_corrected_array.sh` | Submit chromosome jobs through SLURM |
| `06_stream_combine_one_path_corrected.R` | Combine and validate chromosome-level results |
| `06_stream_combine_corrected_array.sh` | Original SLURM launcher for combining the two latent-factor outputs |
| `06_combine_corrected_brain_only.sh` | Separate launcher used for the completed brain-factor combination |

The numbering reflects the original HPC workflow rather than a continuous sequence of GitHub scripts.

### Input validation

The original preflight script checked:

- Corrected genetic covariance matrix dimensions: 10 × 10
- LDSC intercept matrix dimensions: 10 × 10
- Sampling covariance matrix dimensions: 55 × 55
- CAE genetic variance approximately 0.2404654
- Corrected JME genetic variance approximately 0.2534618
- Existence of all 22 chromosome-level SNP inputs
- Availability of an empty corrected-results directory

These checks are implemented in the archived script. The original preflight execution log has not yet been independently reviewed.

### HPC execution

The GWAS was partitioned by autosomal chromosome using a SLURM array.

Each chromosome job was configured with:

- Four CPU cores
- 48 GB memory
- A 24-hour time limit
- A maximum of four concurrent chromosome jobs

The chromosome-level R script saved both latent-factor association results and associated metadata.

### Combination and quality control

Chromosome-level outputs were combined separately for the GGE and brain-volume factors.

The combination script checked:

- Expected chromosome coverage and row counts
- Required output columns
- Chromosome assignment
- Missing SNP identifiers
- Duplicate SNP identifiers within chromosomes
- Reported model errors
- Genome-wide significant SNP counts
- Maximum absolute smoothing statistic
- Distinct warning messages

Results were written sequentially to a temporary CSV file before being renamed to the final output.

## Verified corrected results

The final October 2026 outputs were inspected against their validation reports and HPC combination logs.

| Statistic | GGE factor | Brain-volume factor |
|---|---:|---:|
| Variants analysed | 4,660,537 | 4,660,537 |
| Chromosomes | 22 | 22 |
| Genome-wide significant SNPs | 249 | 608 |
| Duplicate SNP identifiers reported | 0 | 0 |
| Missing SNP identifiers | 0 | 0 |
| Reported model errors | 0 | 0 |
| Maximum absolute Z_smooth | 0.185742 | 0.185742 |

The corrected GGE result supersedes the original August output containing 252 genome-wide significant SNPs.

The chromosome-level validation reports recorded a recurring `lavaan` warning concerning an override of the `se` option. These warnings were retained in the validation files.

## Execution history

The corrected GGE combination completed at the R-script level on 23 October 2026, producing the validated 249-SNP result. Its SLURM wrapper subsequently reported a nonzero exit status.

The initial brain-factor combination attempt on 26 October did not complete. A subsequent job, SLURM ID `1534576`, completed successfully and generated the validated 608-SNP brain-factor output.

The separate brain-factor combination script intentionally omitted gzip compression.

The final archived results are uncompressed CSV files.

## Data preservation

The corrected SNP-level GWAS files and associated validation reports were preserved separately from the GitHub repository.

SHA-256 checksums of the original HPC outputs:

**GGE factor**

`b5e2bc14e6b0524181fd25299e580681c34a1e53cb82fa9472fdd9db3134e4d5`

**Brain-volume factor**

`2e20f2ecded4eaae875dd683e7a363685cfbf7ff1f48254680af422bff56402b`

The complete SNP-level outputs are not committed to GitHub because of their size.

## Reproducibility considerations

These scripts are preserved as historical execution records.

They contain HPC-specific paths, SLURM settings, and assumptions about the installed R environment. They should not be considered portable without adaptation.

In particular:

- The corrected GWAS script requires pre-existing harmonised chromosome-level SNP inputs and the corrected LDSC covariance object.
- The original combination R script uses `data.table` functions without explicitly attaching the package.
- The original combination R script removes pre-existing destination CSV files before writing new results.
- Duplicate SNP checking is performed within each chromosome; no separate global uniqueness check has yet been documented.
- The original scripts should not be rerun against the existing completed HPC output directories.

Additional portable reproduction and validation scripts will be distinguished from these historical scripts.

## Remaining reproducibility work

The following components are not yet fully documented in this directory:

1. Upstream GWAS harmonisation and chromosome partitioning
2. FE–brain multivariate GWAS
3. GGE analysis excluding CAE
4. FE-based brain-factor sensitivity analysis
5. Independent final-output QC and comparison against thesis results

Functional mapping and annotation using FUMA are documented separately in `05_FUMA_Locus_Mapping`.

## Thesis reference

Chapter 3: Genetic overlap between epilepsy and subcortical brain volumes using Genomic Structural Equation Modelling.

Relevant methods: Sections 3.3.4–3.3.5.

Relevant results: Section 3.4.3.
