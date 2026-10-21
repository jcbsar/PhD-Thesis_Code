# Chapter 3 — Genetic overlap between epilepsy and subcortical brain volumes

Chapter 3 combines **two LDSC workflows**:

1. **Standalone Python LDSC** to estimate univariate SNP heritability and pairwise genetic correlations, reported in thesis Figure 3.1 and Appendix Tables A3–A4.
2. **Multivariable LDSC within Genomic SEM** to estimate genetic covariance (S) and sampling covariance (V), subsequently used for exploratory/confirmatory factor models and multivariate GWAS.

These workflows are distinct and their estimates should not be silently mixed. The primary final six-region brain-volume factor excludes brainstem, whereas the descriptive Figure 3.1 correlation analyses include seven brain regions.

**Present release:** [`01_standalone_LDSC/`](01_standalone_LDSC/README.md). Other subfolders will be populated only after the original code, model configurations, and thesis estimates have been cross-checked.
