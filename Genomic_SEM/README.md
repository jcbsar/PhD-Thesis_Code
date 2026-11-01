# Chapter 3 - Genetic Overlap Between Epilepsy and Subcortical Brain Volumes

This chapter investigates the genetic architecture of epilepsy phenotypes and subcortical brain volumes using linkage disequilibrium score regression (LDSC) and Genomic Structural Equation Modelling (Genomic SEM). The analyses aimed to characterise SNP-based heritability, examine genetic correlations between epilepsy subtypes and subcortical brain volumes, investigate underlying latent genetic structures, and identify genetic variants associated with the resulting latent factors.

## Analytical Workflow

The analyses are organised into the following sections:

1. **Standalone LDSC:** Estimation of SNP-based heritability and pairwise genetic correlations across epilepsy phenotypes and subcortical brain volumes.

2. **LDSC within Genomic SEM:** Estimation of genetic covariance matrices (S) and their corresponding sampling covariance matrices (V) for subsequent multivariate analyses.

3. **Exploratory and Confirmatory Factor Analyses (EFA/CFA):** Evaluation of latent genetic structures underlying epilepsy subtypes and subcortical brain volumes.

4. **Multivariate GWAS:** Identification of genetic variants associated with the latent factors using Genomic SEM.

5. **Sensitivity Analyses:** Assessment of alternative factor structures and model specifications.

6. **Functional Annotation:** Identification and annotation of associated genetic loci using FUMA.

7. **Figures and Tables:** Code and supporting outputs for the presentation of analytical results.

## LDSC Analyses

Two LDSC approaches were implemented for distinct analytical purposes.

**Standalone LDSC** was used to estimate SNP-based heritability and pairwise genetic correlations. These analyses included nine epilepsy phenotypes and seven subcortical brain-volume phenotypes for heritability estimation, with genetic correlations examined among seven epilepsy subtypes and seven brain volumes. Results are presented in Figure 3.1 and Appendix Tables A3-A4.

**LDSC within Genomic SEM** was used to estimate genetic covariance and sampling covariance matrices for the factor analyses and multivariate GWAS. The final latent brain-volume factor incorporated six subcortical regions, excluding the brainstem. The corresponding scripts, input specifications and analytical outputs are documented within the relevant subdirectories.

## Analysis Code

[01 - Standalone LDSC](01_standalone_LDSC/README.md)

Additional analytical sections are organised in their respective directories.

## Data and Reproducibility

The analyses used publicly available GWAS summary statistics and European-ancestry LD reference data. The original GWAS summary statistics and external reference datasets are not redistributed in this repository. Individual analysis directories describe the relevant software, input requirements, analytical parameters and execution procedures.
