# PhD Thesis - Analysis Code

This repository contains analytical code and supporting materials associated with a PhD thesis investigating epilepsy and functional dissociative seizures (FDS). The repository is organised by thesis chapter, with accompanying documentation describing the analytical methods, software requirements, input data and outputs.

## Repository Structure

### Chapter 3: Genetic Architecture of Epilepsy and Subcortical Brain Volumes

This chapter investigates genetic relationships between epilepsy phenotypes and subcortical brain volumes using linkage disequilibrium score regression (LDSC) and Genomic Structural Equation Modelling (Genomic SEM).

The analytical workflow comprises:

1. **Standalone LDSC:** Estimation of SNP-based heritability and pairwise genetic correlations between epilepsy phenotypes and subcortical brain volumes.
2. **LDSC within Genomic SEM:** Estimation of genetic covariance and sampling covariance matrices for subsequent multivariate analyses.
3. **Exploratory and Confirmatory Factor Analyses:** Investigation of latent genetic factors underlying epilepsy subtypes and brain volumes.
4. **Multivariate GWAS:** Identification of genetic variants associated with the identified latent factors.
5. **Sensitivity Analyses:** Evaluation of alternative model specifications and robustness of findings.
6. **Functional Annotation:** Characterisation and mapping of associated genetic loci.

Each section contains the relevant code and documentation.

## Reproducibility

Analysis scripts are accompanied, where applicable, by input specifications and summary-level outputs. Scripts use R, Python and shell-based workflows. The original GWAS summary statistics and external LD reference datasets are not distributed in this repository. These must be obtained from their respective sources, subject to applicable data-access requirements.

## Data Availability

Only code and appropriate summary-level analytical outputs are included. Individual-level participant data and restricted datasets are excluded.

## Citation

If using code from this repository, please cite the associated PhD thesis and relevant methodological software publications.
