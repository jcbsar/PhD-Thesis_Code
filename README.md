# PhD Thesis - Analysis Code

This repository contains analytical code and supporting materials associated with a PhD thesis investigating epilepsy and functional dissociative seizures (FDS), with a particular focus on genetic architecture, associated risk factors, and neuroimaging characteristics. The repository is organised according to the chapters of the thesis, with code and supporting documentation provided for the respective analyses.

## Repository Structure

### Chapter 3: Genetic Architecture of Epilepsy and Subcortical Brain Volumes

This chapter investigates the genetic architecture of epilepsy phenotypes and subcortical brain volumes using linkage disequilibrium score regression (LDSC) and Genomic Structural Equation Modelling (Genomic SEM).

The analytical workflow includes:

1. **Standalone LDSC:** Estimation of SNP-based heritability and pairwise genetic correlations.
2. **LDSC within Genomic SEM:** Estimation of genetic covariance and sampling covariance matrices.
3. **Exploratory and Confirmatory Factor Analyses:** Investigation of latent genetic structures underlying epilepsy subtypes and subcortical brain volumes.
4. **Multivariate GWAS:** Identification of genetic variants associated with latent factors.
5. **Sensitivity Analyses:** Assessment of alternative model specifications and robustness of findings.
6. **Functional Annotation:** Identification and characterisation of associated genetic loci.

**[Chapter 3 – Standalone LDSC](Chapter_03_Genomic_SEM/01_standalone_LDSC/README.md)**

Additional analytical code is organised and documented within the corresponding chapter directories.

## Reproducibility

The repository contains R, Python, and shell scripts, together with relevant configuration files and summary-level analytical outputs. Where applicable, individual analysis directories provide information on software requirements, input data, analytical parameters, and execution procedures. The original GWAS summary statistics and external LD reference datasets are not included. These must be obtained from their respective sources, subject to applicable data-access requirements.

## Data Availability

This repository contains analytical code and selected summary-level results. Individual-level participant data and datasets subject to access restrictions are not included.

## Citation

When using or adapting code from this repository, please cite the associated PhD thesis and relevant methodological publications.
