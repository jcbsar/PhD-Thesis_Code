
# South African Pilot Neuroimaging Study

## Overview

This directory contains R scripts documenting the final neuroimaging analyses undertaken as part of a South African pilot study comparing structural brain measures among individuals with genetic generalised epilepsy (GGE), focal epilepsy (FE), functional dissociative seizures (FDS), and unaffected controls.

The study formed part of a doctoral thesis entitled *From genetic risk to brain structure: An integrative investigation of epilepsy and functional dissociative seizures*.

The analyses examined group differences in subcortical brain volumes, thalamic nuclei, hippocampal subfields, amygdala nuclei, cortical thickness, and cortical surface area. The objective was to investigate structural neuroanatomical characteristics across diagnostic groups in a South African cohort.

This repository documents the final 2026 analytical workflow rather than the preliminary analyses undertaken during 2025.

**Reproducibility status:** The scripts have been reconstructed and organised from the original August 2026 analytical source. Numerical equivalence with the final thesis results has not yet been established using the restricted final-study data.

## Software and processing

Structural MRI measures were derived using FreeSurfer version 7.4.1, with quality control informed by ENIGMA neuroimaging protocols.

Statistical analyses were undertaken in R version 4.4.2. The reconstructed scripts primarily use `dplyr` and `emmeans`, together with base R statistical functions.

The repository does not contain raw MRI scans, FreeSurfer segmentation files, participant-level measurements, or restricted clinical information.

## Study sample

The final analytical sample comprised 150 participants:

| Diagnostic group | N |
| --- | ---: |
| Unaffected controls | 61 |
| Genetic generalised epilepsy (GGE) | 21 |
| Focal epilepsy (FE) | 49 |
| Functional dissociative seizures (FDS) | 19 |
| **Total** | **150** |

Two diagnostic grouping specifications were examined:

1. **Four-group specification:** Control, GGE, FE and FDS, with each clinical group compared against controls.
2. **Collapsed specification:** Control, combined epilepsy (GGE and FE) and FDS, with comparisons of epilepsy versus controls, FDS versus controls, and epilepsy versus FDS.

## Statistical analyses

Group differences in regional brain measures were examined using linear regression models, followed by estimated marginal means contrasts.

Benjamini–Hochberg false discovery rate (BH-FDR) correction was applied separately within the specified anatomical outcome families and diagnostic contrasts.

| Analysis | Outcomes | Model covariates | BH-FDR family |
| --- | --- | --- | --- |
| Primary subcortical volumes | 9 bilateral or midline regions | Age, sex, ICV | 9 regions per contrast |
| Thalamic nuclei | 25 bilateral nuclei | Age, sex, ICV | 25 nuclei per contrast |
| Hippocampal subfields | 18 bilateral subfields | Age, sex, ICV | 18 subfields per contrast |
| Amygdala nuclei | 9 bilateral nuclei | Age, sex, ICV | 9 nuclei per contrast |
| Targeted cortical thickness | 12 regions per hemisphere | Age, sex, hemisphere-specific Euler number | 12 regions per hemisphere and contrast |
| Whole-cortex thickness | 34 regions per hemisphere | Age, sex, hemisphere-specific Euler number | 34 regions per hemisphere and contrast |
| Cortical surface area | 34 regions per hemisphere | Age, sex, ICV, hemisphere-specific Euler number | 34 regions per hemisphere and contrast |

ICV denotes intracranial volume.

Bilateral subcortical and subregional measures were calculated by summing the corresponding left and right hemisphere volumes, with brainstem volume treated as a midline measure.

### Sensitivity analyses

A log-transformed lateral ventricular volume analysis was undertaken to examine the influence of distributional characteristics on the reported group differences. The reconstructed code applies BH-FDR across the nine-region sensitivity family; equivalence with the original reported adjustment remains unverified.

Additional right-hemisphere cortical analyses were undertaken following the exclusion of one participant with an extreme right-hemisphere Euler number. This exclusion applied only to the specified sensitivity analyses, not the primary analytical cohort.

Region-specific cortical quality-control decisions were also incorporated into the final analytical dataset.

## Repository structure

```text
Pilot_Neuroimaging_Study/
|
|-- 01_Data_Preparation_QC/
|   |-- 01_prepare_final_dataset.R
|
|-- 02_Statistical_Analyses/
|   |-- 01_primary_subcortical.R
|   |-- 02_subregional_volumes.R
|   |-- 03_cortical_models.R
|
|-- 04_Verification/
|   |-- 01_verify_imaging_results.R
|   |-- README.md
|
|-- .gitignore
|-- README.md
```

### Data preparation and quality control

`01_Data_Preparation_QC/01_prepare_final_dataset.R`

Prepares the final analytical dataset from an authorised local input file, including:

- Diagnostic group definitions and final cohort checks.
- Primary bilateral subcortical volumes.
- Bilateral thalamic, hippocampal and amygdala subregional measures.
- Cortical quality-control adjustments specified in restricted project documentation.
- A private indicator for the right-hemisphere sensitivity exclusion.

The reconstructed preparation workflow requires an authoritative private anatomical variable dictionary and quality-control manifest. These must be reconciled with the original August 2026 workflow before numerical reproducibility can be established.

### Statistical analyses

`02_Statistical_Analyses/01_primary_subcortical.R`

Fits linear regression models for nine primary subcortical volumes across the four-group and collapsed diagnostic specifications. Includes the reconstructed lateral ventricular sensitivity analysis.

Corresponds to Appendix Table A1.

`02_Statistical_Analyses/02_subregional_volumes.R`

Fits regional linear models for 25 thalamic nuclei, 18 hippocampal subfields and nine amygdala nuclei, with BH-FDR correction performed separately for each anatomical family and diagnostic contrast.

Corresponds to Appendix Tables A2, A4 and A5.

`02_Statistical_Analyses/03_cortical_models.R`

Fits targeted cortical thickness, whole-cortex thickness and cortical surface-area models for both diagnostic specifications, including right-hemisphere sensitivity analyses.

Corresponds to Appendix Tables A3 and A6–A8.

### Verification

`04_Verification/01_verify_imaging_results.R`

Checks analytical output structure, expected comparison counts, selected numerical results and the reported numbers of FDR-significant subregional findings.

The verification script does not establish complete numerical equivalence across every comparison in Appendix Tables A1–A8.

See `04_Verification/README.md` for details.

## Restricted data and execution

The analytical scripts are intended for use only by authorised researchers with access to the approved study data and corresponding quality-control documentation.

The preparation script requires the following environment variables:

| Variable | Description |
| --- | --- |
| `ABHGP_INPUT_RDS` | Location of the authorised final-study input dataset |
| `ABHGP_PRIVATE_OUTPUT_DIR` | Private directory for prepared data and analysis outputs |
| `ABHGP_SUBREGION_DICTIONARY` | Private CSV mapping bilateral subregional measures to original FreeSurfer variables |
| `ABHGP_QC_MANIFEST_RDS` | Private RDS containing documented region-specific QC decisions and sensitivity exclusion information |

The preparation script writes `final_prepared_private.rds` to the authorised private output directory.

The primary, subregional and cortical analysis scripts use this prepared dataset. The verification script subsequently reads their aggregate statistical outputs.

Scripts should be executed in the following order:

1. `01_Data_Preparation_QC/01_prepare_final_dataset.R`
2. `02_Statistical_Analyses/01_primary_subcortical.R`
3. `02_Statistical_Analyses/02_subregional_volumes.R`
4. `02_Statistical_Analyses/03_cortical_models.R`
5. `04_Verification/01_verify_imaging_results.R`

The private dictionary and QC manifest are required supporting inputs, not generated by the public preparation script. The public repository therefore cannot be executed end-to-end without these authorised materials.

## Verification status and outstanding limitations

**Status: Reconstructed analytical code; numerical verification pending.**

The public scripts provide a structured implementation of the analytical methods reported in the final thesis, but they have not yet been executed and validated against the restricted final analytical dataset.

The following issues remain relevant:

- The source-to-bilateral mappings for the 52 subregional outcomes must be verified against the authoritative FreeSurfer variable definitions.
- The reconstructed QC manifest must accurately reproduce the original final region-specific exclusions and sensitivity analysis decision.
- The original calculation underlying the reported log-transformed lateral ventricular BH-FDR result requires confirmation.
- Appendix Table A4 mentions Euler number in the hippocampal subfield model note, whereas the thesis Methods and original analytical code specify adjustment for age, sex and ICV. The reported model should be confirmed against archived final analytical outputs.
- Selected numerical checks are implemented, but comprehensive numerical verification of every result in Appendix Tables A1–A8 has not been completed.

The original August 2026 analysis source is retained separately as a historical record and is not represented here as a fully standalone executable pipeline.

## Data availability and ethical considerations

The study was conducted using restricted clinical and structural MRI data from the African Brain Health Genomics Project.

Participant-level clinical information, MRI measurements, participant identifiers, individual quality-control decisions, and private data dictionaries are not included in this repository because of ethical and institutional restrictions.

The code is shared to support transparency regarding analytical methods and statistical procedures. Reproduction of the findings requires authorised access to the original data and documented quality-control workflow.

Aggregate statistical outputs should be reviewed for disclosure suitability before public release.
