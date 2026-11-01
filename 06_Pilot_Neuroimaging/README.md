# Chapter 6 — South African pilot neuroimaging study

## Status
**Reconstruction in progress; not yet validated against restricted final-study data.**

This directory documents the **2026 final-study analyses**, not the 2025 preliminary testing scripts. Source analyses were conducted in R 4.4.2, using measures derived from FreeSurfer 7.4.1. The repository is intended to provide transparent analytical methodology, not an independently executable copy of confidential project data.

## Cohort and comparisons
Final analytical cohort: 150 participants (controls 61; GGE 21; FE 49; FDS 19). Four-level analyses compare each clinical group to controls. Collapsed analyses additionally compare combined epilepsy to controls and to FDS.

## Analysis families (as reported in the final thesis)
| Domain | Outcomes | Covariates | BH-FDR family |
| --- | --- | --- | --- |
| Primary subcortical volumes | Nine bilateral/midline ROIs | Diagnosis, age, sex, ICV | Nine ROIs per contrast |
| Thalamic nuclei | 25 bilateral nuclei | Diagnosis, age, sex, ICV | 25 nuclei per contrast |
| Hippocampal subfields | 18 bilateral subfields | Diagnosis, age, sex, ICV | 18 subfields per contrast |
| Amygdala nuclei | Nine bilateral nuclei | Diagnosis, age, sex, ICV | Nine nuclei per contrast |
| Targeted cortical thickness | 12 ROIs per hemisphere | Diagnosis, age, sex, hemisphere-specific Euler number | 12 ROIs per hemisphere per contrast |
| Whole-cortex thickness | 34 ROIs per hemisphere | Same as targeted thickness | 34 ROIs per hemisphere per contrast |
| Whole-cortex surface area | 34 ROIs per hemisphere | Diagnosis, age, sex, ICV, hemisphere-specific Euler number | 34 ROIs per hemisphere per contrast |

The participant with an extreme right-hemisphere Euler number was excluded from specified **sensitivity analyses only**. QC treatments for implausible region-specific values must be applied using an authorised private QC record, not public participant identifiers.

## Data access and privacy
Participant-level MRI and clinical information is subject to institutional and ethics restrictions and is **not publicly available**. The master file, derived RDS files, individual measurements, participant identifiers, private QC mappings and participant-level diagnostics must not be committed to this repository. Locally authorised users should keep input and outputs outside the Git working tree. The `*.csv` ignore pattern is intentionally conservative; approved aggregate exports may be individually assessed for release.

## Reconstructed scripts
- `01_Data_Preparation_QC/01_prepare_final_dataset.R` — initial, privacy-conscious preparation from an authorised local RDS; this is **not yet the complete final QC pipeline** and should not be used to assert reproduction of published results.
- `02_Statistical_Analyses/01_primary_subcortical.R` — primary 9-ROI models and a reconstructed nine-ROI log-ventricle sensitivity family.
- `02_Statistical_Analyses/02_cortical_models.R` — targeted/whole-cortex thickness and surface-area models; optional private QC flag for right-sided sensitivity analyses.
- Subregional models, figures and numerical result comparisons will be added after source verification.

The local script requires environment variables `ABHGP_INPUT_RDS` and `ABHGP_PRIVATE_OUTPUT_DIR`; the first must point to an authorised final-study data frame, and the second must point to an output directory **outside** the repository.

## Verification issues to resolve
1. The original 2026 monolithic R script has syntax and object-order issues, including an incomplete `mutate()`, malformed directory creation, and incomplete amygdala analyses. It should be preserved **privately** as a historical original.
2. The main script's log-transformed lateral ventricle analysis does not implement the nine-ROI BH-FDR sensitivity family described in the appendix. Do not claim reproduction of q = 0.051 until its calculation is recovered.
3. The note to Appendix Table A4 mentions Euler number for hippocampal subfields, while Chapter 6 Methods and the main source model use age, sex and ICV. Resolve from archived final analytical output before changing models.
4. Figure scripts depend on objects created interactively and need explicit inputs to be standalone.
5. The current preparation script requires the authoritative FreeSurfer variable dictionary, private QC flags, and original final input definition before numerical equivalence can be tested.

## Data availability
The analyses were conducted on restricted clinical and neuroimaging data collected as part of the African Brain Health Genomics Project. Individual-level data and derived imaging measurements are not included due to ethical and institutional restrictions. The analytical code is shared for methodological transparency; reproduction of the results requires authorised access to the original data and documented quality-control decisions.
