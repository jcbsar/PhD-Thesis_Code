
# Verification of Neuroimaging Results

This directory contains the verification script for the final 2026 pilot neuroimaging analyses.

## Verification script

`01_verify_imaging_results.R`

The script compares aggregate statistical outputs with selected findings reported in Appendix Tables A1-A8.

Verification includes:

- Expected numbers of regions and diagnostic comparisons.
- Selected effect estimates and BH-FDR-adjusted p-values.
- Numbers of FDR-significant thalamic nuclei.
- Absence of FDR-significant hippocampal and amygdala subregional findings.
- Selected cortical thickness, surface-area and sensitivity findings.

## Verification status

**Implemented but not yet executed against the final restricted dataset.**

The current checks assess selected reported findings and output structure. They do not establish complete numerical equivalence across every appendix comparison.

Verification requires authorised access to the final analytical dataset and associated private QC documentation. Participant-level data and identifying QC information are not publicly available.
