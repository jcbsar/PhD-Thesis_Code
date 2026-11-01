
# ============================================================
# PILOT NEUROIMAGING STUDY
# 01 - VERIFICATION OF FINAL IMAGING RESULTS

# This script checks:
#   1. Required analytical outputs
#   2. Numbers of comparisons
#   3. Selected estimates and BH-adjusted p-values
#   4. Numbers of FDR-significant subregional results
#
# Only aggregate statistical outputs are read.
# Participant-level data are not accessed.
#
# Numerical verification requires successful execution
# of the original analysis workflow with authorised data.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
})

# ------------------------------------------------------------
# 1. PRIVATE ANALYSIS DIRECTORY
# ------------------------------------------------------------

private_dir <- Sys.getenv("ABHGP_PRIVATE_OUTPUT_DIR")

if (!nzchar(private_dir)) {
  stop("Set ABHGP_PRIVATE_OUTPUT_DIR before verification.")
}

if (!dir.exists(private_dir)) {
  stop("Private analysis directory not found.")
}

primary_file <- file.path(
  private_dir,
  "primary_aggregate_results_private.rds"
)

cortical_file <- file.path(
  private_dir,
  "cortical_aggregate_results_private.rds"
)

subregional_dir <- file.path(
  private_dir,
  "subregional"
)

subregional_files <- c(
  thal_4 = "thalamic_nuclei_four_group_results.csv",
  thal_3 = "thalamic_nuclei_collapsed_group_results.csv",
  hippo_4 = "hippocampal_subfields_four_group_results.csv",
  hippo_3 = "hippocampal_subfields_collapsed_group_results.csv",
  amyg_4 = "amygdala_nuclei_four_group_results.csv",
  amyg_3 = "amygdala_nuclei_collapsed_group_results.csv"
)

subregional_paths <- file.path(
  subregional_dir,
  subregional_files
)

required_files <- c(
  primary_file,
  cortical_file,
  subregional_paths
)

if (!all(file.exists(required_files))) {
  stop(
    "Missing analytical outputs:\n",
    paste(
      required_files[!file.exists(required_files)],
      collapse = "\n"
    )
  )
}

# ------------------------------------------------------------
# 2. LOAD AGGREGATE RESULTS
# ------------------------------------------------------------

primary <- readRDS(primary_file)
cortical <- readRDS(cortical_file)

subregional <- lapply(
  subregional_paths,
  function(path) {
    read.csv(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
)

names(subregional) <- names(subregional_files)

# ------------------------------------------------------------
# 3. REQUIRED RESULTS AND STRUCTURE
# ------------------------------------------------------------

required_primary <- c(
  "primary_4group",
  "primary_3group",
  "primary_3group_log_vent_sensitivity"
)

required_cortical <- c(
  "targeted_thickness_4group",
  "targeted_thickness_3group",
  "whole_thickness_4group",
  "whole_thickness_3group",
  "surface_area_4group",
  "surface_area_3group",
  "targeted_thickness_4group_right_sensitivity",
  "targeted_thickness_3group_right_sensitivity",
  "whole_thickness_4group_right_sensitivity",
  "whole_thickness_3group_right_sensitivity",
  "surface_area_4group_right_sensitivity",
  "surface_area_3group_right_sensitivity"
)

if (!all(required_primary %in% names(primary))) {
  stop("Missing primary analysis result objects.")
}

if (!all(required_cortical %in% names(cortical))) {
  stop("Missing cortical analysis result objects.")
}

all_results <- c(
  primary[required_primary],
  subregional,
  cortical[required_cortical]
)

required_columns <- c(
  "ROI",
  "contrast",
  "estimate",
  "SE",
  "p.value",
  "q_BH"
)

for (nm in names(all_results)) {

  x <- all_results[[nm]]

  if (!is.data.frame(x) ||
      !all(required_columns %in% names(x))) {
    stop("Invalid result structure: ", nm)
  }

  if (anyNA(x[required_columns])) {
    stop("Missing result values in: ", nm)
  }

  numeric_columns <- c(
    "estimate", "SE", "p.value", "q_BH"
  )

  if (!all(vapply(
    x[numeric_columns],
    is.numeric,
    logical(1)
  ))) {
    stop("Non-numeric statistical values in: ", nm)
  }

  if (any(!is.finite(
    as.matrix(x[numeric_columns])
  ))) {
    stop("Non-finite statistical values in: ", nm)
  }

  if (any(
    x$p.value < 0 | x$p.value > 1 |
    x$q_BH < 0 | x$q_BH > 1
  )) {
    stop("Invalid p-values or q-values in: ", nm)
  }

  if (anyDuplicated(
    x[c("ROI", "contrast")]
  )) {
    stop("Duplicated ROI-contrast pairs in: ", nm)
  }
}

# ------------------------------------------------------------
# 4. VERIFY NUMBER OF COMPARISONS
# ------------------------------------------------------------

expected_rows <- c(
  primary_4group = 27,
  primary_3group = 27,
  primary_3group_log_vent_sensitivity = 27,

  thal_4 = 75,
  thal_3 = 75,
  hippo_4 = 54,
  hippo_3 = 54,
  amyg_4 = 27,
  amyg_3 = 27,

  targeted_thickness_4group = 72,
  targeted_thickness_3group = 72,
  whole_thickness_4group = 204,
  whole_thickness_3group = 204,
  surface_area_4group = 204,
  surface_area_3group = 204,

  targeted_thickness_4group_right_sensitivity = 36,
  targeted_thickness_3group_right_sensitivity = 36,
  whole_thickness_4group_right_sensitivity = 102,
  whole_thickness_3group_right_sensitivity = 102,
  surface_area_4group_right_sensitivity = 102,
  surface_area_3group_right_sensitivity = 102
)

dimension_checks <- data.frame(
  analysis = names(expected_rows),
  expected = as.integer(expected_rows),
  observed = vapply(
    all_results[names(expected_rows)],
    nrow,
    integer(1)
  ),
  row.names = NULL
) %>%
  mutate(
    pass = observed == expected
  )

# ------------------------------------------------------------
# 5. NUMERICAL REFERENCE CHECK FUNCTION
# ------------------------------------------------------------

# Reference values below are rounded as reported
# in the final thesis.
#
# Tolerances account for thesis reporting precision.
# They are NOT replacements for full-precision
# comparison with original analytical output.

check_reference <- function(
  results,
  analysis,
  roi,
  comparison,
  expected_estimate,
  expected_q,
  estimate_tolerance = 0.002,
  q_tolerance = 0.002
) {

  x <- results %>%
    filter(
      ROI == roi,
      contrast == comparison
    )

  if (nrow(x) != 1L) {
    stop(
      "Expected exactly one result for: ",
      analysis, " / ", roi, " / ", comparison
    )
  }

  estimate_match <- abs(
    x$estimate - expected_estimate
  ) <= estimate_tolerance

  q_match <- abs(
    x$q_BH - expected_q
  ) <= q_tolerance

  data.frame(
    analysis = analysis,
    ROI = roi,
    contrast = comparison,
    expected_estimate = expected_estimate,
    observed_estimate = x$estimate,
    expected_q = expected_q,
    observed_q = x$q_BH,
    estimate_match = estimate_match,
    q_match = q_match,
    pass = estimate_match & q_match
  )
}

# ------------------------------------------------------------
# 6. PRIMARY SUBCORTICAL VOLUMES
# Appendix Table A1
# ------------------------------------------------------------

primary_checks <- bind_rows(

  check_reference(
    primary$primary_3group,
    "Primary subcortical",
    "Thalamus",
    "Epilepsy vs Control",
    -706.2,
    0.0339,
    estimate_tolerance = 0.1
  ),

  check_reference(
    primary$primary_3group,
    "Primary subcortical",
    "LateralVentricles",
    "Epilepsy vs FDS",
    4244.6,
    0.0384,
    estimate_tolerance = 0.1
  ),

  check_reference(
    primary$primary_3group_log_vent_sensitivity,
    "Ventricular sensitivity",
    "LateralVentricles",
    "Epilepsy vs FDS",
    0.304,
    0.051
  )
)

# ------------------------------------------------------------
# 7. TARGETED CORTICAL THICKNESS
# Appendix Table A3
# ------------------------------------------------------------

targeted_checks <- bind_rows(

  check_reference(
    cortical$targeted_thickness_4group,
    "Targeted thickness",
    "L_superiorfrontal_thickavg",
    "FE vs Control",
    -0.070,
    0.045
  ),

  check_reference(
    cortical$targeted_thickness_4group,
    "Targeted thickness",
    "R_superiorfrontal_thickavg",
    "FE vs Control",
    -0.076,
    0.0048
  ),

  check_reference(
    cortical$targeted_thickness_3group,
    "Targeted thickness",
    "R_superiorfrontal_thickavg",
    "Epilepsy vs Control",
    -0.063,
    0.0132
  )
)

# ------------------------------------------------------------
# 8. WHOLE-CORTEX THICKNESS
# Appendix Table A6
# ------------------------------------------------------------

whole_checks <- check_reference(
  cortical$whole_thickness_4group,
  "Whole-cortex thickness",
  "R_superiorfrontal_thickavg",
  "FE vs Control",
  -0.076,
  0.0136
)

# ------------------------------------------------------------
# 9. CORTICAL SURFACE AREA
# Appendix Table A7
# ------------------------------------------------------------

area_checks <- bind_rows(

  check_reference(
    cortical$surface_area_3group,
    "Surface area",
    "R_parsorbitalis_surfavg",
    "Epilepsy vs Control",
    -62.8,
    0.0272,
    estimate_tolerance = 0.1
  ),

  check_reference(
    cortical$surface_area_3group,
    "Surface area",
    "R_transversetemporal_surfavg",
    "Epilepsy vs Control",
    -26.4,
    0.0272,
    estimate_tolerance = 0.1
  )
)

# ------------------------------------------------------------
# 10. RIGHT-HEMISPHERE SENSITIVITY
# Appendix Table A8
# ------------------------------------------------------------

sensitivity_checks <- bind_rows(

  check_reference(
    cortical$targeted_thickness_4group_right_sensitivity,
    "Right-hemisphere sensitivity",
    "R_superiorfrontal_thickavg",
    "FE vs Control",
    -0.070,
    0.0140
  ),

  check_reference(
    cortical$targeted_thickness_3group_right_sensitivity,
    "Right-hemisphere sensitivity",
    "R_superiorfrontal_thickavg",
    "Epilepsy vs Control",
    -0.059,
    0.0303
  ),

  check_reference(
    cortical$surface_area_3group_right_sensitivity,
    "Right-hemisphere sensitivity",
    "R_parsorbitalis_surfavg",
    "Epilepsy vs Control",
    -64.1,
    0.0461,
    estimate_tolerance = 0.1
  )
)

# ------------------------------------------------------------
# 11. SUBREGIONAL SIGNIFICANCE CHECKS
# Appendix Tables A2, A4 and A5
# ------------------------------------------------------------

count_significant <- function(
  data,
  comparison
) {
  sum(
    data$q_BH[
      data$contrast == comparison
    ] < 0.05
  )
}

subregional_checks <- data.frame(
  analysis = c(
    "Thalamus",
    "Thalamus",
    "Thalamus",
    "Thalamus",
    "Thalamus",
    "Hippocampus",
    "Hippocampus",
    "Amygdala",
    "Amygdala"
  ),

  comparison = c(
    "GGE vs Control",
    "FE vs Control",
    "FDS vs Control",
    "Epilepsy vs Control",
    "Epilepsy vs FDS",
    "All four-group contrasts",
    "All collapsed contrasts",
    "All four-group contrasts",
    "All collapsed contrasts"
  ),

  expected_significant = c(
    16, 1, 0, 11, 0,
    0, 0, 0, 0
  ),

  observed_significant = c(
    count_significant(
      subregional$thal_4,
      "GGE vs Control"
    ),
    count_significant(
      subregional$thal_4,
      "FE vs Control"
    ),
    count_significant(
      subregional$thal_4,
      "FDS vs Control"
    ),
    count_significant(
      subregional$thal_3,
      "Epilepsy vs Control"
    ),
    count_significant(
      subregional$thal_3,
      "Epilepsy vs FDS"
    ),
    sum(subregional$hippo_4$q_BH < 0.05),
    sum(subregional$hippo_3$q_BH < 0.05),
    sum(subregional$amyg_4$q_BH < 0.05),
    sum(subregional$amyg_3$q_BH < 0.05)
  )
) %>%
  mutate(
    pass = expected_significant ==
      observed_significant
  )

# ------------------------------------------------------------
# 12. COMBINE NUMERICAL CHECKS
# ------------------------------------------------------------

numerical_checks <- bind_rows(
  primary_checks,
  targeted_checks,
  whole_checks,
  area_checks,
  sensitivity_checks
)

# ------------------------------------------------------------
# 13. PRINT VERIFICATION RESULTS
# ------------------------------------------------------------

cat("\n====================================\n")
cat("FINAL IMAGING RESULTS VERIFICATION\n")
cat("====================================\n")

cat("\nCOMPARISON COUNT CHECKS\n")
print(dimension_checks)

cat("\nNUMERICAL REFERENCE CHECKS\n")
print(numerical_checks)

cat("\nSUBREGIONAL SIGNIFICANCE CHECKS\n")
print(subregional_checks)

# ------------------------------------------------------------
# 14. OVERALL VERIFICATION STATUS
# ------------------------------------------------------------

all_pass <- (
  all(dimension_checks$pass) &&
  all(numerical_checks$pass) &&
  all(subregional_checks$pass)
)

if (!all_pass) {
  stop(
    "One or more verification checks FAILED. ",
    "Investigate discrepancies before reporting ",
    "successful reproduction."
  )
}

cat("\nAll implemented verification checks PASSED.\n")
cat(
  "This does not constitute full numerical verification ",
  "of every appendix comparison.\n"
)

# ============================================================
# END
# ============================================================
