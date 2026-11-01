
# ============================================================
# PILOT NEUROIMAGING STUDY
# 01 - PRIMARY SUBCORTICAL VOLUME ANALYSES
#
# Final analysis: August 2026
# Thesis: Chapter 6
# Reference: Appendix Table A1
#
# Participant-level data are restricted.
# This script requires authorised local access.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(emmeans)
})

# ------------------------------------------------------------
# 1. RESTRICTED DATA
# ------------------------------------------------------------

private_dir <- Sys.getenv("ABHGP_PRIVATE_OUTPUT_DIR")

if (!nzchar(private_dir)) {
  stop("Set ABHGP_PRIVATE_OUTPUT_DIR before running.")
}

input_file <- file.path(
  private_dir, "final_prepared_private.rds"
)

if (!file.exists(input_file)) {
  stop("Final prepared dataset not found.")
}

dat <- readRDS(input_file)

# ------------------------------------------------------------
# 2. FINAL COHORT AND VARIABLE VALIDATION
# ------------------------------------------------------------

primary_rois <- c(
  "Thalamus",
  "Caudate",
  "Putamen",
  "Pallidum",
  "Hippocampus",
  "Amygdala",
  "Accumbens",
  "LateralVentricles",
  "Brainstem"
)

required_variables <- c(
  "diagnosis4", "diagnosis3",
  "age", "sex", "ICV",
  primary_rois
)

missing_variables <- setdiff(
  required_variables, names(dat)
)

if (length(missing_variables) > 0) {
  stop(
    "Missing variables: ",
    paste(missing_variables, collapse = ", ")
  )
}

dat$diagnosis4 <- factor(
  as.character(dat$diagnosis4),
  levels = c("Control", "GGE", "FE", "FDS")
)

dat$diagnosis3 <- factor(
  as.character(dat$diagnosis3),
  levels = c("Control", "Epilepsy", "FDS")
)

expected_counts <- c(
  Control = 61,
  GGE = 21,
  FE = 49,
  FDS = 19
)

observed_counts <- table(dat$diagnosis4)

stopifnot(
  nrow(dat) == 150L,
  identical(
    as.integer(observed_counts),
    as.integer(expected_counts)
  ),
  !anyNA(dat$diagnosis4),
  !anyNA(dat$diagnosis3)
)

# Ensure collapsed diagnostic coding is consistent.
expected_diagnosis3 <- ifelse(
  dat$diagnosis4 == "Control", "Control",
  ifelse(
    dat$diagnosis4 == "FDS", "FDS", "Epilepsy"
  )
)

stopifnot(
  all(
    as.character(dat$diagnosis3) ==
      expected_diagnosis3
  )
)

# Check that numerical measurements are available.
numeric_variables <- c(
  "age", "ICV", primary_rois
)

stopifnot(
  all(vapply(
    dat[numeric_variables],
    is.numeric,
    logical(1)
  ))
)

dat$sex <- factor(dat$sex)

# ------------------------------------------------------------
# 3. DIAGNOSTIC CONTRASTS
# ------------------------------------------------------------

contrast_spec <- list(

  diagnosis4 = list(
    "GGE vs Control" = c(-1, 1, 0, 0),
    "FE vs Control" = c(-1, 0, 1, 0),
    "FDS vs Control" = c(-1, 0, 0, 1)
  ),

  diagnosis3 = list(
    "Epilepsy vs Control" = c(-1, 1, 0),
    "FDS vs Control" = c(-1, 0, 1),
    "Epilepsy vs FDS" = c(0, 1, -1)
  )
)

# ------------------------------------------------------------
# 4. REGION-WISE LINEAR MODELS
# ------------------------------------------------------------

fit_primary <- function(
  data,
  group_var,
  log_ventricles = FALSE
) {

  stopifnot(
    group_var %in% c("diagnosis4", "diagnosis3")
  )

  results <- lapply(primary_rois, function(roi) {

    response <- roi

    if (
      log_ventricles &&
      roi == "LateralVentricles"
    ) {

      valid_values <- data[[roi]][
        !is.na(data[[roi]])
      ]

      if (
        length(valid_values) == 0 ||
        any(!is.finite(valid_values)) ||
        any(valid_values <= 0)
      ) {
        stop(
          "Ventricular volumes must be positive ",
          "and finite for log transformation."
        )
      }

      response <- "log(LateralVentricles)"
    }

    model_formula <- as.formula(
      paste(
        response,
        "~", group_var,
        "+ age + sex + ICV"
      )
    )

    model <- lm(
      model_formula,
      data = data
    )

    em <- emmeans(
      model,
      specs = group_var
    )

    comparisons <- contrast(
      em,
      method = contrast_spec[[group_var]]
    )

    result <- as.data.frame(
      summary(
        comparisons,
        infer = c(TRUE, TRUE),
        adjust = "none"
      )
    )

    result$ROI <- roi

    result$scale <- if (
      log_ventricles &&
      roi == "LateralVentricles"
    ) {
      "log"
    } else {
      "original"
    }

    result
  })

  bind_rows(results) %>%
    group_by(contrast) %>%
    mutate(
      q_BH = p.adjust(
        p.value,
        method = "BH"
      )
    ) %>%
    ungroup() %>%
    select(
      ROI,
      contrast,
      scale,
      estimate,
      SE,
      df,
      lower.CL,
      upper.CL,
      t.ratio,
      p.value,
      q_BH
    ) %>%
    arrange(contrast, p.value)
}

# ------------------------------------------------------------
# 5. FOUR-GROUP ANALYSIS
# ------------------------------------------------------------

primary_4group <- fit_primary(
  dat,
  group_var = "diagnosis4"
)

# ------------------------------------------------------------
# 6. COLLAPSED THREE-GROUP ANALYSIS
# ------------------------------------------------------------

primary_3group <- fit_primary(
  dat,
  group_var = "diagnosis3"
)

# ------------------------------------------------------------
# 7. LATERAL VENTRICULAR SENSITIVITY
# ------------------------------------------------------------

# The nine-region family is retained.
# Only the ventricular outcome is log-transformed.
#
# This reconstruction follows the thesis-described
# sensitivity correction family, but numerical
# equivalence has not yet been established.

primary_3group_log_vent_sensitivity <- fit_primary(
  dat,
  group_var = "diagnosis3",
  log_ventricles = TRUE
)

# Extract the ventricular sensitivity contrasts.
ventricular_sensitivity <- primary_3group_log_vent_sensitivity %>%
  filter(
    ROI == "LateralVentricles"
  )

# ------------------------------------------------------------
# 8. MODEL OUTPUT CHECKS
# ------------------------------------------------------------

stopifnot(
  nrow(primary_4group) == 27L,
  nrow(primary_3group) == 27L,
  nrow(primary_3group_log_vent_sensitivity) == 27L,
  nrow(ventricular_sensitivity) == 3L
)

# ------------------------------------------------------------
# 9. REFERENCE CHECKS
# ------------------------------------------------------------

# These are comparison targets from Appendix A1.
# They are not substituted into model outputs.

check_reference <- function(
  results, roi, comparison,
  expected_estimate,
  expected_q,
  estimate_tol = 0.1,
  q_tol = 0.001
) {

  observed <- results %>%
    filter(
      ROI == roi,
      contrast == comparison
    )

  if (nrow(observed) != 1L) {
    stop("Expected exactly one reference result.")
  }

  data.frame(
    ROI = roi,
    contrast = comparison,
    estimate_matches = isTRUE(
      abs(observed$estimate - expected_estimate) <=
        estimate_tol
    ),
    q_matches = isTRUE(
      abs(observed$q_BH - expected_q) <=
        q_tol
    )
  )
}

reference_checks <- bind_rows(

  check_reference(
    primary_3group,
    "Thalamus",
    "Epilepsy vs Control",
    -706.2,
    0.0339
  ),

  check_reference(
    primary_3group,
    "LateralVentricles",
    "Epilepsy vs FDS",
    4244.6,
    0.0384
  ),

  check_reference(
    primary_3group_log_vent_sensitivity,
    "LateralVentricles",
    "Epilepsy vs FDS",
    0.304,
    0.051,
    estimate_tol = 0.002,
    q_tol = 0.002
  )
)

print(reference_checks)

if (!all(
  reference_checks$estimate_matches &
  reference_checks$q_matches
)) {
  warning(
    "One or more results differ from Appendix A1. ",
    "Investigate before claiming numerical verification."
  )
}

# ------------------------------------------------------------
# 10. SAVE RESULTS PRIVATELY
# ------------------------------------------------------------

# Contains aggregate model outputs only.
# Files are retained outside the public repository.

saveRDS(
  list(
    primary_4group = primary_4group,
    primary_3group = primary_3group,
    primary_3group_log_vent_sensitivity =
      primary_3group_log_vent_sensitivity,
    ventricular_sensitivity =
      ventricular_sensitivity,
    reference_checks = reference_checks
  ),
  file.path(
    private_dir,
    "primary_aggregate_results_private.rds"
  )
)

message(
  "Primary subcortical analyses completed. ",
  "Review reference checks before validation."
)

# ============================================================
# END
# ============================================================
