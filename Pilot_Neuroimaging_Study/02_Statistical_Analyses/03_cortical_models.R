
# ============================================================
# PILOT NEUROIMAGING STUDY
# 03 - CORTICAL MORPHOLOGY ANALYSES
#
# Final analysis: August 2026
# Thesis: Chapter 6
# Appendix Tables: A3, A6, A7 and A8
#
# Targeted cortical thickness: 12 ROIs/hemisphere
# Whole-cortex thickness: 34 ROIs/hemisphere
# Cortical surface area: 34 ROIs/hemisphere
#
# Data are restricted and not publicly distributed.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(emmeans)
})

# ------------------------------------------------------------
# 1. LOAD RESTRICTED DATA
# ------------------------------------------------------------

private_dir <- Sys.getenv("ABHGP_PRIVATE_OUTPUT_DIR")

if (!nzchar(private_dir)) {
  stop("Set ABHGP_PRIVATE_OUTPUT_DIR.")
}

input_file <- file.path(
  private_dir,
  "final_prepared_private.rds"
)

if (!file.exists(input_file)) {
  stop("Final prepared dataset not found.")
}

dat <- readRDS(input_file)

# ------------------------------------------------------------
# 2. VALIDATE FINAL ANALYTICAL COHORT
# ------------------------------------------------------------

required <- c(
  "diagnosis4", "diagnosis3",
  "age", "sex", "ICV",
  "Euler_LH", "Euler_RH",
  "exclude_right_euler_sensitivity"
)

missing <- setdiff(required, names(dat))

if (length(missing) > 0) {
  stop(
    "Missing required variables: ",
    paste(missing, collapse = ", ")
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

dat$sex <- factor(dat$sex)

expected_counts <- c(
  Control = 61,
  GGE = 21,
  FE = 49,
  FDS = 19
)

stopifnot(
  nrow(dat) == 150L,
  identical(
    as.integer(table(dat$diagnosis4)),
    as.integer(expected_counts)
  ),
  !anyNA(dat$diagnosis4),
  !anyNA(dat$diagnosis3)
)

expected_diagnosis3 <- ifelse(
  dat$diagnosis4 == "Control",
  "Control",
  ifelse(
    dat$diagnosis4 == "FDS",
    "FDS",
    "Epilepsy"
  )
)

stopifnot(
  all(
    as.character(dat$diagnosis3) ==
      expected_diagnosis3
  )
)

numeric_covariates <- c(
  "age", "ICV", "Euler_LH", "Euler_RH"
)

stopifnot(
  all(vapply(
    dat[numeric_covariates],
    is.numeric,
    logical(1)
  ))
)

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
# 4. TARGETED CORTICAL REGIONS
# ------------------------------------------------------------

target_regions <- c(
  "caudalmiddlefrontal",
  "lateralorbitofrontal",
  "medialorbitofrontal",
  "parsopercularis",
  "parsorbitalis",
  "parstriangularis",
  "rostralmiddlefrontal",
  "superiorfrontal",
  "frontalpole",
  "precentral",
  "postcentral",
  "paracentral"
)

stopifnot(length(target_regions) == 12L)

# ------------------------------------------------------------
# 5. CORTICAL MODEL FUNCTION
# ------------------------------------------------------------

run_cortex <- function(
  data,
  suffix,
  group,
  target_only = FALSE,
  area = FALSE,
  hemispheres = c("L", "R")
) {

  stopifnot(
    group %in% names(contrast_spec),
    suffix %in% c("_thickavg", "_surfavg"),
    all(hemispheres %in% c("L", "R"))
  )

  if (target_only && area) {
    stop(
      "Targeted analysis was specified ",
      "for cortical thickness only."
    )
  }

  results <- list()

  for (hemi in hemispheres) {

    roi_cols <- grep(
      paste0("^", hemi, "_.*", suffix, "$"),
      names(data),
      value = TRUE
    )

    if (target_only) {
      roi_cols <- paste0(
        hemi, "_", target_regions, suffix
      )
    }

    expected_n <- if (target_only) 12L else 34L

    if (
      length(roi_cols) != expected_n ||
      !all(roi_cols %in% names(data))
    ) {
      stop(
        "Unexpected ROI set in hemisphere: ",
        hemi
      )
    }

    if (!all(vapply(
      data[roi_cols],
      is.numeric,
      logical(1)
    ))) {
      stop("Cortical ROI variables must be numeric.")
    }

    euler <- if (hemi == "L") {
      "Euler_LH"
    } else {
      "Euler_RH"
    }

    covariates <- if (area) {
      paste("age + sex + ICV +", euler)
    } else {
      paste("age + sex +", euler)
    }

    for (roi in roi_cols) {

      model <- lm(
        as.formula(
          paste(
            roi, "~", group, "+", covariates
          )
        ),
        data = data
      )

      em <- emmeans(
        model,
        specs = group
      )

      comparisons <- contrast(
        em,
        method = contrast_spec[[group]]
      )

      result <- as.data.frame(
        summary(
          comparisons,
          infer = c(TRUE, TRUE),
          adjust = "none"
        )
      )

      result$ROI <- roi
      result$Hemisphere <- if (
        hemi == "L"
      ) "Left" else "Right"

      results[[roi]] <- result
    }
  }

  bind_rows(results) %>%
    group_by(Hemisphere, contrast) %>%
    mutate(
      q_BH = p.adjust(
        p.value,
        method = "BH"
      )
    ) %>%
    ungroup() %>%
    select(
      Hemisphere,
      ROI,
      contrast,
      estimate,
      SE,
      df,
      lower.CL,
      upper.CL,
      t.ratio,
      p.value,
      q_BH
    ) %>%
    arrange(
      Hemisphere,
      contrast,
      p.value
    )
}

# ------------------------------------------------------------
# 6. PRIMARY TARGETED CORTICAL THICKNESS
# Appendix Table A3
# ------------------------------------------------------------

targeted_thickness_4group <- run_cortex(
  dat,
  "_thickavg",
  "diagnosis4",
  target_only = TRUE
)

targeted_thickness_3group <- run_cortex(
  dat,
  "_thickavg",
  "diagnosis3",
  target_only = TRUE
)

# ------------------------------------------------------------
# 7. EXPLORATORY WHOLE-CORTEX THICKNESS
# Appendix Table A6
# ------------------------------------------------------------

whole_thickness_4group <- run_cortex(
  dat,
  "_thickavg",
  "diagnosis4"
)

whole_thickness_3group <- run_cortex(
  dat,
  "_thickavg",
  "diagnosis3"
)

# ------------------------------------------------------------
# 8. EXPLORATORY SURFACE AREA
# Appendix Table A7
# ------------------------------------------------------------

surface_area_4group <- run_cortex(
  dat,
  "_surfavg",
  "diagnosis4",
  area = TRUE
)

surface_area_3group <- run_cortex(
  dat,
  "_surfavg",
  "diagnosis3",
  area = TRUE
)

# ------------------------------------------------------------
# 9. RIGHT-HEMISPHERE SENSITIVITY ANALYSES
# Appendix Table A8
# ------------------------------------------------------------

# The private preparation workflow must assign the
# original QC exclusion decision to this logical flag.
# No participant identifier appears in public code.

qc_flag <- dat$exclude_right_euler_sensitivity

if (!is.logical(qc_flag) || anyNA(qc_flag)) {
  stop("Private sensitivity QC flag must be logical and complete.")
}

if (sum(qc_flag) != 1L) {
  stop("Expected exactly one sensitivity exclusion.")
}

sensitivity_dat <- dat[!qc_flag, , drop = FALSE]

stopifnot(nrow(sensitivity_dat) == 149L)

targeted_thickness_4group_right_sensitivity <- run_cortex(
  sensitivity_dat, "_thickavg", "diagnosis4",
  target_only = TRUE, hemispheres = "R"
)

targeted_thickness_3group_right_sensitivity <- run_cortex(
  sensitivity_dat, "_thickavg", "diagnosis3",
  target_only = TRUE, hemispheres = "R"
)

whole_thickness_4group_right_sensitivity <- run_cortex(
  sensitivity_dat, "_thickavg", "diagnosis4",
  hemispheres = "R"
)

whole_thickness_3group_right_sensitivity <- run_cortex(
  sensitivity_dat, "_thickavg", "diagnosis3",
  hemispheres = "R"
)

surface_area_4group_right_sensitivity <- run_cortex(
  sensitivity_dat, "_surfavg", "diagnosis4",
  area = TRUE, hemispheres = "R"
)

surface_area_3group_right_sensitivity <- run_cortex(
  sensitivity_dat, "_surfavg", "diagnosis3",
  area = TRUE, hemispheres = "R"
)

# ------------------------------------------------------------
# 10. CHECK EXPECTED RESULT DIMENSIONS
# ------------------------------------------------------------

stopifnot(
  nrow(targeted_thickness_4group) == 72L,
  nrow(targeted_thickness_3group) == 72L,
  nrow(whole_thickness_4group) == 204L,
  nrow(whole_thickness_3group) == 204L,
  nrow(surface_area_4group) == 204L,
  nrow(surface_area_3group) == 204L,

  nrow(targeted_thickness_4group_right_sensitivity) == 36L,
  nrow(targeted_thickness_3group_right_sensitivity) == 36L,
  nrow(whole_thickness_4group_right_sensitivity) == 102L,
  nrow(whole_thickness_3group_right_sensitivity) == 102L,
  nrow(surface_area_4group_right_sensitivity) == 102L,
  nrow(surface_area_3group_right_sensitivity) == 102L
)

# ------------------------------------------------------------
# 11. REFERENCE CHECKS
# ------------------------------------------------------------

# Numerical targets are derived from the final thesis.
# These checks do not alter the model estimates.

check_reference <- function(
  data,
  roi,
  comparison,
  estimate,
  q,
  estimate_tol = 0.002,
  q_tol = 0.002
) {

  observed <- data %>%
    filter(
      ROI == roi,
      contrast == comparison
    )

  if (nrow(observed) != 1L) {
    stop("Expected exactly one matching result.")
  }

  data.frame(
    ROI = roi,
    contrast = comparison,
    estimate_matches = isTRUE(
      abs(observed$estimate - estimate) <=
        estimate_tol
    ),
    q_matches = isTRUE(
      abs(observed$q_BH - q) <= q_tol
    )
  )
}

reference_checks <- bind_rows(

  # Table A3 - FE, targeted cortical thickness

  check_reference(
    targeted_thickness_4group,
    "L_superiorfrontal_thickavg",
    "FE vs Control",
    -0.070, 0.045
  ),

  check_reference(
    targeted_thickness_4group,
    "R_superiorfrontal_thickavg",
    "FE vs Control",
    -0.076, 0.0048
  ),

  # Table A3 - combined epilepsy

  check_reference(
    targeted_thickness_3group,
    "R_superiorfrontal_thickavg",
    "Epilepsy vs Control",
    -0.063, 0.0132
  ),

  # Table A6 - whole-cortex thickness

  check_reference(
    whole_thickness_4group,
    "R_superiorfrontal_thickavg",
    "FE vs Control",
    -0.076, 0.0136
  ),

  # Table A7 - surface area

  check_reference(
    surface_area_3group,
    "R_parsorbitalis_surfavg",
    "Epilepsy vs Control",
    -62.843, 0.0272,
    estimate_tol = 0.1
  ),

  check_reference(
    surface_area_3group,
    "R_transversetemporal_surfavg",
    "Epilepsy vs Control",
    -26.439, 0.0272,
    estimate_tol = 0.1
  ),

  # Table A8 - sensitivity analysis

  check_reference(
    targeted_thickness_4group_right_sensitivity,
    "R_superiorfrontal_thickavg",
    "FE vs Control",
    -0.070, 0.0140
  ),

  check_reference(
    targeted_thickness_3group_right_sensitivity,
    "R_superiorfrontal_thickavg",
    "Epilepsy vs Control",
    -0.059, 0.0303
  ),

  check_reference(
    surface_area_3group_right_sensitivity,
    "R_parsorbitalis_surfavg",
    "Epilepsy vs Control",
    -64.122, 0.0461,
    estimate_tol = 0.1
  )
)

print(reference_checks)

if (!all(
  reference_checks$estimate_matches &
    reference_checks$q_matches
)) {
  warning(
    "One or more cortical results differ from ",
    "the final thesis appendix."
  )
}

# ------------------------------------------------------------
# 12. SAVE PRIVATE AGGREGATE OUTPUTS
# ------------------------------------------------------------

results <- list(

  targeted_thickness_4group =
    targeted_thickness_4group,

  targeted_thickness_3group =
    targeted_thickness_3group,

  whole_thickness_4group =
    whole_thickness_4group,

  whole_thickness_3group =
    whole_thickness_3group,

  surface_area_4group =
    surface_area_4group,

  surface_area_3group =
    surface_area_3group,

  targeted_thickness_4group_right_sensitivity =
    targeted_thickness_4group_right_sensitivity,

  targeted_thickness_3group_right_sensitivity =
    targeted_thickness_3group_right_sensitivity,

  whole_thickness_4group_right_sensitivity =
    whole_thickness_4group_right_sensitivity,

  whole_thickness_3group_right_sensitivity =
    whole_thickness_3group_right_sensitivity,

  surface_area_4group_right_sensitivity =
    surface_area_4group_right_sensitivity,

  surface_area_3group_right_sensitivity =
    surface_area_3group_right_sensitivity,

  reference_checks = reference_checks
)

saveRDS(
  results,
  file.path(
    private_dir,
    "cortical_aggregate_results_private.rds"
  )
)

message(
  "Cortical analyses completed. ",
  "Review numerical reference checks."
)

# ============================================================
# END
# ============================================================
