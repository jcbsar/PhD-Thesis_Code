
#!/usr/bin/env Rscript

# ============================================================
# Chapter 3: FE CFA sensitivity analyses
# ============================================================
#
# Reproduces the final FE sensitivity analyses:
#
# 1. Excluding F_Oth
# 2. Excluding F_HS
#
# Both models retain the six-region brain-volume factor.
# The two remaining FE loadings are constrained to equality
# for identification.
#
# Results correspond to the final archived verification
# models and the Chapter 3 thesis results.
#
# Input:
# FE_subtypes_brain_ldsc_July2026.rds
#
# Usage:
# Rscript 07_run_FE_CFA_sensitivity.R <path_to_ldsc_rds>
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1L) {
  stop(
    "Usage: Rscript 07_run_FE_CFA_sensitivity.R ",
    "<path_to_FE_subtypes_brain_ldsc_July2026.rds>"
  )
}

if (!requireNamespace("GenomicSEM", quietly = TRUE)) {
  stop("GenomicSEM is required.")
}

# ------------------------------------------------------------
# 1. Load and validate LDSC input
# ------------------------------------------------------------

input_file <- normalizePath(
  args[[1]],
  mustWork = TRUE
)

fe_ldsc <- readRDS(input_file)

if (is.null(fe_ldsc$S) || is.null(fe_ldsc$V)) {
  stop("Input must contain LDSC S and V matrices.")
}

expected_traits <- c(
  "F_HS",
  "F_Neg",
  "F_Oth",
  "accumbens",
  "amygdala",
  "brainstem",
  "caudate",
  "pallidum",
  "putamen",
  "thalamus"
)

observed_traits <- colnames(fe_ldsc$S)

if (
  length(observed_traits) != length(expected_traits) ||
  !setequal(observed_traits, expected_traits)
) {
  stop("Unexpected phenotype set in FE LDSC input.")
}

cat("LDSC input:", input_file, "\n")
cat(
  "Phenotypes:",
  paste(observed_traits, collapse = ", "),
  "\n"
)

# ------------------------------------------------------------
# 2. Define final sensitivity models
# ------------------------------------------------------------

# Excluding F_Oth:
# Remaining FE loadings constrained to equality.

model_no_FOth <- '
  FE =~ a*F_HS + a*F_Neg

  Brain =~ accumbens + amygdala + caudate +
           pallidum + putamen + thalamus

  FE ~~ Brain
'

# Excluding F_HS:
# Remaining FE loadings constrained to equality.

model_no_FHS <- '
  FE =~ a*F_Neg + a*F_Oth

  Brain =~ accumbens + amygdala + caudate +
           pallidum + putamen + thalamus

  FE ~~ Brain
'

# ------------------------------------------------------------
# 3. Fit both models
# ------------------------------------------------------------

fit_cfa <- function(model) {
  GenomicSEM::usermodel(
    covstruc = fe_ldsc,
    model = model,
    estimation = "DWLS",
    CFIcalc = TRUE,
    std.lv = TRUE,
    imp_cov = FALSE
  )
}

cat("\n===== EXCLUDING F_Oth =====\n")

fit_no_FOth <- fit_cfa(model_no_FOth)

cat("\n===== EXCLUDING F_HS =====\n")

fit_no_FHS <- fit_cfa(model_no_FHS)

fits <- list(
  Without_F_Oth = fit_no_FOth,
  Without_F_HS = fit_no_FHS
)

if (
  any(vapply(
    fits,
    function(x) is.null(x$modelfit) ||
      is.null(x$results),
    logical(1)
  ))
) {
  stop("At least one sensitivity model failed.")
}

# ------------------------------------------------------------
# 4. Extract model statistics and FE-brain correlations
# ------------------------------------------------------------

extract_summary <- function(fit, model_name) {

  mf <- fit$modelfit

  correlation <- subset(
    fit$results,
    lhs == "FE" &
      op == "~~" &
      rhs == "Brain"
  )

  if (nrow(correlation) != 1L) {
    stop(
      "Could not identify FE-brain correlation for ",
      model_name
    )
  }

  loadings <- subset(
    fit$results,
    lhs == "FE" &
      op == "=~"
  )

  if (nrow(loadings) != 2L) {
    stop(
      "Expected exactly two FE loadings for ",
      model_name
    )
  }

  # Validate equality of unstandardised loadings.

  loading_values <- as.numeric(
    loadings$Unstand_Est
  )

  if (
    abs(loading_values[1] - loading_values[2]) > 1e-6
  ) {
    stop(
      "FE loading equality constraint failed for ",
      model_name
    )
  }

  data.frame(
    Model = model_name,
    chisq = as.numeric(mf$chisq[1]),
    df = as.numeric(mf$df[1]),
    p_chisq = as.numeric(mf$p_chisq[1]),
    AIC = as.numeric(mf$AIC[1]),
    CFI = as.numeric(mf$CFI[1]),
    SRMR = as.numeric(mf$SRMR[1]),
    FE_brain_rg = as.numeric(correlation$STD_All[1]),
    FE_brain_SE = as.numeric(
      correlation$STD_Genotype_SE[1]
    ),
    FE_brain_p = as.numeric(correlation$p_value[1]),
    FE_loading_1 = loading_values[1],
    FE_loading_2 = loading_values[2]
  )
}

comparison <- do.call(
  rbind,
  lapply(names(fits), function(name) {
    extract_summary(fits[[name]], name)
  })
)

rownames(comparison) <- NULL

cat("\n===== FE SENSITIVITY RESULTS =====\n")
print(comparison, row.names = FALSE)

# ------------------------------------------------------------
# 5. Validate against final archived verification results
# ------------------------------------------------------------

reference <- data.frame(
  Model = c("Without_F_Oth", "Without_F_HS"),
  chisq = c(64.33789, 67.23327),
  df = c(20, 20),
  AIC = c(96.33789, 99.23327),
  CFI = c(0.8423788, 0.8456436),
  SRMR = c(0.07223118, 0.1071396),
  FE_brain_rg = c(-0.009587434, -0.06318293),
  FE_brain_SE = c(0.1375355884, 0.1920549818),
  FE_brain_p = c(0.9616861, 0.4553026)
)

tolerances <- c(
  chisq = 0.02,
  df = 0,
  AIC = 0.02,
  CFI = 0.002,
  SRMR = 0.002,
  FE_brain_rg = 0.002,
  FE_brain_SE = 0.002,
  FE_brain_p = 0.002
)

if (!identical(comparison$Model, reference$Model)) {
  stop("Sensitivity model order does not match reference.")
}

for (metric in names(tolerances)) {

  observed <- comparison[[metric]]
  expected <- reference[[metric]]

  if (
    any(!is.finite(observed)) ||
    any(
      abs(observed - expected) >
        tolerances[[metric]]
    )
  ) {
    stop(
      "Reproduction check failed for: ",
      metric
    )
  }
}

cat("\nFE sensitivity fit checks: PASSED\n")
cat("FE loading equality checks: PASSED\n")

# ------------------------------------------------------------
# 6. Prepare output directory
# ------------------------------------------------------------

script_arg <- grep(
  "^--file=",
  commandArgs(trailingOnly = FALSE),
  value = TRUE
)

if (length(script_arg) != 1L) {
  stop("Unable to determine script location.")
}

script_file <- sub(
  "^--file=",
  "",
  script_arg[[1]]
)

section_dir <- normalizePath(
  file.path(dirname(script_file), ".."),
  mustWork = TRUE
)

output_dir <- file.path(
  section_dir,
  "outputs",
  "FE_CFA_sensitivity"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 7. Save results and fitted models
# ------------------------------------------------------------

write.csv(
  comparison,
  file.path(
    output_dir,
    "FE_CFA_sensitivity_summary.csv"
  ),
  row.names = FALSE
)

saveRDS(
  fits,
  file.path(
    output_dir,
    "FE_CFA_sensitivity_models.rds"
  )
)

for (name in names(fits)) {

  write.csv(
    fits[[name]]$results,
    file.path(
      output_dir,
      paste0(name, "_parameters.csv")
    ),
    row.names = FALSE
  )
}

writeLines(
  c(
    "MODEL EXCLUDING F_Oth",
    model_no_FOth,
    "",
    "MODEL EXCLUDING F_HS",
    model_no_FHS
  ),
  file.path(
    output_dir,
    "FE_CFA_sensitivity_model_syntax.txt"
  )
)

# ------------------------------------------------------------
# 8. Save provenance and session information
# ------------------------------------------------------------

provenance <- c(
  paste("LDSC input:", input_file),
  paste(
    "Input MD5:",
    unname(tools::md5sum(input_file))
  ),
  "",
  "Chapter 3: FE CFA sensitivity analyses",
  "",
  "Final verification reference:",
  "FE_CFA_final_verification/",
  "Without_F_Oth_model.rds",
  "Without_F_HS_model.rds",
  "",
  "Both sensitivity models exclude brainstem.",
  "The six-region brain-volume factor is retained.",
  "",
  "The two remaining FE indicators have",
  "equal unstandardised factor loadings.",
  "",
  "Estimation: DWLS",
  "Latent variances fixed to one.",
  "GenomicSEM usermodel(); CFIcalc = TRUE;",
  "imp_cov = FALSE.",
  "",
  "Interpretation:",
  "Neither sensitivity analysis identified a",
  "statistically significant FE-brain correlation.",
  "Some loading standard errors were very large,",
  "indicating parameter instability.",
  "",
  "These scripts reproduce the final",
  "equality-constrained sensitivity analyses,",
  "not the earlier unconstrained August models."
)

writeLines(
  provenance,
  file.path(
    output_dir,
    "FE_CFA_sensitivity_provenance.txt"
  )
)

writeLines(
  capture.output(sessionInfo()),
  file.path(
    output_dir,
    "FE_CFA_sensitivity_session_info.txt"
  )
)

cat("\nResults saved to:", output_dir, "\n")
cat("FE CFA sensitivity analyses completed.\n")
