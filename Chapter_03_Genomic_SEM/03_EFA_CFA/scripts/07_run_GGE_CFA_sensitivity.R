
#!/usr/bin/env Rscript

# ============================================================
# Chapter 3: GGE CFA sensitivity analysis excluding CAE
# ============================================================
#
# Reproduces the GGE CFA sensitivity analysis excluding CAE.
#
# GGE indicators:
#   GTCS, JAE, JME
#
# Brain indicators:
#   accumbens, amygdala, caudate,
#   pallidum, putamen, thalamus
#
# CAE and brainstem are excluded from the fitted model.
#
# Original model specification:
#   GGE_no_CAE.R / model_GGE_no_CAE.R
#
# Corrected LDSC input:
#   GGE_ldsc_corrected_prevalence_Aug2026.rds
#
# Usage:
# Rscript 08_run_GGE_CFA_sensitivity.R <path_to_ldsc_rds>
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1L) {
  stop(
    "Usage: Rscript 08_run_GGE_CFA_sensitivity.R ",
    "<path_to_GGE_ldsc_corrected_prevalence_Aug2026.rds>"
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

gge_ldsc <- readRDS(input_file)

if (is.null(gge_ldsc$S) || is.null(gge_ldsc$V)) {
  stop("Input must contain LDSC S and V matrices.")
}

expected_traits <- c(
  "CAE",
  "GTCS",
  "JAE",
  "JME",
  "accumbens",
  "amygdala",
  "brainstem",
  "caudate",
  "pallidum",
  "putamen",
  "thalamus"
)

observed_traits <- colnames(gge_ldsc$S)

if (
  length(observed_traits) != length(expected_traits) ||
  !setequal(observed_traits, expected_traits)
) {
  stop("Unexpected phenotype set in GGE LDSC input.")
}

# Verify the corrected JME liability-scale heritability.
jme_index <- match("JME", observed_traits)
jme_h2 <- diag(gge_ldsc$S)[jme_index]

if (
  !is.finite(jme_h2) ||
  abs(jme_h2 - 0.2534618) > 0.001
) {
  stop(
    "JME heritability does not match the corrected ",
    "August LDSC input."
  )
}

cat("LDSC input:", input_file, "\n")
cat(
  "Phenotypes:",
  paste(observed_traits, collapse = ", "),
  "\n"
)
cat("JME corrected h2:", jme_h2, "\n")

# ------------------------------------------------------------
# 2. Define GGE sensitivity model excluding CAE
# ------------------------------------------------------------

model_no_CAE <- '
  GGE =~ GTCS + JAE + JME

  Brain =~ accumbens + amygdala + caudate +
           pallidum + putamen + thalamus

  GGE ~~ Brain
'

# No residual variance constraint is required.

# ------------------------------------------------------------
# 3. Fit the CFA model
# ------------------------------------------------------------

cat("\n===== GGE CFA EXCLUDING CAE =====\n")

fit_no_CAE <- GenomicSEM::usermodel(
  covstruc = gge_ldsc,
  model = model_no_CAE,
  estimation = "DWLS",
  CFIcalc = TRUE,
  std.lv = TRUE,
  imp_cov = FALSE
)

if (
  is.null(fit_no_CAE$modelfit) ||
  is.null(fit_no_CAE$results)
) {
  stop("GGE sensitivity CFA did not return results.")
}

cat("\nMODEL FIT:\n")
print(fit_no_CAE$modelfit)

# ------------------------------------------------------------
# 4. Extract latent-factor correlation
# ------------------------------------------------------------

correlation <- subset(
  fit_no_CAE$results,
  lhs == "GGE" &
    op == "~~" &
    rhs == "Brain"
)

if (nrow(correlation) != 1L) {
  stop("Could not identify the GGE-brain correlation.")
}

cat("\nGGE-BRAIN CORRELATION:\n")
print(correlation)

# ------------------------------------------------------------
# 5. Check model admissibility
# ------------------------------------------------------------

residuals <- subset(
  fit_no_CAE$results,
  op == "~~" &
    lhs == rhs &
    !(lhs %in% c("GGE", "Brain"))
)

if (nrow(residuals) != 9L) {
  stop("Unexpected number of residual variances.")
}

if (
  any(!is.finite(residuals$Unstand_Est)) ||
  any(residuals$Unstand_Est < 0)
) {
  stop("Negative or invalid residual variance detected.")
}

rg <- as.numeric(correlation$STD_All[1])

if (!is.finite(rg) || abs(rg) > 1) {
  stop("GGE-brain correlation is inadmissible.")
}

cat("\nModel admissibility checks: PASSED\n")

# ------------------------------------------------------------
# 6. Assemble sensitivity results
# ------------------------------------------------------------

mf <- fit_no_CAE$modelfit

summary_table <- data.frame(
  Model = "GGE_without_CAE",
  chisq = as.numeric(mf$chisq[1]),
  df = as.numeric(mf$df[1]),
  p_chisq = as.numeric(mf$p_chisq[1]),
  AIC = as.numeric(mf$AIC[1]),
  CFI = as.numeric(mf$CFI[1]),
  SRMR = as.numeric(mf$SRMR[1]),
  GGE_brain_rg = rg,
  GGE_brain_SE = as.numeric(
    correlation$STD_Genotype_SE[1]
  ),
  GGE_brain_p = as.numeric(
    correlation$p_value[1]
  )
)

cat("\n===== GGE SENSITIVITY RESULTS =====\n")
print(summary_table, row.names = FALSE)

# ------------------------------------------------------------
# 7. Validate independently reproduced statistics
# ------------------------------------------------------------

reference <- c(
  chisq = 78.21217,
  df = 26,
  AIC = 116.2122,
  CFI = 0.8831704,
  SRMR = 0.0835687,
  GGE_brain_rg = -0.181212,
  GGE_brain_SE = 0.0786081364,
  GGE_brain_p = 0.02115206
)

tolerances <- c(
  chisq = 0.02,
  df = 0,
  AIC = 0.02,
  CFI = 0.002,
  SRMR = 0.002,
  GGE_brain_rg = 0.002,
  GGE_brain_SE = 0.002,
  GGE_brain_p = 0.002
)

for (metric in names(reference)) {

  observed <- summary_table[[metric]][1]
  expected <- reference[[metric]]

  if (
    !is.finite(observed) ||
    abs(observed - expected) > tolerances[[metric]]
  ) {
    stop(
      "Reproduction check failed for: ",
      metric
    )
  }
}

cat("\nGGE sensitivity fit checks: PASSED\n")

# ------------------------------------------------------------
# 8. Prepare output directory
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
  "GGE_CFA_sensitivity"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 9. Save model results
# ------------------------------------------------------------

write.csv(
  summary_table,
  file.path(
    output_dir,
    "GGE_CFA_no_CAE_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  fit_no_CAE$results,
  file.path(
    output_dir,
    "GGE_CFA_no_CAE_parameters.csv"
  ),
  row.names = FALSE
)

write.csv(
  residuals,
  file.path(
    output_dir,
    "GGE_CFA_no_CAE_residuals.csv"
  ),
  row.names = FALSE
)

saveRDS(
  fit_no_CAE,
  file.path(
    output_dir,
    "GGE_CFA_no_CAE_model.rds"
  )
)

writeLines(
  model_no_CAE,
  file.path(
    output_dir,
    "GGE_CFA_no_CAE_model_syntax.txt"
  )
)

# ------------------------------------------------------------
# 10. Save provenance and session information
# ------------------------------------------------------------

provenance <- c(
  paste("LDSC input:", input_file),
  paste(
    "Input MD5:",
    unname(tools::md5sum(input_file))
  ),
  paste("Corrected JME h2:", jme_h2),
  "",
  "Chapter 3: GGE CFA sensitivity excluding CAE",
  "",
  "Original model specification recovered from:",
  "GGE_no_CAE.R",
  "model_GGE_no_CAE.R",
  "",
  "Corrected August LDSC covariance object used.",
  "",
  "GGE indicators:",
  "GTCS, JAE, JME",
  "",
  "Brain indicators:",
  "accumbens, amygdala, caudate,",
  "pallidum, putamen, thalamus",
  "",
  "CAE and brainstem excluded.",
  "No additional residual variance constraints.",
  "",
  "Estimation: DWLS",
  "Latent variances fixed to one.",
  "GenomicSEM usermodel(); CFIcalc = TRUE;",
  "imp_cov = FALSE.",
  "",
  "All estimated indicator residual variances",
  "were positive.",
  "",
  "The GGE-brain genetic correlation remained",
  "negative after excluding CAE.",
  "",
  "Model fit did not meet the prespecified",
  "CFI >= 0.90 or SRMR <= 0.08 criteria.",
  "",
  "The retained model and sensitivity model",
  "use different phenotype sets; their AIC",
  "values should not be directly compared.",
  "",
  "The three-factor brain-volume analysis",
  "is reproduced separately in",
  "04_compare_GGE_CFA_models.R."
)

writeLines(
  provenance,
  file.path(
    output_dir,
    "GGE_CFA_no_CAE_provenance.txt"
  )
)

writeLines(
  capture.output(sessionInfo()),
  file.path(
    output_dir,
    "GGE_CFA_no_CAE_session_info.txt"
  )
)

cat("\nResults saved to:", output_dir, "\n")
cat("GGE CFA sensitivity analysis completed.\n")

