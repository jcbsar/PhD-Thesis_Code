
#!/usr/bin/env Rscript

# Reproduce the GGE CFA model comparisons reported
# in Chapter 3, Appendix Table A8.
#
# Input: GGE_ldsc_corrected_prevalence_Aug2026.rds
#
# Models:
# 1. One-factor epilepsy and brain-volume model
# 2. Retained two-factor GGE and brain-volume model
# 3. Three-factor GGE, striatal/thalamic and limbic model
#
# Brainstem is excluded from all three specifications.
# The CAE residual variance is fixed to zero.
#
# This script does not rerun LDSC or modify its inputs.

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1L) {
  stop(
    "Usage: Rscript 04_compare_GGE_CFA_models.R ",
    "<corrected_GGE_LDSC_object.rds>"
  )
}

if (!requireNamespace("GenomicSEM", quietly = TRUE)) {
  stop("The GenomicSEM package is required.")
}

input_file <- normalizePath(args[[1]], mustWork = TRUE)
covstruc <- readRDS(input_file)

# Validate the August LDSC input.
expected_traits <- c(
  "CAE", "GTCS", "JAE", "JME",
  "accumbens", "amygdala", "brainstem",
  "caudate", "pallidum", "putamen", "thalamus"
)

observed_traits <- colnames(covstruc$S)

if (length(observed_traits) != length(expected_traits) ||
    !setequal(observed_traits, expected_traits)) {
  stop("Unexpected phenotype set in the GGE LDSC object.")
}

jme_index <- match("JME", observed_traits)
jme_h2 <- covstruc$S[jme_index, jme_index]

if (!is.finite(jme_h2) ||
    abs(jme_h2 - 0.2534618) > 0.001) {
  stop("LDSC input does not match the corrected August JME estimate.")
}

cat("LDSC input:", input_file, "\n")
cat("JME heritability:", jme_h2, "\n")

# Model specifications verified against Appendix Table A8.
models <- list(

  One_factor = '
    Common =~ CAE + GTCS + JAE + JME +
              accumbens + amygdala + caudate +
              pallidum + putamen + thalamus

    CAE ~~ 0*CAE
  ',

  Two_factor = '
    GGE =~ CAE + GTCS + JAE + JME

    Brain =~ accumbens + amygdala + caudate +
             pallidum + putamen + thalamus

    CAE ~~ 0*CAE

    GGE ~~ Brain
  ',

  Three_factor = '
    GGE =~ CAE + GTCS + JAE + JME

    StrThal =~ caudate + pallidum + putamen + thalamus

    Limbic =~ amygdala + accumbens

    CAE ~~ 0*CAE
  '
)

# Fit all three CFA specifications.
fits <- list()
comparison <- list()

for (model_name in names(models)) {

  cat("\nFitting:", model_name, "\n")

  fit <- GenomicSEM::usermodel(
    covstruc = covstruc,
    model = models[[model_name]],
    estimation = "DWLS",
    CFIcalc = TRUE,
    std.lv = TRUE,
    imp_cov = FALSE
  )

  fits[[model_name]] <- fit

  stats <- fit$modelfit

  comparison[[model_name]] <- data.frame(
    Model = model_name,
    Chi_square = as.numeric(stats$chisq[1]),
    df = as.numeric(stats$df[1]),
    p_value = as.numeric(stats$p_chisq[1]),
    AIC = as.numeric(stats$AIC[1]),
    CFI = as.numeric(stats$CFI[1]),
    SRMR = as.numeric(stats$SRMR[1])
  )

  print(stats)
}

comparison_table <- do.call(
  rbind,
  comparison
)

rownames(comparison_table) <- NULL

cat("\nGGE CFA MODEL COMPARISON\n")
print(comparison_table)

# Identify the script directory.
script_arg <- grep(
  "^--file=",
  commandArgs(trailingOnly = FALSE),
  value = TRUE
)

if (length(script_arg) != 1L) {
  stop("Run this script using Rscript.")
}

script_file <- sub("^--file=", "", script_arg[[1]])

section_dir <- normalizePath(
  file.path(dirname(script_file), ".."),
  mustWork = TRUE
)

output_dir <- file.path(
  section_dir,
  "outputs",
  "GGE_CFA_comparison"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Save comparison statistics.
write.csv(
  comparison_table,
  file.path(output_dir, "GGE_CFA_model_comparison.csv"),
  row.names = FALSE
)

# Save model objects and specifications.
saveRDS(
  fits,
  file.path(output_dir, "GGE_CFA_comparison_fits.rds")
)

for (model_name in names(models)) {
  writeLines(
    models[[model_name]],
    file.path(
      output_dir,
      paste0(model_name, "_syntax.txt")
    )
  )
}

# Record input provenance.
writeLines(
  c(
    paste("LDSC input:", input_file),
    paste("Input MD5:", unname(tools::md5sum(input_file))),
    paste("JME heritability:", jme_h2),
    paste("Phenotypes:", paste(observed_traits, collapse = ", "))
  ),
  file.path(output_dir, "GGE_CFA_comparison_provenance.txt")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "GGE_CFA_comparison_session_info.txt")
)

cat("\nResults saved to:", output_dir, "\n")
cat("GGE CFA model comparison completed.\n")
