
#!/usr/bin/env Rscript

# Reproduce the retained two-factor GGE confirmatory factor model
# reported in Chapter 3, Appendix Table A9.
#
# Input: July-derived six-region LDSC covariance object.
# This script does not rerun LDSC or modify covariance inputs.

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1L) {
  stop(
    "Usage: Rscript 01_fit_GGE_retained_CFA.R ",
    "ldsc_GGE_brain_6region.rds"
  )
}

if (!requireNamespace("GenomicSEM", quietly = TRUE)) {
  stop("The GenomicSEM package is required.")
}

input_file <- normalizePath(args[[1]], mustWork = TRUE)
covstruc <- readRDS(input_file)

# Confirm the expected phenotype set.
expected_traits <- c(
  "CAE", "GTCS", "JAE", "JME",
  "accumbens", "amygdala", "caudate",
  "pallidum", "putamen", "thalamus"
)

if (!setequal(colnames(covstruc$S), expected_traits)) {
  stop("Unexpected phenotypes in the GGE covariance object.")
}

# Retained two-factor CFA specification.
model_GGE <- '
  GGE =~ CAE + GTCS + JAE + JME

  Brain =~ accumbens +
           amygdala +
           caudate +
           pallidum +
           putamen +
           thalamus

  CAE ~~ 0*CAE

  GGE ~~ Brain
'

# Fit the retained model.
gge_fit <- GenomicSEM::usermodel(
  covstruc = covstruc,
  model = model_GGE,
  estimation = "DWLS",
  CFIcalc = TRUE,
  std.lv = TRUE,
  imp_cov = FALSE
)

# Display results.
cat("\nRetained GGE CFA model fit:\n")
print(gge_fit$modelfit)

cat("\nParameter estimates:\n")
print(gge_fit$results)

# Identify the directory containing this script.
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

output_dir <- file.path(section_dir, "outputs", "GGE")
dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Save the fitted model and results.
saveRDS(
  gge_fit,
  file.path(output_dir, "GGE_brain_retained_fit.rds")
)

write.csv(
  gge_fit$modelfit,
  file.path(output_dir, "GGE_brain_retained_model_fit.csv"),
  row.names = TRUE
)

write.csv(
  gge_fit$results,
  file.path(output_dir, "GGE_brain_retained_parameter_estimates.csv"),
  row.names = FALSE
)

writeLines(
  model_GGE,
  file.path(output_dir, "GGE_brain_retained_model_syntax.txt")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "GGE_brain_retained_session_info.txt")
)

cat("\nResults saved to:", output_dir, "\n")
cat("GGE retained CFA completed.\n")
