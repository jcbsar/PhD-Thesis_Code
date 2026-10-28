
#!/usr/bin/env Rscript

# Reproduce the retained two-factor GGE confirmatory factor model
# reported in Chapter 3, Appendix Tables A8 and A9.
#
# Input: Corrected LDSC covariance object using population
# prevalences of 0.00035 for CAE and 0.00015 for JME.
#
# The archived corrected LDSC object may contain brainstem,
# which is excluded from the retained six-region CFA model.
#
# This script does not rerun LDSC or modify covariance inputs.

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1L) {
  stop(
    "Usage: Rscript 01_fit_GGE_retained_CFA.R ",
    "<corrected_GGE_LDSC_object.rds>"
  )
}

if (!requireNamespace("GenomicSEM", quietly = TRUE)) {
  stop("The GenomicSEM package is required.")
}

input_file <- normalizePath(args[[1]], mustWork = TRUE)
covstruc <- readRDS(input_file)

# Confirm the phenotypes required by the retained model.
expected_traits <- c(
  "CAE", "GTCS", "JAE", "JME",
  "accumbens", "amygdala", "caudate",
  "pallidum", "putamen", "thalamus"
)

observed_traits <- colnames(covstruc$S)

valid_input <- (
  setequal(observed_traits, expected_traits) ||
  setequal(observed_traits, c(expected_traits, "brainstem"))
)

if (!valid_input) {
  stop(
    "Unexpected phenotype set in the GGE covariance object: ",
    paste(observed_traits, collapse = ", ")
  )
}

cat("LDSC input:", input_file, "\n")
cat("Input phenotypes:", paste(observed_traits, collapse = ", "), "\n")

if ("brainstem" %in% observed_traits) {
  cat("Brainstem is not included in the retained CFA model.\n")
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

# Record the exact LDSC input used for this reproduction.
writeLines(
  c(
    paste("LDSC input:", input_file),
    paste("Input MD5:", unname(tools::md5sum(input_file))),
    paste("Phenotypes:", paste(observed_traits, collapse = ", "))
  ),
  file.path(output_dir, "GGE_brain_retained_input_provenance.txt")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "GGE_brain_retained_session_info.txt")
)

cat("\nResults saved to:", output_dir, "\n")
cat("GGE retained CFA completed.\n")

