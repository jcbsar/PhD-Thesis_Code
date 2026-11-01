
#!/usr/bin/env Rscript

# Reproduce the retained two-factor FE confirmatory factor model
# reported in Chapter 3, Appendix Tables A14 and A15.
#
# Input: FE LDSC covariance object containing three FE subtypes
# and seven subcortical brain volumes.
#
# Brainstem is excluded from the retained six-region CFA model.
#
# This script does not rerun LDSC or modify covariance inputs.

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1L) {
  stop(
    "Usage: Rscript 02_fit_FE_retained_CFA.R ",
    "<FE_LDSC_object.rds>"
  )
}

if (!requireNamespace("GenomicSEM", quietly = TRUE)) {
  stop("The GenomicSEM package is required.")
}

input_file <- normalizePath(args[[1]], mustWork = TRUE)
covstruc <- readRDS(input_file)

# Confirm the expected phenotype set.
expected_traits <- c(
  "F_HS", "F_Neg", "F_Oth",
  "accumbens", "amygdala", "brainstem",
  "caudate", "pallidum", "putamen", "thalamus"
)

observed_traits <- colnames(covstruc$S)

if (!setequal(observed_traits, expected_traits) ||
    length(observed_traits) != length(expected_traits)) {
  stop(
    "Unexpected phenotype set in the FE covariance object: ",
    paste(observed_traits, collapse = ", ")
  )
}

cat("LDSC input:", input_file, "\n")
cat("Input phenotypes:", paste(observed_traits, collapse = ", "), "\n")
cat("Brainstem is excluded from the retained CFA model.\n")

# Retained two-factor CFA specification.
model_FE <- '
  FE =~ F_HS + F_Neg + F_Oth

  Brain =~ accumbens +
           amygdala +
           caudate +
           pallidum +
           putamen +
           thalamus

  FE ~~ Brain
'

# Fit the retained model.
fe_fit <- GenomicSEM::usermodel(
  covstruc = covstruc,
  model = model_FE,
  estimation = "DWLS",
  CFIcalc = TRUE,
  std.lv = TRUE,
  imp_cov = FALSE
)

# Display results.
cat("\nRetained FE CFA model fit:\n")
print(fe_fit$modelfit)

cat("\nParameter estimates:\n")
print(fe_fit$results)

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

output_dir <- file.path(section_dir, "outputs", "FE")

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Save the fitted model and results.
saveRDS(
  fe_fit,
  file.path(output_dir, "FE_brain_retained_fit.rds")
)

write.csv(
  fe_fit$modelfit,
  file.path(output_dir, "FE_brain_retained_model_fit.csv"),
  row.names = TRUE
)

write.csv(
  fe_fit$results,
  file.path(output_dir, "FE_brain_retained_parameter_estimates.csv"),
  row.names = FALSE
)

writeLines(
  model_FE,
  file.path(output_dir, "FE_brain_retained_model_syntax.txt")
)

# Record the exact LDSC input used for reproduction.
writeLines(
  c(
    paste("LDSC input:", input_file),
    paste("Input MD5:", unname(tools::md5sum(input_file))),
    paste("Phenotypes:", paste(observed_traits, collapse = ", "))
  ),
  file.path(output_dir, "FE_brain_retained_input_provenance.txt")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "FE_brain_retained_session_info.txt")
)

cat("\nResults saved to:", output_dir, "\n")
cat("FE retained CFA completed.\n")
