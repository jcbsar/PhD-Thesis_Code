
#!/usr/bin/env Rscript

# Reproduce the FE exploratory factor analyses reported
# in Chapter 3, Appendix Tables A11-A13.
#
# Input: FE_subtypes_brain_ldsc_July2026.rds
#
# The input contains three FE subtypes and seven
# subcortical brain-volume phenotypes.
#
# EFA uses a nearPD-smoothed LDSC genetic
# covariance matrix and promax rotation.
#
# This script does not rerun LDSC.

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1L) {
  stop(
    "Usage: Rscript 05_run_FE_EFA.R ",
    "<FE_LDSC_object.rds>"
  )
}

if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop("The Matrix package is required.")
}

input_file <- normalizePath(args[[1]], mustWork = TRUE)
covstruc <- readRDS(input_file)

# Validate the expected phenotype set.
expected_traits <- c(
  "F_HS", "F_Neg", "F_Oth",
  "accumbens", "amygdala", "brainstem",
  "caudate", "pallidum", "putamen", "thalamus"
)

observed_traits <- colnames(covstruc$S)

if (
  length(observed_traits) != length(expected_traits) ||
  !setequal(observed_traits, expected_traits)
) {
  stop("Unexpected phenotype set in FE LDSC object.")
}

cat("LDSC input:", input_file, "\n")
cat("Phenotypes:", paste(observed_traits, collapse = ", "), "\n")

# Smooth the genetic covariance matrix for EFA.
Ssmooth <- as.matrix(
  Matrix::nearPD(
    covstruc$S,
    corr = FALSE
  )$mat
)

dimnames(Ssmooth) <- list(
  observed_traits,
  observed_traits
)

smoothing_adjustment <- max(
  abs(Ssmooth - covstruc$S)
)

minimum_eigenvalue <- min(
  eigen(
    Ssmooth,
    symmetric = TRUE,
    only.values = TRUE
  )$values
)

cat(
  "Maximum smoothing adjustment:",
  smoothing_adjustment, "\n"
)

cat(
  "Minimum smoothed eigenvalue:",
  minimum_eigenvalue, "\n"
)

# Fit two-, three- and four-factor EFA solutions.
efa_models <- list()

for (k in 2:4) {

  fit <- stats::factanal(
    covmat = Ssmooth,
    factors = k,
    rotation = "promax"
  )

  efa_models[[paste0("EFA_", k)]] <- fit

  cat(
    "\n===== FE EFA:",
    k,
    "FACTOR(S) =====\n"
  )

  print(
    fit$loadings,
    cutoff = 0.20
  )

  cat("\nUniquenesses:\n")
  print(round(fit$uniquenesses, 3))
}

# Determine the output directory.
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
  "FE_EFA"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Save smoothed covariance matrix and model objects.
saveRDS(
  Ssmooth,
  file.path(output_dir, "FE_EFA_smoothed_S.rds")
)

saveRDS(
  efa_models,
  file.path(output_dir, "FE_EFA_models.rds")
)

# Export loadings, variance summaries and uniquenesses.
for (k in 2:4) {

  fit <- efa_models[[paste0("EFA_", k)]]

  # Remove the special loadings class for CSV export.
  L <- unclass(fit$loadings)

  loadings_df <- data.frame(
    Phenotype = rownames(L),
    as.data.frame.matrix(L),
    check.names = FALSE
  )

  write.csv(
    loadings_df,
    file.path(
      output_dir,
      paste0("FE_EFA_", k, "_loadings.csv")
    ),
    row.names = FALSE
  )

  ss <- colSums(L^2)
  proportion <- ss / nrow(L)

  variance_summary <- data.frame(
    Factor = colnames(L),
    SS_Loadings = ss,
    Proportion_Variance = proportion,
    Cumulative_Variance = cumsum(proportion)
  )

  write.csv(
    variance_summary,
    file.path(
      output_dir,
      paste0("FE_EFA_", k, "_variance.csv")
    ),
    row.names = FALSE
  )

  uniqueness_df <- data.frame(
    Phenotype = names(fit$uniquenesses),
    Uniqueness = as.numeric(fit$uniquenesses)
  )

  write.csv(
    uniqueness_df,
    file.path(
      output_dir,
      paste0("FE_EFA_", k, "_uniquenesses.csv")
    ),
    row.names = FALSE
  )
}

# Record input provenance and smoothing diagnostics.
writeLines(
  c(
    paste("LDSC input:", input_file),
    paste(
      "Input MD5:",
      unname(tools::md5sum(input_file))
    ),
    paste(
      "Phenotypes:",
      paste(observed_traits, collapse = ", ")
    ),
    paste(
      "Maximum smoothing adjustment:",
      smoothing_adjustment
    ),
    paste(
      "Minimum smoothed eigenvalue:",
      minimum_eigenvalue
    )
  ),
  file.path(output_dir, "FE_EFA_provenance.txt")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "FE_EFA_session_info.txt")
)

cat("\nResults saved to:", output_dir, "\n")
cat("FE EFA completed.\n")

