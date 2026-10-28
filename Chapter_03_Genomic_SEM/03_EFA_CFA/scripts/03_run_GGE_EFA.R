
#!/usr/bin/env Rscript

# Reproduce the GGE exploratory factor analyses reported
# in Chapter 3, Appendix Tables A5–A7.
#
# Input: GGE_ldsc_corrected_prevalence_Aug2026.rds
#
# The input contains four GGE subtypes and seven
# subcortical brain-volume phenotypes.
#
# EFA is performed using a nearPD-smoothed LDSC
# genetic covariance matrix.

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1L) {
  stop(
    "Usage: Rscript 03_run_GGE_EFA.R ",
    "<corrected_GGE_LDSC_object.rds>"
  )
}

if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop("The Matrix package is required.")
}

input_file <- normalizePath(args[[1]], mustWork = TRUE)
covstruc <- readRDS(input_file)

expected_traits <- c(
  "CAE", "GTCS", "JAE", "JME",
  "accumbens", "amygdala", "brainstem",
  "caudate", "pallidum", "putamen", "thalamus"
)

observed_traits <- colnames(covstruc$S)

if (length(observed_traits) != length(expected_traits) ||
    !setequal(observed_traits, expected_traits)) {
  stop("Unexpected phenotype set in GGE LDSC object.")
}

# Check the corrected August LDSC specification.
jme_index <- match("JME", observed_traits)
jme_h2 <- covstruc$S[jme_index, jme_index]

if (!is.finite(jme_h2) ||
    abs(jme_h2 - 0.2534618) > 0.001) {
  stop(
    "Input does not match the verified corrected ",
    "August JME heritability."
  )
}

cat("LDSC input:", input_file, "\n")
cat("JME heritability:", jme_h2, "\n")

# Smooth the genetic covariance matrix for EFA.
Ssmooth <- as.matrix(
  Matrix::nearPD(
    covstruc$S,
    corr = FALSE
  )$mat
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

# Fit the one- to four-factor EFA solutions.
efa_models <- list()

for (k in 1:4) {

  rotation_method <- if (k == 1) {
    "none"
  } else {
    "promax"
  }

  fit <- stats::factanal(
    covmat = Ssmooth,
    factors = k,
    rotation = rotation_method
  )

  efa_models[[paste0("EFA_", k)]] <- fit

  cat(
    "\n===== GGE EFA:",
    k,
    "FACTOR(S) =====\n"
  )

  print(
    fit$loadings,
    cutoff = if (k == 1) 0 else 0.20,
    sort = TRUE
  )
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

script_file <- sub(
  "^--file=", "",
  script_arg[[1]]
)

section_dir <- normalizePath(
  file.path(dirname(script_file), ".."),
  mustWork = TRUE
)

output_dir <- file.path(
  section_dir,
  "outputs",
  "GGE_EFA"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Save smoothed covariance matrix and models.
saveRDS(
  Ssmooth,
  file.path(output_dir, "GGE_EFA_smoothed_S.rds")
)

saveRDS(
  efa_models,
  file.path(output_dir, "GGE_EFA_models.rds")
)

# Save loadings and explained variance for each model.
for (k in 1:4) {

  fit <- efa_models[[paste0("EFA_", k)]]

  L <- as.matrix(fit$loadings)

  write.csv(
    data.frame(
      Phenotype = rownames(L),
      L,
      check.names = FALSE
    ),
    file.path(
      output_dir,
      paste0("GGE_EFA_", k, "_loadings.csv")
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
      paste0("GGE_EFA_", k, "_variance.csv")
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
    paste("JME heritability:", jme_h2),
    paste(
      "Maximum smoothing adjustment:",
      smoothing_adjustment
    ),
    paste(
      "Minimum smoothed eigenvalue:",
      minimum_eigenvalue
    )
  ),
  file.path(output_dir, "GGE_EFA_provenance.txt")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "GGE_EFA_session_info.txt")
)

cat("\nResults saved to:", output_dir, "\n")
cat("GGE EFA completed.\n")
