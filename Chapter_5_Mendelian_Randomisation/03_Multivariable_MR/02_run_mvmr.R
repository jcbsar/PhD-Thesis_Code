
# Chapter 5: Multivariable Mendelian Randomisation
# Educational attainment and type 2 diabetes -> epilepsy
#
# R 4.4.2
# TwoSampleMR 0.7.5
#
# Final analysis: joint LD-clumped instrument set
# IVW MVMR, no intercept, instrument_specific = FALSE
#
# Requires:
#   - Six harmonised_MVMR_data_*.rds files
#   - 09_MVMR_results_all_outcomes.csv
#
# Set MVMR_DATA_DIR before running this script.

library(TwoSampleMR)

# 1. Input directory ---------------------------------------

data_dir <- Sys.getenv("MVMR_DATA_DIR", unset = "")

if (!nzchar(data_dir) || !dir.exists(data_dir)) {
  stop("Set MVMR_DATA_DIR to the final JOINT_CLUMPED archive.")
}

# 2. Epilepsy outcomes -------------------------------------

outcomes <- c(
  "All_Epilepsy",
  "Focal_Epilepsy",
  "GGE_Epilepsy",
  "FocalHS_Epilepsy",
  "FocalLesionNegative_Epilepsy",
  "JME_Epilepsy"
)

# 3. Original results for validation -----------------------

reference <- read.csv(
  file.path(data_dir, "09_MVMR_results_all_outcomes.csv"),
  stringsAsFactors = FALSE
)

same <- function(actual, expected) {
  isTRUE(all.equal(
    as.numeric(actual),
    as.numeric(expected),
    tolerance = 1e-8,
    check.attributes = FALSE
  ))
}

# 4. Reproduce MVMR models --------------------------------

results <- lapply(outcomes, function(outcome) {

  input_file <- file.path(
    data_dir,
    paste0("harmonised_MVMR_data_", outcome, ".rds")
  )

  if (!file.exists(input_file)) {
    stop("Missing input: ", input_file)
  }

  dat <- readRDS(input_file)
  n_snps <- nrow(dat$exposure_beta)

  fit <- TwoSampleMR::mv_multiple(
    dat,
    intercept = FALSE,
    instrument_specific = FALSE,
    plots = FALSE
  )

  actual <- fit$result

  expected <- reference[
    reference$outcome_label == outcome, ,
    drop = FALSE
  ]

  # Match rows by exposure, independent of row order
  idx <- match(expected$exposure, actual$exposure)

  if (nrow(expected) != 2L || anyNA(idx)) {
    stop("Exposure mismatch for ", outcome)
  }

  actual <- actual[idx, , drop = FALSE]

  # Validate against the original final analysis
  checks <- c(
    SNPs = all(
      n_snps == expected$n_instruments_in_MVMR_model
    ),
    beta = same(actual$b, expected$b),
    SE = same(actual$se, expected$se),
    pvalue = same(actual$pval, expected$pval),
    OR = same(exp(actual$b), expected$OR)
  )

  if (!all(checks)) {
    stop(
      "Validation failed for ", outcome, ": ",
      paste(names(checks)[!checks], collapse = ", ")
    )
  }

  message(
    outcome, ": ", n_snps,
    " SNPs | All validation checks passed"
  )

  data.frame(
    outcome = outcome,
    exposure = as.character(actual$exposure),
    nsnp = n_snps,
    beta = actual$b,
    se = actual$se,
    pvalue = actual$pval,
    OR = exp(actual$b),
    OR_lci95 = exp(actual$b - 1.96 * actual$se),
    OR_uci95 = exp(actual$b + 1.96 * actual$se),
    bonferroni_significant = actual$pval < 0.05 / 12,
    row.names = NULL
  )
})

# 5. Combine and export -----------------------------------

results <- do.call(rbind, results)

output_dir <- "outputs"
dir.create(output_dir, showWarnings = FALSE)

write.csv(
  results,
  file.path(output_dir, "MVMR_results_reproduced.csv"),
  row.names = FALSE
)

cat("\nAll six MVMR models reproduced successfully.\n")
cat("12 exposure-effect estimates validated.\n\n")

print(results, row.names = FALSE)

# Preserve software environment for reproducibility
writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "MVMR_sessionInfo.txt")
)
