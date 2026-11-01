
# ============================================================
# CHAPTER 5: FORWARD UNIVARIABLE MENDELIAN RANDOMISATION
# ============================================================
#
# Thesis:
# From genetic risk to brain structure:
# An integrative investigation of epilepsy and
# functional dissociative seizures
#
# Purpose:
# Reproduce the forward two-sample univariable MR analyses
# using archived, clumped exposure instruments and cleaned
# ILAE 2023 epilepsy GWAS datasets.
#
# Exposures (n = 4):
#   Educational attainment
#   Type 2 diabetes
#   Stroke
#   Major depressive disorder
#
# Epilepsy outcomes (n = 10):
#   All epilepsy, FE, GGE, F-HS, F-Neg, F-Oth,
#   GTCSA, JAE, JME, CAE
#
# MR methods:
#   Inverse variance weighted (IVW)
#   MR-Egger
#   Weighted median
#   MR-RAPS
#
# Total:
#   40 exposure-outcome comparisons
#   160 MR estimates
#
# Requirements:
#   R
#   TwoSampleMR
#   Archived input .rds datasets
#
# Input datasets are read-only.
# This script does not regenerate GWAS summary statistics
# or perform the original LD clumping.
#
# ============================================================


# ------------------------------------------------------------
# 1. LOAD PACKAGE
# ------------------------------------------------------------

if (!requireNamespace("TwoSampleMR", quietly = TRUE)) {
  stop(
    "The TwoSampleMR package is required. ",
    "Please install it before running this script."
  )
}


# ------------------------------------------------------------
# 2. CONFIGURE DATA DIRECTORY
# ------------------------------------------------------------

# Set MR_DATA_DIR to the location containing the original
# exposure instrument and cleaned epilepsy outcome datasets.
#
# Example in an RStudio Console:
#
# Sys.setenv(
#   MR_DATA_DIR = "D:/External_backup/PHD/5. Two sample MR"
# )
#
# The original GWAS datasets and .rds files are not
# distributed with the GitHub repository.

root <- Sys.getenv("MR_DATA_DIR", unset = "")

if (!nzchar(root)) {
  stop(
    "MR_DATA_DIR is not set.\n",
    "Set it to the directory containing the original ",
    "MR input datasets."
  )
}

root <- normalizePath(
  root,
  winslash = "/",
  mustWork = TRUE
)

output_dir <- file.path(
  root,
  "Results",
  "Univariable_MR_Reproduction"
)

cat("\n=== DATA CONFIGURATION ===\n")
cat("Input directory:", root, "\n")
cat("Output directory:", output_dir, "\n")


# ------------------------------------------------------------
# 3. DEFINE EXPOSURE INSTRUMENT DATASETS
# ------------------------------------------------------------

exposure_paths <- c(

  Educational_attainment = file.path(
    root,
    "educational_attainment_ssgac_2022_eur",
    "EA4_additive_excl_23andMe.txt",
    "Education_EUR_Instruments_Clumped_P8.rds"
  ),

  T2D = file.path(
    root,
    "type2_diabetes_diagram_2024_eur",
    "T2D_EUR_Suzuki2024_Instruments_Clumped_P5e8.rds"
  ),

  Stroke = file.path(
    root,
    "gigastroke_2022_eur",
    "Stroke_EUR_Instruments_Clumped.rds"
  ),

  MDD = file.path(
    root,
    "mdd_pgc_2025_eur",
    "pgc-mdd2025_no23andMe-noUKBB_eur_v3-49-24-11.tsv",
    "MDD_EUR_Instruments_Clumped_P5e8.rds"
  )
)

# Explicitly retain exposure labels
exposure_paths <- setNames(
  unname(exposure_paths),
  c(
    "Educational_attainment",
    "T2D",
    "Stroke",
    "MDD"
  )
)


# ------------------------------------------------------------
# 4. DEFINE EPILEPSY OUTCOME DATASETS
# ------------------------------------------------------------

outcome_files <- c(

  All_epilepsy =
    "ILAE_All_Epilepsy_Cleaned.rds",

  FE =
    "ILAE_Focal_Epilepsy_Cleaned.rds",

  GGE =
    "ILAE_GGE_Epilepsy_Cleaned.rds",

  Focal_HS =
    "ILAE_Focal_HS_Epilepsy_Cleaned.rds",

  Focal_lesion_negative =
    "ILAE_Focal_Lesion_Negative_Epilepsy_Cleaned.rds",

  Focal_other_lesion =
    "ILAE_Focal_Other_Lesion_Epilepsy_Cleaned.rds",

  GTCSA =
    "ILAE_GTCS_Epilepsy_Cleaned.rds",

  JAE =
    "ILAE_JAE_Epilepsy_Cleaned.rds",

  JME =
    "ILAE_JME_Epilepsy_Cleaned.rds",

  CAE =
    "ILAE_CAE_Epilepsy_Cleaned.rds"
)

# Preserve outcome labels explicitly
outcome_paths <- setNames(
  file.path(root, unname(outcome_files)),
  names(outcome_files)
)


# ------------------------------------------------------------
# 5. VALIDATE CONFIGURATION AND INPUT FILES
# ------------------------------------------------------------

cat("\n=== INPUT VALIDATION ===\n")

cat(
  "Number of exposures:",
  length(exposure_paths),
  "\n"
)

cat(
  "Number of outcomes:",
  length(outcome_paths),
  "\n"
)

if (
  length(exposure_paths) != 4L ||
  length(outcome_paths) != 10L ||
  is.null(names(exposure_paths)) ||
  is.null(names(outcome_paths)) ||
  anyNA(names(exposure_paths)) ||
  anyNA(names(outcome_paths)) ||
  any(!nzchar(names(exposure_paths))) ||
  any(!nzchar(names(outcome_paths)))
) {
  stop("Exposure or outcome configuration is invalid.")
}

all_paths <- c(
  exposure_paths,
  outcome_paths
)

missing_paths <- all_paths[
  !file.exists(all_paths)
]

if (length(missing_paths) > 0L) {
  stop(
    "The following input datasets are missing:\n",
    paste(
      missing_paths,
      collapse = "\n"
    )
  )
}

cat("All 14 input datasets found.\n")


# ------------------------------------------------------------
# 6. CREATE OUTPUT DIRECTORY
# ------------------------------------------------------------

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 7. LOAD EXPOSURE INSTRUMENTS
# ------------------------------------------------------------

exposures <- lapply(
  exposure_paths,
  readRDS
)

required_exposure_columns <- c(
  "SNP",
  "beta.exposure",
  "se.exposure",
  "effect_allele.exposure",
  "other_allele.exposure",
  "id.exposure"
)

cat("\n=== EXPOSURE INSTRUMENTS ===\n")

for (exposure_name in names(exposures)) {

  dat <- exposures[[exposure_name]]

  if (
    !is.data.frame(dat) ||
    !all(required_exposure_columns %in% names(dat))
  ) {
    stop(
      "Invalid exposure dataset: ",
      exposure_name
    )
  }

  if (anyDuplicated(dat$SNP)) {
    stop(
      "Duplicate SNPs detected in exposure: ",
      exposure_name
    )
  }

  cat(
    exposure_name,
    ":",
    nrow(dat),
    "instruments\n"
  )
}


# ------------------------------------------------------------
# 8. COLLECT ALL EXPOSURE SNPs
# ------------------------------------------------------------

instrument_snps <- unique(
  unlist(
    lapply(
      exposures,
      function(dat) as.character(dat$SNP)
    ),
    use.names = FALSE
  )
)

cat(
  "\nUnique SNPs across all exposures:",
  length(instrument_snps),
  "\n"
)


# ------------------------------------------------------------
# 9. DEFINE MR METHODS
# ------------------------------------------------------------

mr_methods <- c(
  "mr_ivw",
  "mr_egger_regression",
  "mr_weighted_median",
  "mr_raps"
)


# ------------------------------------------------------------
# 10. INITIALISE RESULTS
# ------------------------------------------------------------

results_list <- list()
log_list <- list()

pair_number <- 0L


# ------------------------------------------------------------
# 11. RUN FORWARD UNIVARIABLE MR
# ------------------------------------------------------------

# Process one epilepsy GWAS dataset at a time.
# This avoids retaining all ten large outcome datasets
# simultaneously in memory.

for (outcome_name in names(outcome_paths)) {

  cat(
    "\n====================================\n",
    "OUTCOME:", outcome_name, "\n",
    "====================================\n"
  )

  # ----------------------------------------------------------
  # Load epilepsy GWAS
  # ----------------------------------------------------------

  outcome_data <- readRDS(
    outcome_paths[[outcome_name]]
  )

  required_outcome_columns <- c(
    "SNP",
    "beta.outcome",
    "se.outcome",
    "effect_allele.outcome",
    "other_allele.outcome",
    "id.outcome"
  )

  if (
    !is.data.frame(outcome_data) ||
    !all(required_outcome_columns %in% names(outcome_data))
  ) {
    stop(
      "Invalid outcome dataset: ",
      outcome_name
    )
  }

  # ----------------------------------------------------------
  # Extract SNPs required by any exposure
  # ----------------------------------------------------------

  outcome_subset <- outcome_data[
    outcome_data$SNP %in% instrument_snps,
    ,
    drop = FALSE
  ]

  # Free large GWAS object from memory
  rm(outcome_data)
  invisible(gc())

  if (nrow(outcome_subset) == 0L) {
    stop(
      "No matching SNPs found for outcome: ",
      outcome_name
    )
  }

  # ----------------------------------------------------------
  # Analyse each exposure against this outcome
  # ----------------------------------------------------------

  for (exposure_name in names(exposures)) {

    pair_number <- pair_number + 1L

    cat(
      "\nComparison", pair_number, "of 40:",
      exposure_name, "->", outcome_name, "\n"
    )

    exposure_data <- exposures[[exposure_name]]

    # Match outcome SNPs to exposure instruments
    matched_outcome <- outcome_subset[
      outcome_subset$SNP %in% exposure_data$SNP,
      ,
      drop = FALSE
    ]

    if (nrow(matched_outcome) == 0L) {
      stop(
        "No overlapping SNPs for ",
        exposure_name, " -> ", outcome_name
      )
    }

    # --------------------------------------------------------
    # Harmonise exposure and outcome alleles
    # --------------------------------------------------------

    harmonised <- TwoSampleMR::harmonise_data(
      exposure_dat = exposure_data,
      outcome_dat = matched_outcome,
      action = 2
    )

    if (
      !is.data.frame(harmonised) ||
      !("mr_keep" %in% names(harmonised))
    ) {
      stop(
        "Harmonisation failed for ",
        exposure_name, " -> ", outcome_name
      )
    }

    # Retain only eligible harmonised instruments
    harmonised <- harmonised[
      !is.na(harmonised$mr_keep) &
        harmonised$mr_keep,
      ,
      drop = FALSE
    ]

    nsnp <- length(
      unique(harmonised$SNP)
    )

    if (nsnp < 3L) {
      stop(
        "Fewer than 3 harmonised SNPs for ",
        exposure_name, " -> ", outcome_name,
        " (n = ", nsnp, ")"
      )
    }

    # --------------------------------------------------------
    # Estimate MR effects
    # --------------------------------------------------------

    # Fixed seed to make bootstrap-dependent estimates
    # reproducible within the same computational environment.
    #
    # Historical weighted-median standard errors and p-values
    # may differ because original bootstrap random states
    # were not preserved.

    set.seed(
      20261029L + pair_number
    )

    result <- TwoSampleMR::mr(
      harmonised,
      method_list = mr_methods
    )

    result <- TwoSampleMR::generate_odds_ratios(
      result
    )

    # Require all four intended methods
    if (
      !is.data.frame(result) ||
      nrow(result) != 4L ||
      anyNA(result$b) ||
      anyNA(result$se) ||
      anyNA(result$pval)
    ) {
      stop(
        "Incomplete MR results for ",
        exposure_name, " -> ", outcome_name
      )
    }

    # Add consistent labels
    result$exposure_label <- exposure_name
    result$outcome_label <- outcome_name

    # Store estimated effects
    results_list[[pair_number]] <- result

    # --------------------------------------------------------
    # Record instrument flow
    # --------------------------------------------------------

    log_list[[pair_number]] <- data.frame(
      exposure = exposure_name,
      outcome = outcome_name,
      n_exposure_instruments = nrow(exposure_data),
      n_matching_outcome_records = nrow(matched_outcome),
      n_harmonised_instruments = nsnp,
      n_methods = nrow(result),
      stringsAsFactors = FALSE
    )

    cat(
      "Retained SNPs:", nsnp,
      "| MR methods:", nrow(result),
      "| Success\n"
    )
  }

  # Release outcome subset before next GWAS
  rm(outcome_subset)
  invisible(gc())
}


# ------------------------------------------------------------
# 12. VALIDATE COMPLETENESS
# ------------------------------------------------------------

if (
  length(results_list) != 40L ||
  length(log_list) != 40L
) {
  stop(
    "Expected 40 completed comparisons, but found ",
    length(results_list)
  )
}

all_results <- do.call(
  rbind,
  results_list
)

analysis_log <- do.call(
  rbind,
  log_list
)

if (
  !is.data.frame(all_results) ||
  !is.data.frame(analysis_log) ||
  nrow(all_results) != 160L ||
  nrow(analysis_log) != 40L
) {
  stop(
    "Validation failed: expected 160 MR estimates ",
    "and 40 analysis-log records."
  )
}

# Every exposure-outcome comparison must contain
# exactly four MR method estimates.

pair_counts <- table(
  all_results$exposure_label,
  all_results$outcome_label
)

if (
  !identical(dim(pair_counts), c(4L, 10L)) ||
  any(pair_counts != 4L)
) {
  stop(
    "Validation failed: incomplete exposure-outcome ",
    "method combinations."
  )
}

cat(
  "\nAll 40 comparisons passed completeness checks.\n"
)


# ------------------------------------------------------------
# 13. SAVE REPRODUCED MR RESULTS
# ------------------------------------------------------------

output_columns <- c(
  "exposure_label",
  "outcome_label",
  "method",
  "nsnp",
  "b",
  "se",
  "pval",
  "or",
  "or_lci95",
  "or_uci95"
)

write.csv(
  all_results[
    ,
    output_columns,
    drop = FALSE
  ],
  file.path(
    output_dir,
    "01_univariable_MR_results.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 14. SAVE INSTRUMENT-FLOW LOG
# ------------------------------------------------------------

write.csv(
  analysis_log,
  file.path(
    output_dir,
    "02_univariable_MR_run_log.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 15. SAVE SOFTWARE ENVIRONMENT
# ------------------------------------------------------------

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    output_dir,
    "03_sessionInfo.txt"
  )
)


# ------------------------------------------------------------
# 16. FINAL SUMMARY
# ------------------------------------------------------------

cat("\n====================================\n")
cat("UNIVARIABLE MR REPRODUCTION COMPLETE\n")
cat("====================================\n")

cat(
  "Exposure-outcome comparisons:",
  nrow(analysis_log),
  "\n"
)

cat(
  "Total MR estimates:",
  nrow(all_results),
  "\n"
)

cat(
  "Results saved in:\n",
  output_dir,
  "\n"
)

cat("\nResults by exposure and MR method:\n")

print(
  table(
    all_results$exposure_label,
    all_results$method
  )
)

cat("\nPrimary IVW results:\n")

print(
  all_results[
    all_results$method == "Inverse variance weighted",
    c(
      "exposure_label",
      "outcome_label",
      "nsnp",
      "or",
      "or_lci95",
      "or_uci95",
      "pval"
    ),
    drop = FALSE
  ],
  row.names = FALSE
)

cat("\nAnalysis completed successfully.\n")
