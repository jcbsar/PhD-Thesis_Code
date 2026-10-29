
# ============================================================
# CHAPTER 5: UNIVARIABLE MR DIAGNOSTICS
# Cochran's Q and MR-Egger intercept tests
# ============================================================
#
# Reproduces diagnostic analyses for:
#   4 exposures x 10 epilepsy outcomes = 40 comparisons
#
# Outputs:
#   80 Cochran's Q estimates (IVW and MR-Egger)
#   40 MR-Egger intercept estimates
#   40 instrument-flow records
#
# Uses archived clumped exposure instruments and
# cleaned epilepsy GWAS datasets.
#
# Does not perform MR-PRESSO.
# ============================================================


# ------------------------------------------------------------
# 1. REQUIRE PACKAGE
# ------------------------------------------------------------

if (!requireNamespace("TwoSampleMR", quietly = TRUE)) {
  stop("TwoSampleMR must be installed.")
}


# ------------------------------------------------------------
# 2. CONFIGURE DIRECTORIES
# ------------------------------------------------------------

root <- Sys.getenv("MR_DATA_DIR", unset = "")

if (!nzchar(root)) {
  stop("MR_DATA_DIR is not set.")
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


# ------------------------------------------------------------
# 3. EXPOSURE DATASETS
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

exposure_paths <- setNames(
  unname(exposure_paths),
  c("Educational_attainment", "T2D", "Stroke", "MDD")
)


# ------------------------------------------------------------
# 4. EPILEPSY OUTCOME DATASETS
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

# Preserve names explicitly
outcome_paths <- setNames(
  file.path(root, unname(outcome_files)),
  names(outcome_files)
)


# ------------------------------------------------------------
# 5. VALIDATE INPUTS
# ------------------------------------------------------------

stopifnot(
  length(exposure_paths) == 4L,
  length(outcome_paths) == 10L,
  !is.null(names(exposure_paths)),
  !is.null(names(outcome_paths)),
  !anyDuplicated(names(exposure_paths)),
  !anyDuplicated(names(outcome_paths))
)

all_paths <- c(exposure_paths, outcome_paths)

missing_files <- all_paths[!file.exists(all_paths)]

if (length(missing_files) > 0L) {
  stop(
    "Missing input files:\n",
    paste(missing_files, collapse = "\n")
  )
}

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("\nAll 14 input datasets found.\n")


# ------------------------------------------------------------
# 6. LOAD EXPOSURE INSTRUMENTS
# ------------------------------------------------------------

exposures <- lapply(exposure_paths, readRDS)

instrument_snps <- unique(
  unlist(
    lapply(
      exposures,
      function(x) as.character(x$SNP)
    ),
    use.names = FALSE
  )
)


# ------------------------------------------------------------
# 7. INITIALISE RESULTS
# ------------------------------------------------------------

heterogeneity_results <- list()
pleiotropy_results <- list()
diagnostic_logs <- list()

pair_number <- 0L


# ------------------------------------------------------------
# 8. ANALYSE ALL 40 COMPARISONS
# ------------------------------------------------------------

# Load one large epilepsy outcome dataset at a time.

for (outcome_name in names(outcome_paths)) {

  cat("\n=== OUTCOME:", outcome_name, "===\n")

  outcome_data <- readRDS(
    outcome_paths[[outcome_name]]
  )

  # Extract outcome records matching any exposure SNP
  outcome_subset <- outcome_data[
    outcome_data$SNP %in% instrument_snps,
    ,
    drop = FALSE
  ]

  rm(outcome_data)
  invisible(gc())

  for (exposure_name in names(exposures)) {

    pair_number <- pair_number + 1L

    exposure_data <- exposures[[exposure_name]]

    matched_outcome <- outcome_subset[
      outcome_subset$SNP %in% exposure_data$SNP,
      ,
      drop = FALSE
    ]

    if (nrow(matched_outcome) == 0L) {
      stop(
        "No SNP overlap: ",
        exposure_name, " -> ", outcome_name
      )
    }

    # --------------------------------------------------------
    # Harmonise alleles
    # --------------------------------------------------------

    harmonised <- TwoSampleMR::harmonise_data(
      exposure_dat = exposure_data,
      outcome_dat = matched_outcome,
      action = 2
    )

    harmonised <- harmonised[
      !is.na(harmonised$mr_keep) &
        harmonised$mr_keep,
      ,
      drop = FALSE
    ]

    nsnp <- length(unique(harmonised$SNP))

    if (nsnp < 3L) {
      stop(
        "Insufficient SNPs: ",
        exposure_name, " -> ", outcome_name
      )
    }

    # --------------------------------------------------------
    # Cochran's Q heterogeneity tests
    # --------------------------------------------------------

    q_results <- TwoSampleMR::mr_heterogeneity(
      harmonised,
      method_list = c(
        "mr_ivw",
        "mr_egger_regression"
      )
    )

    # --------------------------------------------------------
    # MR-Egger intercept test
    # --------------------------------------------------------

    egger_results <- TwoSampleMR::mr_pleiotropy_test(
      harmonised
    )

    # --------------------------------------------------------
    # Validate results
    # --------------------------------------------------------

    if (
      !is.data.frame(q_results) ||
      nrow(q_results) != 2L ||
      !setequal(
        q_results$method,
        c("Inverse variance weighted", "MR Egger")
      ) ||
      anyNA(q_results$Q) ||
      anyNA(q_results$Q_pval)
    ) {
      stop(
        "Incomplete heterogeneity results: ",
        exposure_name, " -> ", outcome_name
      )
    }

    if (
      !is.data.frame(egger_results) ||
      nrow(egger_results) != 1L ||
      anyNA(egger_results$egger_intercept) ||
      anyNA(egger_results$se) ||
      anyNA(egger_results$pval)
    ) {
      stop(
        "Incomplete Egger intercept result: ",
        exposure_name, " -> ", outcome_name
      )
    }

    # --------------------------------------------------------
    # Add consistent analysis labels
    # --------------------------------------------------------

    q_results$exposure_label <- exposure_name
    q_results$outcome_label <- outcome_name

    egger_results$exposure_label <- exposure_name
    egger_results$outcome_label <- outcome_name

    heterogeneity_results[[pair_number]] <- q_results
    pleiotropy_results[[pair_number]] <- egger_results

    diagnostic_logs[[pair_number]] <- data.frame(
      exposure = exposure_name,
      outcome = outcome_name,
      n_exposure_instruments = nrow(exposure_data),
      n_matched_outcome_records = nrow(matched_outcome),
      n_harmonised_instruments = nsnp,
      status = "Success",
      stringsAsFactors = FALSE
    )

    cat(
      exposure_name,
      "->",
      outcome_name,
      "| SNPs:",
      nsnp,
      "| Q tests: 2",
      "| Egger intercept: 1",
      "| Success\n"
    )
  }

  rm(outcome_subset)
  invisible(gc())
}


# ------------------------------------------------------------
# 9. COMBINE RESULTS
# ------------------------------------------------------------

if (
  length(heterogeneity_results) != 40L ||
  length(pleiotropy_results) != 40L ||
  length(diagnostic_logs) != 40L
) {
  stop("Incomplete diagnostics: expected 40 comparisons.")
}

all_heterogeneity <- do.call(
  rbind,
  heterogeneity_results
)

all_pleiotropy <- do.call(
  rbind,
  pleiotropy_results
)

all_logs <- do.call(
  rbind,
  diagnostic_logs
)


# ------------------------------------------------------------
# 10. VERIFY COMPLETENESS
# ------------------------------------------------------------

if (
  nrow(all_heterogeneity) != 80L ||
  nrow(all_pleiotropy) != 40L ||
  nrow(all_logs) != 40L
) {
  stop(
    "Unexpected output counts. Expected ",
    "80 Q estimates, 40 intercepts, 40 logs."
  )
}

q_counts <- table(
  all_heterogeneity$exposure_label,
  all_heterogeneity$outcome_label
)

egger_counts <- table(
  all_pleiotropy$exposure_label,
  all_pleiotropy$outcome_label
)

if (
  !identical(dim(q_counts), c(4L, 10L)) ||
  !identical(dim(egger_counts), c(4L, 10L)) ||
  any(q_counts != 2L) ||
  any(egger_counts != 1L)
) {
  stop("One or more diagnostic comparisons are incomplete.")
}


# ------------------------------------------------------------
# 11. SAVE DIAGNOSTIC RESULTS
# ------------------------------------------------------------

write.csv(
  all_heterogeneity[
    ,
    c(
      "exposure_label",
      "outcome_label",
      "method",
      "Q",
      "Q_df",
      "Q_pval"
    ),
    drop = FALSE
  ],
  file.path(
    output_dir,
    "04_heterogeneity_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  all_pleiotropy[
    ,
    c(
      "exposure_label",
      "outcome_label",
      "egger_intercept",
      "se",
      "pval"
    ),
    drop = FALSE
  ],
  file.path(
    output_dir,
    "05_egger_intercept_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  all_logs,
  file.path(
    output_dir,
    "06_diagnostics_run_log.csv"
  ),
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(
    output_dir,
    "07_diagnostics_sessionInfo.txt"
  )
)


# ------------------------------------------------------------
# 12. FINAL SUMMARY
# ------------------------------------------------------------

cat("\n========================================\n")
cat("MR DIAGNOSTIC REPRODUCTION COMPLETE\n")
cat("========================================\n")

cat(
  "Exposure-outcome comparisons:",
  nrow(all_logs),
  "\n"
)

cat(
  "Cochran's Q estimates:",
  nrow(all_heterogeneity),
  "\n"
)

cat(
  "MR-Egger intercept estimates:",
  nrow(all_pleiotropy),
  "\n"
)

cat(
  "Results saved in:\n",
  output_dir,
  "\n"
)

cat("\nQ tests by exposure:\n")

print(
  table(
    all_heterogeneity$exposure_label,
    all_heterogeneity$method
  )
)

cat("\nMR-Egger intercept tests by exposure:\n")

print(
  table(
    all_pleiotropy$exposure_label
  )
)

cat("\nDiagnostics completed successfully.\n")
