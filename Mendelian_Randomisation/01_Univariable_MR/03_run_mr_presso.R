
# ============================================================
# CHAPTER 5: FORWARD UNIVARIABLE MR-PRESSO
# ============================================================
#
# Four exposures and ten epilepsy outcomes (40 comparisons).
#
# Original analysis settings:
#   OUTLIERtest     = TRUE
#   DISTORTIONtest  = TRUE
#   NbDistribution  = 10000
#   SignifThreshold = 0.05
#
# Uses archived clumped exposure instruments and cleaned
# ILAE 2023 epilepsy GWAS outcome datasets.
#
# Reconstructed workflow:
#   - TwoSampleMR allele harmonisation (action = 2)
#   - MR-PRESSO global, outlier and distortion tests
#   - Fixed, comparison-specific random seeds
#   - Individual saved results for restartable analyses
#   - Consolidated global-test results and session information
#
# Original input datasets are never modified.
# ============================================================


# ------------------------------------------------------------
# 1. REQUIRE PACKAGES
# ------------------------------------------------------------

if (!requireNamespace("TwoSampleMR", quietly = TRUE)) {
  stop("TwoSampleMR is required.")
}

if (!requireNamespace("MRPRESSO", quietly = TRUE)) {
  stop("MRPRESSO is required.")
}


# ------------------------------------------------------------
# 2. CONFIGURE DIRECTORIES
# ------------------------------------------------------------

# Set MR_DATA_DIR before running.
#
# Example:
# Sys.setenv(MR_DATA_DIR = "/path/to/MR_data")

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
  "Univariable_MR_Reproduction",
  "MR_PRESSO"
)

raw_dir <- file.path(
  output_dir,
  "raw_results"
)

dir.create(
  raw_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 3. EXPOSURE INSTRUMENT PATHS
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


# ------------------------------------------------------------
# 4. EPILEPSY OUTCOME PATHS
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
  !anyDuplicated(names(exposure_paths)),
  !anyDuplicated(names(outcome_paths))
)

all_paths <- c(
  exposure_paths,
  outcome_paths
)

missing_files <- all_paths[
  !file.exists(all_paths)
]

if (length(missing_files) > 0L) {
  stop(
    "Missing input files:\n",
    paste(missing_files, collapse = "\n")
  )
}

cat("\nAll 14 input datasets found.\n")


# ------------------------------------------------------------
# 6. SELECT EXPOSURES
# ------------------------------------------------------------

# Optional environment variable:
#
# Sys.setenv(MR_PRESSO_EXPOSURES = "MDD")
#
# Or run all exposures:
#
# Sys.unsetenv("MR_PRESSO_EXPOSURES")

selected <- Sys.getenv(
  "MR_PRESSO_EXPOSURES",
  unset = ""
)

selected_names <- if (nzchar(selected)) {

  trimws(
    strsplit(selected, ",", fixed = TRUE)[[1]]
  )

} else {

  names(exposure_paths)

}

if (
  length(selected_names) == 0L ||
  anyDuplicated(selected_names) ||
  !all(selected_names %in% names(exposure_paths))
) {
  stop("Invalid MR_PRESSO_EXPOSURES selection.")
}

cat(
  "Selected exposures:",
  paste(selected_names, collapse = ", "),
  "\n"
)


# ------------------------------------------------------------
# 7. LOAD EXPOSURE INSTRUMENTS
# ------------------------------------------------------------

exposures <- lapply(
  exposure_paths[selected_names],
  readRDS
)

required_exposure_columns <- c(
  "SNP",
  "beta.exposure",
  "se.exposure",
  "effect_allele.exposure",
  "other_allele.exposure"
)

for (exposure_name in names(exposures)) {

  dat <- exposures[[exposure_name]]

  if (
    !is.data.frame(dat) ||
    !all(required_exposure_columns %in% names(dat)) ||
    anyDuplicated(dat$SNP)
  ) {
    stop(
      "Invalid exposure dataset: ",
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

instrument_snps <- unique(
  unlist(
    lapply(
      exposures,
      function(dat) as.character(dat$SNP)
    ),
    use.names = FALSE
  )
)


# ------------------------------------------------------------
# 8. ORIGINAL MR-PRESSO SETTINGS
# ------------------------------------------------------------

nb_distribution <- 10000L

significance_threshold <- 0.05

# Comparison-specific seeds, independent of run order.
# The formula matches the validated local analysis.

all_exposure_names <- c(
  "Educational_attainment",
  "T2D",
  "Stroke",
  "MDD"
)

seed_for_pair <- function(
  exposure_name,
  outcome_name
) {

  e <- match(
    exposure_name,
    all_exposure_names
  )

  o <- match(
    outcome_name,
    names(outcome_paths)
  )

  if (is.na(e) || is.na(o)) {
    stop("Unknown exposure or outcome.")
  }

  as.integer(
    20261029L + 10L * (o - 1L) + e
  )
}


# ------------------------------------------------------------
# 9. VALIDATE A SAVED ANALYSIS
# ------------------------------------------------------------

# Saved analyses are reused only when their metadata
# and MR-PRESSO global-test component are valid.

read_validated_result <- function(
  path,
  exposure_name,
  outcome_name
) {

  saved <- readRDS(path)

  expected_seed <- seed_for_pair(
    exposure_name,
    outcome_name
  )

  global <- saved$result[["MR-PRESSO results"]][["Global Test"]]

  if (
    !identical(saved$exposure, exposure_name) ||
    !identical(saved$outcome, outcome_name) ||
    !identical(saved$seed, expected_seed) ||
    !identical(
      saved$nb_distribution,
      nb_distribution
    ) ||
    !isTRUE(all.equal(
      saved$significance_threshold,
      significance_threshold
    )) ||
    length(saved$nsnp) != 1L ||
    is.na(saved$nsnp) ||
    saved$nsnp < 4L ||
    is.null(global) ||
    is.null(global$RSSobs) ||
    is.null(global$Pvalue)
  ) {
    stop(
      "Invalid or incompatible saved analysis:\n",
      path
    )
  }

  saved
}


# ------------------------------------------------------------
# 10. RUN MR-PRESSO ANALYSES
# ------------------------------------------------------------

for (outcome_name in names(outcome_paths)) {

  cat("\n====================================\n")
  cat("OUTCOME:", outcome_name, "\n")
  cat("====================================\n")

  expected_files <- file.path(
    raw_dir,
    paste0(
      names(exposures),
      "__",
      outcome_name,
      ".rds"
    )
  )

  # If every selected analysis is already complete,
  # validate the cached results without loading the
  # large epilepsy GWAS dataset again.

  if (all(file.exists(expected_files))) {

    for (exposure_name in names(exposures)) {

      raw_file <- file.path(
        raw_dir,
        paste0(
          exposure_name,
          "__",
          outcome_name,
          ".rds"
        )
      )

      read_validated_result(
        raw_file,
        exposure_name,
        outcome_name
      )

      cat(
        "Already completed:",
        exposure_name,
        "->",
        outcome_name,
        "\n"
      )
    }

    next
  }

  # ----------------------------------------------------------
  # Load one large GWAS outcome at a time
  # ----------------------------------------------------------

  outcome_full <- readRDS(
    outcome_paths[[outcome_name]]
  )

  outcome_subset <- outcome_full[
    outcome_full$SNP %in% instrument_snps,
    ,
    drop = FALSE
  ]

  rm(outcome_full)
  invisible(gc())

  # ----------------------------------------------------------
  # Analyse each selected exposure
  # ----------------------------------------------------------

  for (exposure_name in names(exposures)) {

    raw_file <- file.path(
      raw_dir,
      paste0(
        exposure_name,
        "__",
        outcome_name,
        ".rds"
      )
    )

    seed <- seed_for_pair(
      exposure_name,
      outcome_name
    )

    # Reuse previously completed result
    if (file.exists(raw_file)) {

      read_validated_result(
        raw_file,
        exposure_name,
        outcome_name
      )

      cat(
        "Already completed:",
        exposure_name,
        "->",
        outcome_name,
        "\n"
      )

      next
    }

    exposure <- exposures[[exposure_name]]

    matched_outcome <- outcome_subset[
      outcome_subset$SNP %in% exposure$SNP,
      ,
      drop = FALSE
    ]

    if (nrow(matched_outcome) == 0L) {
      stop(
        "No overlapping SNPs: ",
        exposure_name,
        " -> ",
        outcome_name
      )
    }

    # --------------------------------------------------------
    # Harmonise exposure and outcome alleles
    # --------------------------------------------------------

    dat <- TwoSampleMR::harmonise_data(
      exposure_dat = exposure,
      outcome_dat = matched_outcome,
      action = 2
    )

    dat <- dat[
      !is.na(dat$mr_keep) & dat$mr_keep,
      ,
      drop = FALSE
    ]

    nsnp <- length(unique(dat$SNP))

    if (
      nsnp < 4L ||
      nsnp != nrow(dat)
    ) {
      stop(
        "Invalid harmonised instrument count: ",
        exposure_name,
        " -> ",
        outcome_name
      )
    }

    # --------------------------------------------------------
    # Validate effect and standard error columns
    # --------------------------------------------------------

    required_numeric <- c(
      "beta.exposure",
      "beta.outcome",
      "se.exposure",
      "se.outcome"
    )

    if (
      !all(required_numeric %in% names(dat)) ||
      any(!is.finite(as.matrix(
        dat[, required_numeric, drop = FALSE]
      ))) ||
      any(dat$se.exposure <= 0) ||
      any(dat$se.outcome <= 0)
    ) {
      stop(
        "Invalid MR-PRESSO input values: ",
        exposure_name,
        " -> ",
        outcome_name
      )
    }

    # Preserve SNP identifiers in result row names
    rownames(dat) <- dat$SNP

    cat(
      "\nRunning:",
      exposure_name,
      "->",
      outcome_name,
      "| SNPs:",
      nsnp,
      "\n"
    )

    # --------------------------------------------------------
    # Run MR-PRESSO
    # --------------------------------------------------------

    set.seed(seed)

    presso <- MRPRESSO::mr_presso(
      BetaOutcome = "beta.outcome",
      BetaExposure = "beta.exposure",
      SdOutcome = "se.outcome",
      SdExposure = "se.exposure",
      OUTLIERtest = TRUE,
      DISTORTIONtest = TRUE,
      data = dat,
      NbDistribution = nb_distribution,
      SignifThreshold = significance_threshold
    )

    global <- presso[["MR-PRESSO results"]][["Global Test"]]

    if (
      is.null(global) ||
      is.null(global$RSSobs) ||
      is.null(global$Pvalue)
    ) {
      stop(
        "Missing MR-PRESSO global-test result: ",
        exposure_name,
        " -> ",
        outcome_name
      )
    }

    # --------------------------------------------------------
    # Preserve the complete MR-PRESSO result
    # --------------------------------------------------------

    # This retains:
    #   Global Test
    #   Outlier Test, where available
    #   Distortion Test, where available
    #   Main MR results, including corrected estimates
    #   where those estimates are available

    saved_result <- list(
      exposure = exposure_name,
      outcome = outcome_name,
      nsnp = nsnp,
      seed = seed,
      nb_distribution = nb_distribution,
      significance_threshold = significance_threshold,
      result = presso
    )

    saveRDS(
      saved_result,
      raw_file
    )

    cat(
      "Global RSS:",
      global$RSSobs,
      "| Global p:",
      as.character(global$Pvalue),
      "| Saved\n"
    )

    rm(
      dat,
      matched_outcome,
      presso,
      saved_result
    )
  }

  rm(outcome_subset)
  invisible(gc())
}


# ------------------------------------------------------------
# 11. COMBINE COMPLETED RESULTS
# ------------------------------------------------------------

# Read only the exact files expected for this run,
# avoiding duplicates or unrelated RDS files.

expected_pairs <- expand.grid(
  exposure = selected_names,
  outcome = names(outcome_paths),
  stringsAsFactors = FALSE
)

expected_pairs$file <- file.path(
  raw_dir,
  paste0(
    expected_pairs$exposure,
    "__",
    expected_pairs$outcome,
    ".rds"
  )
)

missing_results <- expected_pairs$file[
  !file.exists(expected_pairs$file)
]

if (length(missing_results) > 0L) {
  stop(
    "Missing MR-PRESSO analyses:\n",
    paste(missing_results, collapse = "\n")
  )
}

global_results <- do.call(
  rbind,
  lapply(seq_len(nrow(expected_pairs)), function(i) {

    exposure_name <- expected_pairs$exposure[i]
    outcome_name <- expected_pairs$outcome[i]

    saved <- read_validated_result(
      expected_pairs$file[i],
      exposure_name,
      outcome_name
    )

    # Correct R syntax for nested list extraction
    global <- saved$result[["MR-PRESSO results"]][["Global Test"]]

    data.frame(
      exposure_label = saved$exposure,
      outcome_label = saved$outcome,
      nsnp = saved$nsnp,
      rssobs = as.numeric(global$RSSobs),
      global_pval = as.character(global$Pvalue),
      seed = saved$seed,
      nb_distribution = saved$nb_distribution,
      stringsAsFactors = FALSE
    )
  })
)


# ------------------------------------------------------------
# 12. VALIDATE COMPLETENESS
# ------------------------------------------------------------

expected_count <- length(selected_names) * 10L

result_keys <- paste(
  global_results$exposure_label,
  global_results$outcome_label,
  sep = "|"
)

if (
  nrow(global_results) != expected_count ||
  anyDuplicated(result_keys) ||
  any(table(global_results$exposure_label) != 10L)
) {
  stop(
    "MR-PRESSO results failed completeness validation."
  )
}

cat(
  "\nCompleted MR-PRESSO comparisons:",
  nrow(global_results),
  "\n"
)


# ------------------------------------------------------------
# 13. SAVE CONSOLIDATED RESULTS
# ------------------------------------------------------------

# For partial runs, use a separate filename to avoid
# overwriting the complete 40-comparison results.

output_filename <- if (length(selected_names) == 4L) {

  "08_MR_PRESSO_global_tests.csv"

} else {

  paste0(
    "08_MR_PRESSO_global_tests_",
    paste(selected_names, collapse = "_"),
    ".csv"
  )

}

write.csv(
  global_results,
  file.path(
    output_dir,
    output_filename
  ),
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(
    output_dir,
    "09_MR_PRESSO_sessionInfo.txt"
  )
)


# ------------------------------------------------------------
# 14. FINAL SUMMARY
# ------------------------------------------------------------

cat("\n====================================\n")
cat("MR-PRESSO RUN COMPLETE\n")
cat("====================================\n")

cat(
  "Completed comparisons:",
  nrow(global_results),
  "of",
  expected_count,
  "\n"
)

cat("\nGlobal tests by exposure:\n")

print(
  table(global_results$exposure_label)
)

cat("\nMR-PRESSO GLOBAL TEST RESULTS:\n")

print(
  global_results[
    ,
    c(
      "exposure_label",
      "outcome_label",
      "nsnp",
      "rssobs",
      "global_pval"
    ),
    drop = FALSE
  ],
  row.names = FALSE
)

cat("\nResults directory:\n", output_dir, "\n")

cat("\nMR-PRESSO reproduction completed successfully.\n")
