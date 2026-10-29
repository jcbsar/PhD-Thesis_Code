
# ============================================================
# CHAPTER 5: FORWARD UNIVARIABLE MR-PRESSO
# ============================================================
# 4 exposures x 10 epilepsy outcomes
# Original settings:
#   OUTLIERtest = TRUE
#   DISTORTIONtest = TRUE
#   NbDistribution = 10000
#   SignifThreshold = 0.05
#
# Inputs: archived clumped exposure instruments and
#         cleaned ILAE 2023 outcome GWAS datasets
#
# Outputs:
#   MR-PRESSO global-test CSV
#   Individual raw result RDS files, preserving global,
#   outlier, distortion and MR-estimate results
#   Session information
#
# Fixed seeds are added for reproducible reruns.
# ============================================================

if (!requireNamespace("TwoSampleMR", quietly = TRUE) ||
    !requireNamespace("MRPRESSO", quietly = TRUE)) {
  stop("TwoSampleMR and MRPRESSO are required.")
}

root <- Sys.getenv("MR_DATA_DIR", unset = "")

if (!nzchar(root)) {
  stop("Set MR_DATA_DIR before running.")
}

root <- normalizePath(root, winslash = "/", mustWork = TRUE)

output_dir <- file.path(
  root, "Results", "Univariable_MR_Reproduction",
  "MR_PRESSO"
)

raw_dir <- file.path(output_dir, "raw_results")

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# INPUT PATHS
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

outcome_files <- c(
  All_epilepsy = "ILAE_All_Epilepsy_Cleaned.rds",
  FE = "ILAE_Focal_Epilepsy_Cleaned.rds",
  GGE = "ILAE_GGE_Epilepsy_Cleaned.rds",
  Focal_HS = "ILAE_Focal_HS_Epilepsy_Cleaned.rds",
  Focal_lesion_negative =
    "ILAE_Focal_Lesion_Negative_Epilepsy_Cleaned.rds",
  Focal_other_lesion =
    "ILAE_Focal_Other_Lesion_Epilepsy_Cleaned.rds",
  GTCSA = "ILAE_GTCS_Epilepsy_Cleaned.rds",
  JAE = "ILAE_JAE_Epilepsy_Cleaned.rds",
  JME = "ILAE_JME_Epilepsy_Cleaned.rds",
  CAE = "ILAE_CAE_Epilepsy_Cleaned.rds"
)

outcome_paths <- setNames(
  file.path(root, unname(outcome_files)),
  names(outcome_files)
)

all_paths <- c(exposure_paths, outcome_paths)

if (any(!file.exists(all_paths))) {
  stop(
    "Missing files:\n",
    paste(all_paths[!file.exists(all_paths)],
          collapse = "\n")
  )
}

# ------------------------------------------------------------
# SELECT EXPOSURES
# ------------------------------------------------------------
# Optional: run one exposure first.
# MR_PRESSO_EXPOSURES = "MDD" runs only MDD.
# If unset, all four exposures are processed.

selected <- Sys.getenv("MR_PRESSO_EXPOSURES", unset = "")

selected_names <- if (nzchar(selected)) {
  trimws(strsplit(selected, ",", fixed = TRUE)[[1]])
} else {
  names(exposure_paths)
}

if (
  anyDuplicated(selected_names) ||
  !all(selected_names %in% names(exposure_paths))
) {
  stop("Invalid MR_PRESSO_EXPOSURES selection.")
}

exposure_paths <- exposure_paths[selected_names]
exposures <- lapply(exposure_paths, readRDS)

instrument_snps <- unique(
  unlist(
    lapply(exposures, function(x) as.character(x$SNP)),
    use.names = FALSE
  )
)

# ------------------------------------------------------------
# MR-PRESSO SETTINGS
# ------------------------------------------------------------

nb_distribution <- 10000L
significance_threshold <- 0.05

# Fixed seed determined by each exposure-outcome combination
all_exposure_names <- c(
  "Educational_attainment", "T2D", "Stroke", "MDD"
)

seed_for_pair <- function(exposure_name, outcome_name) {
  e <- match(exposure_name, all_exposure_names)
  o <- match(outcome_name, names(outcome_paths))
  as.integer(20261029L + 10L * (o - 1L) + e)
}

# ------------------------------------------------------------
# RUN ANALYSES
# ------------------------------------------------------------

for (outcome_name in names(outcome_paths)) {

  cat("\n=== OUTCOME:", outcome_name, "===\n")

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

  for (exposure_name in names(exposures)) {

    raw_file <- file.path(
      raw_dir,
      paste0(exposure_name, "__", outcome_name, ".rds")
    )

    seed <- seed_for_pair(exposure_name, outcome_name)

    # Reuse completed analyses when the stored settings match.
    if (file.exists(raw_file)) {
      existing <- readRDS(raw_file)

      if (
        !identical(existing$seed, seed) ||
        !identical(existing$nb_distribution, nb_distribution)
      ) {
        stop("Cached analysis uses different settings: ", raw_file)
      }

      cat(
        "Already completed:",
        exposure_name, "->", outcome_name, "\n"
      )

      next
    }

    exposure <- exposures[[exposure_name]]

    matched <- outcome_subset[
      outcome_subset$SNP %in% exposure$SNP,
      ,
      drop = FALSE
    ]

    if (nrow(matched) == 0L) {
      stop("No matching SNPs: ", exposure_name, " -> ", outcome_name)
    }

    dat <- TwoSampleMR::harmonise_data(
      exposure_dat = exposure,
      outcome_dat = matched,
      action = 2
    )

    dat <- dat[
      !is.na(dat$mr_keep) & dat$mr_keep,
      ,
      drop = FALSE
    ]

    nsnp <- length(unique(dat$SNP))

    if (nsnp < 4L || nsnp != nrow(dat)) {
      stop(
        "Invalid number of instruments: ",
        exposure_name, " -> ", outcome_name
      )
    }

    if (
      any(!is.finite(dat$beta.exposure)) ||
      any(!is.finite(dat$se.exposure)) ||
      any(!is.finite(dat$beta.outcome)) ||
      any(!is.finite(dat$se.outcome)) ||
      any(dat$se.exposure <= 0) ||
      any(dat$se.outcome <= 0)
    ) {
      stop("Invalid effect or SE values: ",
           exposure_name, " -> ", outcome_name)
    }

    # Preserve SNP identities in outlier-test row names
    rownames(dat) <- dat$SNP

    cat(
      "\nRunning:", exposure_name, "->", outcome_name,
      "| SNPs:", nsnp, "\n"
    )

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
      stop("Missing MR-PRESSO global test result.")
    }

    # Preserve the full MR-PRESSO object.
    # It includes outlier, distortion and corrected MR results
    # when these are available.
    saveRDS(
      list(
        exposure = exposure_name,
        outcome = outcome_name,
        nsnp = nsnp,
        seed = seed,
        nb_distribution = nb_distribution,
        significance_threshold = significance_threshold,
        result = presso
      ),
      raw_file
    )

    cat(
      "Global RSS:", global$RSSobs,
      "| Global p:", as.character(global$Pvalue),
      "| Saved\n"
    )

    rm(dat, matched, presso)
  }

  rm(outcome_subset)
  invisible(gc())
}

# ------------------------------------------------------------
# COMBINE COMPLETED GLOBAL TESTS
# ------------------------------------------------------------

raw_files <- list.files(
  raw_dir,
  pattern = "\\.rds$",
  full.names = TRUE
)

if (length(raw_files) == 0L) {
  stop("No completed MR-PRESSO analyses found.")
}

global_results <- do.call(
  rbind,
  lapply(raw_files, function(path) {

    saved <- readRDS(path)

    global <- saved$result[
      ["MR-PRESSO results"]
    ][["Global Test"]]

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

if (anyDuplicated(
  paste(global_results$exposure_label,
        global_results$outcome_label)
)) {
  stop("Duplicate MR-PRESSO results detected.")
}

write.csv(
  global_results,
  file.path(output_dir, "08_MR_PRESSO_global_tests.csv"),
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "09_MR_PRESSO_sessionInfo.txt")
)

cat("\n====================================\n")
cat("MR-PRESSO RUN COMPLETE\n")
cat("====================================\n")

cat(
  "Completed comparisons stored:",
  nrow(global_results), "of 40\n"
)

cat("Results directory:\n", output_dir, "\n")

print(
  global_results[
    ,
    c(
      "exposure_label",
      "outcome_label",
      "nsnp",
      "rssobs",
      "global_pval"
    )
  ],
  row.names = FALSE
)
