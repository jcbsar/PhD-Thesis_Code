
#!/usr/bin/env Rscript

# ============================================================
# Chapter 3 — Multivariate GWAS
# 04. Independently validate combined Genomic SEM results
#
# Usage:
#   Rscript 04_validate_results.R ROOT ANALYSIS
#
# ANALYSIS: GGE | FE | GGE_no_CAE
#
# Verifies:
#   - Complete SNP and chromosome coverage
#   - Exact SNP identity against harmonised inputs
#   - Both factor paths and SNP metadata
#   - Valid association P-values and no model errors
#   - Genome-wide significant SNP counts
#   - Q-SNP results and smoothing diagnostics
#
# Writes reports to ROOT/validation/ANALYSIS/.
# Does not modify GWAS result files.
# ============================================================

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2L) {
  stop("Usage: Rscript 04_validate_results.R ROOT ANALYSIS")
}

root <- normalizePath(args[[1]], mustWork = TRUE)
analysis <- args[[2]]

# ------------------------------------------------------------
# 1. Define the verified reference results
# ------------------------------------------------------------

specifications <- list(

  GGE = list(
    input = "GGE_brain",
    factors = c("GGE", "Brain"),
    total = 4660537L,
    expected_gws = c(GGE = 249L, Brain = 608L)
  ),

  FE = list(
    input = "FE_brain",
    factors = c("FE", "Brain"),
    total = 4178950L,
    expected_gws = c(FE = 0L, Brain = 594L)
  ),

  GGE_no_CAE = list(
    input = "GGE_brain",
    factors = c("GGE", "Brain"),
    total = 4660537L,
    expected_gws = c(GGE = 68L, Brain = NA_integer_)
  )
)

if (!analysis %in% names(specifications)) {
  stop("Analysis must be GGE, FE or GGE_no_CAE.")
}

cfg <- specifications[[analysis]]

result_dir <- file.path(
  root, "final_consolidated", analysis
)

split_dir <- file.path(
  root, "input_by_chr", cfg$input
)

summary_file <- file.path(
  split_dir,
  paste0(
    if (analysis == "FE") "FE" else "GGE",
    "_split_summary.csv"
  )
)

if (!dir.exists(result_dir) ||
    !file.exists(summary_file)) {
  stop("Required combined results or input summary missing.")
}

# ------------------------------------------------------------
# 2. Validate chromosome reference counts
# ------------------------------------------------------------

split_summary <- fread(summary_file)
setorder(split_summary, chromosome)

if (
  nrow(split_summary) != 22L ||
  !identical(as.integer(split_summary$chromosome), 1:22) ||
  anyNA(split_summary$rows) ||
  sum(split_summary$rows) != cfg$total
) {
  stop("Invalid input chromosome summary.")
}

# ------------------------------------------------------------
# 3. Define required output columns
# ------------------------------------------------------------

required_columns <- c(
  "SNP", "CHR", "BP", "A1", "A2",
  "lhs", "op", "rhs",
  "est", "SE", "Z_Estimate",
  "Pval_Estimate",
  "Q_SNP", "Q_SNP_pval",
  "error", "Z_smooth"
)

summary_list <- list()
reference_metadata <- NULL

cat("\n========================================\n")
cat("Validating analysis: ", analysis, "\n", sep = "")
cat("========================================\n")

# ------------------------------------------------------------
# 4. Independently validate each factor
# ------------------------------------------------------------

for (factor in cfg$factors) {

  result_file <- file.path(
    result_dir,
    paste0(
      analysis, "_", factor,
      "_userGWAS_full.csv"
    )
  )

  if (!file.exists(result_file)) {
    stop("Missing combined result: ", result_file)
  }

  header <- names(fread(result_file, nrows = 0L))

  missing_columns <- setdiff(required_columns, header)

  if (length(missing_columns) > 0L) {
    stop(
      "Missing result columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  tab <- fread(
    result_file,
    select = required_columns,
    showProgress = TRUE
  )

  # SNP count and uniqueness
  if (
    nrow(tab) != cfg$total ||
    anyNA(tab$SNP) ||
    any(as.character(tab$SNP) == "") ||
    uniqueN(tab$SNP) != cfg$total
  ) {
    stop(factor, ": SNP count or uniqueness check failed.")
  }

  # Chromosome and coordinate validity
  if (
    anyNA(tab$CHR) ||
    any(!tab$CHR %in% 1:22) ||
    anyNA(tab$BP) ||
    any(tab$BP < 1)
  ) {
    stop(factor, ": invalid chromosome coordinates.")
  }

  chromosome_counts <- tabulate(
    as.integer(tab$CHR),
    nbins = 22L
  )

  if (!identical(
    as.integer(chromosome_counts),
    as.integer(split_summary$rows)
  )) {
    stop(factor, ": chromosome SNP counts do not match inputs.")
  }

  # Factor-path identity
  if (
    anyNA(tab$lhs) ||
    anyNA(tab$op) ||
    anyNA(tab$rhs) ||
    any(as.character(tab$lhs) != factor) ||
    any(as.character(tab$op) != "~") ||
    any(as.character(tab$rhs) != "SNP")
  ) {
    stop(factor, ": incorrect factor-path labels.")
  }

  # Verify alleles are present
  if (
    anyNA(tab$A1) ||
    anyNA(tab$A2) ||
    any(as.character(tab$A1) == "") ||
    any(as.character(tab$A2) == "")
  ) {
    stop(factor, ": missing allele information.")
  }

  # Verify that the two factors refer to identical variants
  # with identical chromosome, position and allele metadata.
  metadata <- tab[
    ,
    .(
      SNP = as.character(SNP),
      CHR = as.integer(CHR),
      BP = as.integer(BP),
      A1 = as.character(A1),
      A2 = as.character(A2)
    )
  ]

  if (is.null(reference_metadata)) {

    reference_metadata <- metadata

  } else if (!identical(metadata, reference_metadata)) {

    stop(
      factor,
      ": variant identity/order or SNP metadata ",
      "differs between the two factor results."
    )
  }

  # Association P-values
  if (
    !is.numeric(tab$Pval_Estimate) ||
    anyNA(tab$Pval_Estimate) ||
    any(!is.finite(tab$Pval_Estimate)) ||
    any(tab$Pval_Estimate < 0 |
        tab$Pval_Estimate > 1)
  ) {
    stop(factor, ": invalid association P-values.")
  }

  # Model errors
  errors <- trimws(as.character(tab$error))

  has_error <- !is.na(errors) &
    !(errors %in% c("", "0"))

  if (any(has_error)) {
    stop(
      factor, ": ",
      sum(has_error),
      " SNPs have model errors."
    )
  }

  # Q-SNP heterogeneity P-values
  qsnp_p <- suppressWarnings(
    as.numeric(tab$Q_SNP_pval)
  )

  if (any(
    !is.na(qsnp_p) &
      (!is.finite(qsnp_p) |
         qsnp_p < 0 |
         qsnp_p > 1)
  )) {
    stop(factor, ": invalid Q-SNP P-values.")
  }

  # Genome-wide significance
  observed_gws <- sum(
    tab$Pval_Estimate < 5e-8
  )

  expected_gws <- cfg$expected_gws[[factor]]

  if (
    !is.na(expected_gws) &&
    observed_gws != expected_gws
  ) {
    stop(
      factor,
      ": expected ", expected_gws,
      " genome-wide significant SNPs; observed ",
      observed_gws,
      ". Investigate before reporting replication."
    )
  }

  # Smoothing diagnostic
  smooth <- suppressWarnings(
    as.numeric(tab$Z_smooth)
  )

  max_smooth <- if (all(is.na(smooth))) {
    NA_real_
  } else {
    max(abs(smooth), na.rm = TRUE)
  }

  summary_list[[factor]] <- data.table(
    analysis = analysis,
    factor = factor,
    n_snps = nrow(tab),
    n_chromosomes = uniqueN(tab$CHR),
    n_gws_snps = observed_gws,
    expected_gws_snps = expected_gws,
    n_model_errors = sum(has_error),
    n_qsnp_tested = sum(!is.na(qsnp_p)),
    n_qsnp_nominal_p_lt_0_05 = sum(
      qsnp_p < 0.05, na.rm = TRUE
    ),
    max_abs_Z_smooth = max_smooth,
    result_md5 = unname(tools::md5sum(result_file))
  )

  cat(
    "PASS: ", factor,
    " — ", format(nrow(tab), big.mark = ","),
    " SNPs; ", observed_gws,
    " genome-wide significant\n",
    sep = ""
  )

  rm(tab, metadata)
  invisible(gc(verbose = FALSE))
}

# ------------------------------------------------------------
# 5. Compare variant identities with original inputs
# ------------------------------------------------------------

cat("\nChecking SNP identities against chromosome inputs...\n")

for (chr in 1:22) {

  input_file <- file.path(
    split_dir,
    sprintf("chr%02d_sumstats.rds", chr)
  )

  if (!file.exists(input_file)) {
    stop("Missing original chromosome input: ", input_file)
  }

  input_data <- readRDS(input_file)

  input_snps <- as.character(input_data$SNP)

  result_snps <- reference_metadata[
    CHR == chr,
    SNP
  ]

  if (
    length(input_snps) != length(result_snps) ||
    !setequal(input_snps, result_snps)
  ) {
    stop(
      "Chromosome ", chr,
      ": result SNP identities differ from inputs."
    )
  }

  rm(input_data, input_snps, result_snps)
  invisible(gc(verbose = FALSE))
}

cat("PASS: SNP identities match the harmonised inputs.\n")

# ------------------------------------------------------------
# 6. Save the independent validation report
# ------------------------------------------------------------

validation_root <- file.path(root, "validation")
validation_dir <- file.path(validation_root, analysis)

if (file.exists(validation_dir) ||
    dir.exists(validation_dir)) {
  stop(
    "Validation directory already exists: ",
    validation_dir,
    "\nRefusing to overwrite."
  )
}

dir.create(
  validation_root,
  recursive = TRUE,
  showWarnings = FALSE
)

stage_dir <- tempfile(
  pattern = paste0(".", analysis, "_validation_"),
  tmpdir = validation_root
)

if (!dir.create(stage_dir)) {
  stop("Unable to create validation staging directory.")
}

validation_summary <- rbindlist(summary_list)

fwrite(
  validation_summary,
  file.path(stage_dir, "validation_summary.csv")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(stage_dir, "session_info.txt")
)

if (file.exists(validation_dir) ||
    dir.exists(validation_dir) ||
    !file.rename(stage_dir, validation_dir)) {
  stop("Unable to publish validation report.")
}

cat("\n========================================\n")
cat("PASS: INDEPENDENT GWAS VALIDATION\n")
cat("========================================\n")

print(validation_summary)

cat(
  "\nReport directory: ",
  validation_dir,
  "\n",
  sep = ""
)
