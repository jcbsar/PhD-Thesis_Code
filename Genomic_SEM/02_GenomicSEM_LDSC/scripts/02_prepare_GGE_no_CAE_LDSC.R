
#!/usr/bin/env Rscript

# ============================================================
# Chapter 3 — Multivariable LDSC for Genomic SEM
# 02. Prepare GGE LDSC covariance object excluding CAE
#
# Purpose:
# Independently estimate the nine-phenotype LDSC covariance
# structure used in the GGE CAE-exclusion sensitivity GWAS.
#
# Phenotypes:
#   GGE: GTCS, JAE, JME
#   Brain: accumbens, amygdala, caudate,
#          pallidum, putamen, thalamus
#
# Usage:
# Rscript 02_prepare_GGE_no_CAE_LDSC.R ROOT_GSEM SOURCE_ROOT
#
# Optional:
# Rscript 02_prepare_GGE_no_CAE_LDSC.R ROOT_GSEM SOURCE_ROOT --check-only
#
# The script does not modify or re-harmonise SNP statistics.
# CAE SNP columns are excluded in 02_run_userGWAS.R.
# ============================================================

suppressPackageStartupMessages({
  library(GenomicSEM)
})

# ------------------------------------------------------------
# 1. Parse arguments
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

check_only <- length(args) == 3L &&
  identical(args[[3]], "--check-only")

if (
  !length(args) %in% c(2L, 3L) ||
  (length(args) == 3L && !check_only)
) {
  stop(
    "Usage: Rscript 02_prepare_GGE_no_CAE_LDSC.R ",
    "ROOT_GSEM SOURCE_ROOT [--check-only]"
  )
}

root_gsem <- normalizePath(
  args[[1]],
  mustWork = TRUE
)

source_root <- normalizePath(
  args[[2]],
  mustWork = TRUE
)

input_dir <- file.path(root_gsem, "input")
munged_dir <- file.path(source_root, "munged")
ld_dir <- file.path(
  source_root, "ldsc_refdata", "eur_w_ld_chr"
)

if (!dir.exists(input_dir)) {
  stop("Missing output directory: ", input_dir)
}

# ------------------------------------------------------------
# 2. Define the nine phenotypes and source files
# ------------------------------------------------------------

trait_names <- c(
  "GTCS",
  "JAE",
  "JME",
  "accumbens",
  "amygdala",
  "caudate",
  "pallidum",
  "putamen",
  "thalamus"
)

trait_filenames <- c(
  "GTCS_cau_21_Jan_25.hm3.sumstats.gz",
  "JAE_cau_21_Jan_25.hm3.sumstats.gz",
  "JME_cau_21_Jan_25.hm3.sumstats.gz",
  "accumbens_21_Jan_25.hm3.sumstats.gz",
  "amygdala_21_Jan_25.hm3.sumstats.gz",
  "caudate_21_Jan_25.hm3.sumstats.gz",
  "pallidum_21_Jan_25.hm3.sumstats.gz",
  "putamen_21_Jan_25.hm3.sumstats.gz",
  "thalamus_21_Jan_25.hm3.sumstats.gz"
)

trait_files <- file.path(
  munged_dir,
  trait_filenames
)

# Preserve the prevalence assumptions from the
# original CAE-exclusion LDSC analysis.

sample_prev <- c(
  0.50, 0.50, 0.50,
  NA, NA, NA, NA, NA, NA
)

population_prev <- c(
  0.00020, 0.00015, 0.00015,
  NA, NA, NA, NA, NA, NA
)

stopifnot(
  length(trait_names) == 9L,
  length(trait_files) == 9L,
  length(sample_prev) == 9L,
  length(population_prev) == 9L
)

# ------------------------------------------------------------
# 3. Validate external LDSC inputs
# ------------------------------------------------------------

missing_traits <- trait_files[
  !file.exists(trait_files)
]

if (length(missing_traits) > 0L) {
  stop(
    "Missing munged summary statistics:\n",
    paste(missing_traits, collapse = "\n")
  )
}

if (!dir.exists(ld_dir)) {
  stop("Missing European LD-score directory: ", ld_dir)
}

output_file <- file.path(
  input_dir,
  "ldsc_GGE_no_CAE_brain_6region.rds"
)

session_file <- file.path(
  input_dir,
  "ldsc_GGE_no_CAE_brain_6region_session_info.txt"
)

cat("\n========================================\n")
cat("GGE CAE-exclusion LDSC preparation\n")
cat("========================================\n")

cat("Phenotypes: ", length(trait_names), "\n", sep = "")
cat("GGE indicators: GTCS, JAE, JME\n")
cat("Brain indicators: 6 subcortical regions\n")
cat("LD reference: ", ld_dir, "\n", sep = "")
cat("Output: ", output_file, "\n", sep = "")

if (check_only) {
  cat("\nPASS: Required input files and directories found.\n")
  cat("LDSC was not executed.\n")
  quit(save = "no", status = 0L)
}

if (file.exists(output_file) ||
    file.exists(session_file)) {
  stop(
    "Output already exists. Refusing to overwrite ",
    "the original LDSC analysis."
  )
}

# ------------------------------------------------------------
# 4. Estimate the nine-phenotype LDSC structure
# ------------------------------------------------------------

cat("\nRunning nine-phenotype LDSC...\n")

ldsc_result <- GenomicSEM::ldsc(
  traits = trait_files,
  sample.prev = sample_prev,
  population.prev = population_prev,
  ld = paste0(ld_dir, "/"),
  wld = paste0(ld_dir, "/"),
  trait.names = trait_names
)

# ------------------------------------------------------------
# 5. Validate the resulting covariance objects
# ------------------------------------------------------------

if (
  is.null(ldsc_result$S) ||
  is.null(ldsc_result$V) ||
  is.null(ldsc_result$I)
) {
  stop("LDSC output is missing S, V or I.")
}

if (!identical(dim(ldsc_result$S), c(9L, 9L))) {
  stop("Unexpected S matrix dimensions.")
}

if (!identical(dim(ldsc_result$I), c(9L, 9L))) {
  stop("Unexpected I matrix dimensions.")
}

if (!identical(dim(ldsc_result$V), c(45L, 45L))) {
  stop("Unexpected V matrix dimensions.")
}

if (
  any(!is.finite(ldsc_result$S)) ||
  any(!is.finite(ldsc_result$V)) ||
  any(!is.finite(ldsc_result$I))
) {
  stop("Non-finite values detected in LDSC matrices.")
}

cat("\nPASS: LDSC covariance object validated.\n")
cat("S matrix: 9 x 9\n")
cat("V matrix: 45 x 45\n")
cat("I matrix: 9 x 9\n")

# ------------------------------------------------------------
# 6. Save output and analysis provenance
# ------------------------------------------------------------

temporary_output <- tempfile(
  pattern = ".ldsc_no_CAE_",
  tmpdir = input_dir,
  fileext = ".rds"
)

temporary_session <- tempfile(
  pattern = ".ldsc_no_CAE_session_",
  tmpdir = input_dir,
  fileext = ".txt"
)

saveRDS(
  ldsc_result,
  temporary_output,
  compress = "gzip"
)

provenance <- c(
  "Chapter 3: GGE CAE-exclusion LDSC analysis",
  "",
  "Phenotypes:",
  paste(trait_names, collapse = ", "),
  "",
  "Munged summary statistics:",
  trait_files,
  "",
  "Sample prevalences:",
  paste(sample_prev, collapse = ", "),
  "",
  "Population prevalences:",
  paste(population_prev, collapse = ", "),
  "",
  paste("LD reference:", ld_dir),
  paste("Weighted LD reference:", ld_dir),
  "",
  "R session information:",
  capture.output(sessionInfo())
)

writeLines(provenance, temporary_session)

if (file.exists(output_file) ||
    file.exists(session_file)) {
  stop("An output file appeared during processing.")
}

if (!file.rename(temporary_output, output_file)) {
  stop("Failed to save LDSC covariance object.")
}

if (!file.rename(temporary_session, session_file)) {
  stop("Failed to save LDSC session information.")
}

cat("\n========================================\n")
cat("PASS: CAE-EXCLUSION LDSC COMPLETED\n")
cat("========================================\n")

cat("LDSC object: ", output_file, "\n", sep = "")
cat("Session record: ", session_file, "\n", sep = "")
