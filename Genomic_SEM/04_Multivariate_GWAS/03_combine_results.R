
#!/usr/bin/env Rscript

# ============================================================
# Chapter 3 — Multivariate GWAS
# 03. Combine chromosome-level Genomic SEM results
#
# Usage:
#   Rscript 03_combine_results.R ROOT ANALYSIS
#
# ANALYSIS: GGE | FE | GGE_no_CAE
#
# Reads results produced by 02_run_userGWAS.R.
# Checks chromosome and chunk completeness, SNP identities,
# factor labels, row counts and model errors.
#
# Outputs:
#   final_consolidated/ANALYSIS/
#     ANALYSIS_FACTOR_userGWAS_full.csv
#     ANALYSIS_combination_summary.csv
#     ANALYSIS_factor_summary.csv
#     session_info.txt
#
# Never overwrites completed output directories.
# ============================================================

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2L) {
  stop("Usage: Rscript 03_combine_results.R ROOT ANALYSIS")
}

root <- normalizePath(args[[1]], mustWork = TRUE)
analysis <- args[[2]]

if (!analysis %in% c("GGE", "FE", "GGE_no_CAE")) {
  stop("Analysis must be GGE, FE or GGE_no_CAE.")
}

factors <- if (analysis == "FE") {
  c("FE", "Brain")
} else {
  c("GGE", "Brain")
}

input_group <- if (analysis == "FE") {
  "FE_brain"
} else {
  "GGE_brain"
}

expected_total <- if (analysis == "FE") {
  4178950L
} else {
  4660537L
}

result_dir <- file.path(
  root, "results", "consolidated", analysis
)

split_summary_file <- file.path(
  root,
  "input_by_chr",
  input_group,
  paste0(
    if (analysis == "FE") "FE" else "GGE",
    "_split_summary.csv"
  )
)

if (!dir.exists(result_dir)) {
  stop("Missing result directory: ", result_dir)
}

if (!file.exists(split_summary_file)) {
  stop("Missing chromosome summary: ", split_summary_file)
}

# ------------------------------------------------------------
# 1. Validate original chromosome counts
# ------------------------------------------------------------

split_summary <- fread(split_summary_file)
setorder(split_summary, chromosome)

if (
  nrow(split_summary) != 22L ||
  !identical(as.integer(split_summary$chromosome), 1:22) ||
  anyNA(split_summary$rows) ||
  sum(split_summary$rows) != expected_total
) {
  stop("Chromosome summary does not match expected input.")
}

expected_chr_rows <- split_summary$rows

# ------------------------------------------------------------
# 2. Discover and validate all completed chunks
# ------------------------------------------------------------

pattern <- paste0(
  "^chr([0-9]{2})_part([0-9]{2})",
  "_of([0-9]{2})_paths\\.rds$"
)

files <- list.files(
  result_dir,
  pattern = pattern,
  full.names = TRUE
)

if (length(files) == 0L) {
  stop("No consolidated chromosome results found.")
}

file_index <- rbindlist(lapply(files, function(file) {

  match <- regmatches(
    basename(file),
    regexec(pattern, basename(file))
  )[[1]]

  data.table(
    chromosome = as.integer(match[[2]]),
    part = as.integer(match[[3]]),
    parts = as.integer(match[[4]]),
    file = file
  )
}))

setorder(file_index, chromosome, part)

if (any(!file_index$chromosome %in% 1:22)) {
  stop("Invalid chromosome number in result filenames.")
}

for (chr in 1:22) {

  x <- file_index[chromosome == chr]

  if (nrow(x) == 0L) {
    stop("Missing results for chromosome ", chr)
  }

  if (
    uniqueN(x$parts) != 1L ||
    nrow(x) != x$parts[[1]] ||
    !setequal(x$part, seq_len(x$parts[[1]]))
  ) {
    stop(
      "Missing, duplicated or inconsistent chunks ",
      "on chromosome ", chr
    )
  }
}

cat(
  "\nAnalysis: ", analysis,
  "\nChromosomes: 22",
  "\nResult chunks: ", nrow(file_index),
  "\nExpected SNPs: ", expected_total,
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 3. Define exact expected Genomic SEM output columns
# ------------------------------------------------------------

expected_columns <- c(
  "SNP", "CHR", "BP", "MAF", "A1", "A2",
  "lhs", "op", "rhs", "free", "label",
  "est", "SE", "Z_Estimate", "Pval_Estimate",
  "chisq", "chisq_df", "chisq_pval", "AIC",
  "Q_SNP", "Q_SNP_df", "Q_SNP_pval",
  "error", "warning", "Z_smooth"
)

# ------------------------------------------------------------
# 4. Create a protected staging directory
# ------------------------------------------------------------

final_root <- file.path(root, "final_consolidated")
final_dir <- file.path(final_root, analysis)

if (file.exists(final_dir) || dir.exists(final_dir)) {
  stop(
    "Final output directory already exists: ",
    final_dir,
    "\nRefusing to overwrite."
  )
}

dir.create(
  final_root,
  recursive = TRUE,
  showWarnings = FALSE
)

stage_dir <- tempfile(
  pattern = paste0(".", analysis, "_staging_"),
  tmpdir = final_root
)

if (!dir.create(stage_dir)) {
  stop("Unable to create staging directory.")
}

output_files <- setNames(
  file.path(
    stage_dir,
    paste0(
      analysis, "_", factors,
      "_userGWAS_full.csv"
    )
  ),
  factors
)

# ------------------------------------------------------------
# 5. Combine results in chromosome and chunk order
# ------------------------------------------------------------

summary_rows <- list()
snp_chunks <- list()
row_totals <- setNames(integer(length(factors)), factors)

for (chr in 1:22) {

  chr_files <- file_index[chromosome == chr]
  observed_chr_rows <- 0L

  for (j in seq_len(nrow(chr_files))) {

    entry <- chr_files[j]
    object <- readRDS(entry$file)

    if (
      !identical(object$analysis, analysis) ||
      !identical(as.integer(object$chromosome), chr) ||
      !identical(as.integer(object$part), entry$part) ||
      !identical(as.integer(object$parts), entry$parts)
    ) {
      stop("Inconsistent result metadata: ", entry$file)
    }

    if (
      !is.list(object$results) ||
      !identical(names(object$results), factors)
    ) {
      stop("Unexpected factor result structure: ", entry$file)
    }

    chunk_rows <- as.integer(object$chunk_snps)

    if (length(chunk_rows) != 1L ||
        is.na(chunk_rows) ||
        chunk_rows < 1L) {
      stop("Invalid chunk SNP count: ", entry$file)
    }

    reference_ids <- as.character(
      object$results[[factors[[1]]]]$SNP
    )

    if (
      length(reference_ids) != chunk_rows ||
      anyNA(reference_ids) ||
      any(reference_ids == "") ||
      uniqueN(reference_ids) != chunk_rows
    ) {
      stop("Invalid SNP identifiers: ", entry$file)
    }

    snp_chunks[[length(snp_chunks) + 1L]] <- reference_ids
    observed_chr_rows <- observed_chr_rows + chunk_rows

    for (factor in factors) {

      tab <- as.data.table(object$results[[factor]])

      if (
        !identical(names(tab), expected_columns) ||
        nrow(tab) != chunk_rows
      ) {
        stop(
          "Unexpected result columns or row count: ",
          entry$file, " / ", factor
        )
      }

      if (
        !identical(as.character(tab$SNP), reference_ids) ||
        anyNA(tab$CHR) ||
        any(as.integer(tab$CHR) != chr) ||
        anyNA(tab$lhs) ||
        any(as.character(tab$lhs) != factor) ||
        anyNA(tab$op) ||
        any(as.character(tab$op) != "~") ||
        anyNA(tab$rhs) ||
        any(as.character(tab$rhs) != "SNP")
      ) {
        stop(
          "SNP or factor-path mismatch: ",
          entry$file, " / ", factor
        )
      }

      errors <- trimws(as.character(tab$error))

      has_error <- !is.na(errors) &
        !(errors %in% c("", "0"))

      if (any(has_error)) {
        stop(
          factor, ": ",
          sum(has_error),
          " model errors in ", entry$file
        )
      }

      if (
        !is.numeric(tab$Pval_Estimate) ||
        anyNA(tab$Pval_Estimate) ||
        any(!is.finite(tab$Pval_Estimate)) ||
        any(tab$Pval_Estimate < 0 |
            tab$Pval_Estimate > 1)
      ) {
        stop("Invalid association P-values: ", entry$file)
      }

      append_file <- file.exists(output_files[[factor]])

      fwrite(
        tab,
        file = output_files[[factor]],
        append = append_file,
        col.names = !append_file,
        quote = "auto"
      )

      smooth <- suppressWarnings(
        as.numeric(tab$Z_smooth)
      )

      maximum_smooth <- if (all(is.na(smooth))) {
        NA_real_
      } else {
        max(abs(smooth), na.rm = TRUE)
      }

      summary_rows[[length(summary_rows) + 1L]] <- data.table(
        analysis = analysis,
        factor = factor,
        chromosome = chr,
        part = entry$part,
        parts = entry$parts,
        rows = nrow(tab),
        gws_snps = sum(tab$Pval_Estimate < 5e-8),
        maximum_abs_Z_smooth = maximum_smooth
      )

      row_totals[[factor]] <- row_totals[[factor]] + nrow(tab)

      rm(tab)
    }

    rm(object, reference_ids)
    invisible(gc(verbose = FALSE))
  }

  if (observed_chr_rows != expected_chr_rows[[chr]]) {
    stop(
      "Chromosome ", chr,
      ": expected ", expected_chr_rows[[chr]],
      " SNPs; observed ", observed_chr_rows
    )
  }

  cat(
    "PASS: chromosome ",
    sprintf("%02d", chr),
    " — ", observed_chr_rows,
    " SNPs\n",
    sep = ""
  )
}

# ------------------------------------------------------------
# 6. Validate genome-wide SNP uniqueness
# ------------------------------------------------------------

all_snps <- unlist(
  snp_chunks,
  use.names = FALSE
)

if (
  length(all_snps) != expected_total ||
  uniqueN(all_snps) != expected_total
) {
  stop(
    "Genome-wide SNP count or uniqueness check failed."
  )
}

if (any(row_totals != expected_total)) {
  stop("Factor result row totals do not match expected SNPs.")
}

# ------------------------------------------------------------
# 7. Save QC and provenance
# ------------------------------------------------------------

combination_summary <- rbindlist(summary_rows)

factor_summary <- combination_summary[
  ,
  .(
    n_snps = sum(rows),
    n_gws_snps = sum(gws_snps),
    max_abs_Z_smooth = if (all(is.na(maximum_abs_Z_smooth))) {
      NA_real_
    } else {
      max(maximum_abs_Z_smooth, na.rm = TRUE)
    }
  ),
  by = .(analysis, factor)
]

fwrite(
  combination_summary,
  file.path(
    stage_dir,
    paste0(analysis, "_combination_summary.csv")
  )
)

fwrite(
  factor_summary,
  file.path(
    stage_dir,
    paste0(analysis, "_factor_summary.csv")
  )
)

writeLines(
  capture.output(sessionInfo()),
  file.path(stage_dir, "session_info.txt")
)

# ------------------------------------------------------------
# 8. Publish completed outputs together
# ------------------------------------------------------------

if (file.exists(final_dir) || dir.exists(final_dir)) {
  stop("Final directory appeared during processing.")
}

if (!file.rename(stage_dir, final_dir)) {
  stop(
    "Failed to publish completed outputs. ",
    "Staged files remain at: ", stage_dir
  )
}

cat("\n========================================\n")
cat("PASS: GWAS RESULTS COMBINED AND VALIDATED\n")
cat("========================================\n")

print(factor_summary)

cat("\nFinal output directory:\n", final_dir, "\n")
