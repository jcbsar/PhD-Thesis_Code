
#!/usr/bin/env Rscript

# ============================================================
# Chapter 3 — Multivariate GWAS
# 01. Validate and partition harmonised SNP summary statistics
#
# Inputs:
#   GGE_brain_sumstats.rds
#   FE_brain_sumstats.rds
#
# Outputs:
#   GGE_brain/chr01_sumstats.rds ... chr22_sumstats.rds
#   FE_brain/chr01_sumstats.rds  ... chr22_sumstats.rds
#   Chromosome-level summary and session information
#
# Usage:
# Rscript 01_prepare_inputs.R <GGE_input> <FE_input> <output_dir>
#
# Optional validation-only mode:
# Rscript 01_prepare_inputs.R <GGE_input> <FE_input> <output_dir> --check-only
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

if (
  !length(args) %in% c(3L, 4L) ||
  (length(args) == 4L && args[[4]] != "--check-only")
) {
  stop(
    "Usage: Rscript 01_prepare_inputs.R ",
    "<GGE_input.rds> <FE_input.rds> <output_dir> ",
    "[--check-only]"
  )
}

check_only <- length(args) == 4L

input_files <- c(
  GGE = args[[1]],
  FE = args[[2]]
)

output_root <- args[[3]]

expected_rows <- c(
  GGE = 4660537L,
  FE = 4178950L
)

brain_traits <- c(
  "accumbens",
  "amygdala",
  "caudate",
  "pallidum",
  "putamen",
  "thalamus"
)

analysis_traits <- list(
  GGE = c("CAE", "GTCS", "JAE", "JME", brain_traits),
  FE = c("F_HS", "F_Neg", "F_Oth", brain_traits)
)

# ------------------------------------------------------------
# 1. Validate source files
# ------------------------------------------------------------

if (!all(file.exists(input_files))) {
  stop(
    "Missing input files:\n",
    paste(
      input_files[!file.exists(input_files)],
      collapse = "\n"
    )
  )
}

input_files <- normalizePath(
  input_files,
  winslash = "/",
  mustWork = TRUE
)

# ------------------------------------------------------------
# 2. Validate an individual harmonised SNP dataset
# ------------------------------------------------------------

validate_input <- function(dat, analysis) {

  required_columns <- c(
    "SNP", "CHR", "BP", "MAF", "A1", "A2"
  )

  trait_columns <- unlist(
    lapply(
      analysis_traits[[analysis]],
      function(trait) {
        c(
          paste0("beta.", trait),
          paste0("se.", trait)
        )
      }
    ),
    use.names = FALSE
  )

  required_columns <- c(
    required_columns,
    trait_columns
  )

  missing_columns <- setdiff(
    required_columns,
    names(dat)
  )

  if (length(missing_columns) > 0L) {
    stop(
      analysis,
      " is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  if (nrow(dat) != expected_rows[[analysis]]) {
    stop(
      analysis,
      ": expected ",
      expected_rows[[analysis]],
      " SNPs; observed ",
      nrow(dat)
    )
  }

  chromosome_values <- suppressWarnings(
    as.integer(as.character(dat$CHR))
  )

  if (
    anyNA(chromosome_values) ||
    any(!chromosome_values %in% 1:22)
  ) {
    stop(
      analysis,
      ": invalid or non-autosomal chromosome values."
    )
  }

  set(dat, j = "CHR", value = chromosome_values)

  if (
    anyNA(dat$SNP) ||
    any(dat$SNP == "")
  ) {
    stop(analysis, ": missing SNP identifiers.")
  }

  if (uniqueN(dat$SNP) != nrow(dat)) {
    stop(analysis, ": duplicated SNP identifiers.")
  }

  if (
    anyNA(dat$BP) ||
    anyNA(dat$A1) ||
    anyNA(dat$A2) ||
    anyNA(dat$MAF)
  ) {
    stop(analysis, ": missing SNP metadata.")
  }

  cat(
    "PASS: ", analysis, "\n",
    "  SNPs: ", format(nrow(dat), big.mark = ","), "\n",
    "  Traits: ", length(analysis_traits[[analysis]]), "\n",
    "  Chromosomes: ", uniqueN(dat$CHR), "\n",
    "  SNP columns: ", ncol(dat), "\n",
    sep = ""
  )

  invisible(TRUE)
}

# ------------------------------------------------------------
# 3. Process GGE and FE separately
# ------------------------------------------------------------

for (analysis in names(input_files)) {

  cat("\n========================================\n")
  cat("Analysis: ", analysis, "\n", sep = "")
  cat("========================================\n")

  input_file <- input_files[[analysis]]

  cat("Reading: ", input_file, "\n", sep = "")

  dat <- as.data.table(
    readRDS(input_file)
  )

  validate_input(dat, analysis)

  if (check_only) {
    rm(dat)
    invisible(gc())
    next
  }

  analysis_dir <- file.path(
    output_root,
    paste0(analysis, "_brain")
  )

  # Never overwrite existing chromosome outputs.
  if (
    dir.exists(analysis_dir) &&
    length(list.files(analysis_dir, all.files = TRUE)) > 0L
  ) {
    stop(
      "Output directory is not empty:\n",
      analysis_dir,
      "\nUse a new output directory."
    )
  }

  dir.create(
    analysis_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  chromosome_summaries <- vector(
    "list",
    22L
  )

  # ----------------------------------------------------------
  # 4. Partition by autosomal chromosome
  # ----------------------------------------------------------

  for (chr in 1:22) {

    chromosome_data <- dat[CHR == chr]

    if (nrow(chromosome_data) == 0L) {
      stop(
        analysis,
        ": no SNPs on chromosome ",
        chr
      )
    }

    output_file <- file.path(
      analysis_dir,
      sprintf("chr%02d_sumstats.rds", chr)
    )

    if (file.exists(output_file)) {
      stop("Refusing to overwrite: ", output_file)
    }

    saveRDS(
      chromosome_data,
      output_file,
      compress = "gzip"
    )

    chromosome_summaries[[chr]] <- data.table(
      analysis = analysis,
      chromosome = chr,
      rows = nrow(chromosome_data),
      unique_snps = uniqueN(chromosome_data$SNP),
      output_file = output_file,
      size_mb = round(
        file.info(output_file)$size / 1024^2,
        2
      )
    )

    cat(
      analysis,
      " chromosome ",
      sprintf("%02d", chr),
      ": ",
      format(nrow(chromosome_data), big.mark = ","),
      " SNPs\n",
      sep = ""
    )

    rm(chromosome_data)
    invisible(gc())
  }

  # ----------------------------------------------------------
  # 5. Verify and save chromosome summary
  # ----------------------------------------------------------

  summary_dt <- rbindlist(
    chromosome_summaries
  )

  if (
    sum(summary_dt$rows) !=
    expected_rows[[analysis]]
  ) {
    stop(
      analysis,
      ": chromosome totals do not match input SNP count."
    )
  }

  fwrite(
    summary_dt,
    file.path(
      analysis_dir,
      paste0(analysis, "_split_summary.csv")
    )
  )

  cat(
    "\nPASS: ",
    analysis,
    " chromosome partitioning completed.\n",
    sep = ""
  )

  rm(dat, summary_dt)
  invisible(gc())
}

# ------------------------------------------------------------
# 6. Record software environment
# ------------------------------------------------------------

if (!check_only) {

  writeLines(
    capture.output(sessionInfo()),
    file.path(
      output_root,
      "split_sumstats_session_info.txt"
    )
  )

}

cat("\n========================================\n")

if (check_only) {
  cat("INPUT VALIDATION COMPLETED\n")
} else {
  cat("CHROMOSOME PARTITIONING COMPLETED\n")
}

cat("========================================\n")
