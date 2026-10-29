
#!/usr/bin/env Rscript

# ============================================================
# Chapter 3 — Multivariate GWAS
# 02. Run Genomic SEM userGWAS by chromosome
#
# Analyses:
#   GGE         Corrected primary GGE + brain model
#   FE          Primary FE + brain model
#   GGE_no_CAE  GGE sensitivity model excluding CAE
#
# Usage:
# Rscript 02_run_userGWAS.R ROOT ANALYSIS CHR [CORES]
# Rscript 02_run_userGWAS.R ROOT ANALYSIS CHR CORES PART PARTS
#
# Append --check-only to validate inputs without running GWAS.
#
# Each job writes one RDS containing the two factor results.
# Results are kept separate from historical HPC outputs.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicSEM)
})

# ------------------------------------------------------------
# 1. Command-line arguments
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

check_only <- length(args) > 0L &&
  identical(tail(args, 1L), "--check-only")

if (check_only) {
  args <- head(args, -1L)
}

if (!(length(args) %in% c(3L, 4L, 6L))) {
  stop(
    "Usage: Rscript 02_run_userGWAS.R ",
    "ROOT ANALYSIS CHR [CORES] [PART PARTS] ",
    "[--check-only]"
  )
}

root <- normalizePath(args[[1]], mustWork = TRUE)
analysis <- args[[2]]
chr <- suppressWarnings(as.integer(args[[3]]))

if (!analysis %in% c("GGE", "FE", "GGE_no_CAE")) {
  stop("Analysis must be GGE, FE or GGE_no_CAE.")
}

if (is.na(chr) || chr < 1L || chr > 22L) {
  stop("Chromosome must be an integer from 1 to 22.")
}

cores <- if (length(args) >= 4L) {
  suppressWarnings(as.integer(args[[4]]))
} else if (analysis == "GGE_no_CAE") {
  1L
} else {
  4L
}

part <- if (length(args) == 6L) {
  suppressWarnings(as.integer(args[[5]]))
} else {
  1L
}

parts <- if (length(args) == 6L) {
  suppressWarnings(as.integer(args[[6]]))
} else {
  1L
}

if (anyNA(c(cores, part, parts)) ||
    cores < 1L || parts < 1L ||
    part < 1L || part > parts) {
  stop("Invalid CPU or chunk arguments.")
}

# Preserve the serial approach used in the CAE sensitivity
# analysis; primary models can use multiple CPUs.
if (analysis == "GGE_no_CAE" && cores != 1L) {
  stop("GGE_no_CAE must run with one CPU.")
}

# ------------------------------------------------------------
# 2. Define the original models
# ------------------------------------------------------------

brain_traits <- c(
  "accumbens",
  "amygdala",
  "caudate",
  "pallidum",
  "putamen",
  "thalamus"
)

brain_model <- paste(
  "Brain =~ accumbens + amygdala + caudate +",
  "pallidum + putamen + thalamus"
)

gge_model <- paste(
  "GGE =~ CAE + GTCS + JAE + JME",
  brain_model,
  "CAE ~~ 0*CAE",
  "GGE ~~ Brain",
  "GGE ~ SNP",
  "Brain ~ SNP",
  sep = "\n"
)

fe_model <- paste(
  "FE =~ F_HS + F_Neg + F_Oth",
  brain_model,
  "FE ~~ Brain",
  "FE ~ SNP",
  "Brain ~ SNP",
  sep = "\n"
)

no_cae_model <- paste(
  "GGE =~ GTCS + JAE + JME",
  brain_model,
  "GGE ~~ Brain",
  "GGE ~ SNP",
  "Brain ~ SNP",
  sep = "\n"
)

# ------------------------------------------------------------
# 3. Analysis-specific configuration
# ------------------------------------------------------------

config <- switch(
  analysis,

  GGE = list(
    input = "GGE_brain",
    ldsc = "ldsc_GGE_brain_6region_corrected.rds",
    traits = c(
      "CAE", "GTCS", "JAE", "JME", brain_traits
    ),
    model = gge_model,
    paths = c("GGE~SNP", "Brain~SNP"),
    labels = c("GGE", "Brain")
  ),

  FE = list(
    input = "FE_brain",
    ldsc = "ldsc_FE_brain_6region.rds",
    traits = c(
      "F_HS", "F_Neg", "F_Oth", brain_traits
    ),
    model = fe_model,
    paths = c("FE~SNP", "Brain~SNP"),
    labels = c("FE", "Brain")
  ),

  GGE_no_CAE = list(
    input = "GGE_brain",
    ldsc = "ldsc_GGE_no_CAE_brain_6region.rds",
    traits = c(
      "GTCS", "JAE", "JME", brain_traits
    ),
    model = no_cae_model,
    paths = c("GGE~SNP", "Brain~SNP"),
    labels = c("GGE", "Brain")
  )
)

# ------------------------------------------------------------
# 4. Load and validate inputs
# ------------------------------------------------------------

ldsc_file <- file.path(
  root, "input", config$ldsc
)

snp_file <- file.path(
  root,
  "input_by_chr",
  config$input,
  sprintf("chr%02d_sumstats.rds", chr)
)

if (!file.exists(ldsc_file)) {
  stop("Missing LDSC input: ", ldsc_file)
}

if (!file.exists(snp_file)) {
  stop("Missing SNP input: ", snp_file)
}

covstruc <- readRDS(ldsc_file)

k <- length(config$traits)

if (
  is.null(covstruc$S) ||
  is.null(covstruc$V) ||
  is.null(covstruc$I) ||
  !identical(dim(covstruc$S), c(k, k)) ||
  !identical(dim(covstruc$I), c(k, k)) ||
  !identical(dim(covstruc$V), c(k * (k + 1L) / 2L,
                                k * (k + 1L) / 2L))
) {
  stop("LDSC covariance object has unexpected dimensions.")
}

snps <- as.data.table(readRDS(snp_file))

if (analysis == "GGE_no_CAE") {
  cae_columns <- c("beta.CAE", "se.CAE")

  if (!all(cae_columns %in% names(snps))) {
    stop("CAE columns are missing from GGE input.")
  }

  snps[, (cae_columns) := NULL]
}

required <- c(
  "SNP", "CHR", "BP", "A1", "A2", "MAF",
  unlist(
    lapply(
      config$traits,
      function(x) c(paste0("beta.", x),
                    paste0("se.", x))
    ),
    use.names = FALSE
  )
)

missing <- setdiff(required, names(snps))

if (length(missing) > 0L) {
  stop(
    "Missing SNP columns: ",
    paste(missing, collapse = ", ")
  )
}

if (
  nrow(snps) == 0L ||
  anyNA(snps$CHR) ||
  any(snps$CHR != chr) ||
  anyNA(snps$SNP) ||
  uniqueN(snps$SNP) != nrow(snps)
) {
  stop("Invalid chromosome SNP data.")
}

# ------------------------------------------------------------
# 5. Select the requested chunk
# ------------------------------------------------------------

total_chr_snps <- nrow(snps)

if (parts > total_chr_snps) {
  stop("More chunks requested than available SNPs.")
}

first <- floor((part - 1L) * total_chr_snps / parts) + 1L
last <- floor(part * total_chr_snps / parts)

snps <- snps[seq.int(first, last)]

cat("\nAnalysis: ", analysis, "\n", sep = "")
cat("Chromosome: ", chr, "\n", sep = "")
cat("Chunk: ", part, "/", parts, "\n", sep = "")
cat("Chunk SNPs: ", nrow(snps), "\n", sep = "")
cat("CPUs: ", cores, "\n", sep = "")
cat("LDSC input: ", ldsc_file, "\n", sep = "")

if (check_only) {
  cat("\nPASS: Inputs validated; GWAS not executed.\n")
  quit(save = "no", status = 0L)
}

# ------------------------------------------------------------
# 6. Protect existing results
# ------------------------------------------------------------

output_dir <- file.path(
  root, "results", "consolidated", analysis
)

output_file <- file.path(
  output_dir,
  sprintf(
    "chr%02d_part%02d_of%02d_paths.rds",
    chr, part, parts
  )
)

if (file.exists(output_file)) {
  stop("Refusing to overwrite existing result: ", output_file)
}

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 7. Run the multivariate GWAS
# ------------------------------------------------------------

started <- Sys.time()

results <- GenomicSEM::userGWAS(
  covstruc = covstruc,
  SNPs = snps,
  estimation = "DWLS",
  model = config$model,
  printwarn = TRUE,
  sub = config$paths,
  cores = cores,
  toler = FALSE,
  SNPSE = FALSE,
  parallel = cores > 1L,
  GC = "standard",
  MPI = FALSE,
  smooth_check = TRUE,
  TWAS = FALSE,
  std.lv = FALSE,
  fix_measurement = TRUE,
  Q_SNP = TRUE
)

if (!is.list(results) ||
    length(results) != length(config$paths)) {
  stop("Unexpected userGWAS result structure.")
}

# ------------------------------------------------------------
# 8. Validate both factor outputs
# ------------------------------------------------------------

names(results) <- config$labels

for (label in config$labels) {

  tab <- as.data.table(results[[label]])

  if (!all(c("SNP", "CHR", "BP") %in% names(tab))) {
    stop(label, ": missing essential result columns.")
  }

  if (nrow(tab) != nrow(snps)) {
    stop(label, ": result row count mismatch.")
  }

  if (
    anyNA(tab$SNP) ||
    uniqueN(tab$SNP) != nrow(tab) ||
    !setequal(
      as.character(tab$SNP),
      as.character(snps$SNP)
    )
  ) {
    stop(label, ": SNP identifiers do not match input.")
  }

  if ("error" %in% names(tab)) {

  errors <- trimws(as.character(tab$error))

  has_error <- !is.na(errors) &
    !(errors %in% c("", "0"))

  if (any(has_error)) {
    stop(
      label,
      ": ",
      sum(has_error),
      " SNPs have model errors."
    )
  }
}

  results[[label]] <- tab
}

# ------------------------------------------------------------
# 9. Save results and provenance
# ------------------------------------------------------------

output <- list(
  analysis = analysis,
  chromosome = chr,
  part = part,
  parts = parts,
  chromosome_total_snps = total_chr_snps,
  chunk_snps = nrow(snps),
  ldsc_file = ldsc_file,
  snp_file = snp_file,
  traits = config$traits,
  model = config$model,
  paths = config$paths,
  settings = list(
    estimation = "DWLS",
    cores = cores,
    parallel = cores > 1L,
    GC = "standard",
    smooth_check = TRUE,
    fix_measurement = TRUE,
    Q_SNP = TRUE,
    toler = FALSE,
    SNPSE = FALSE,
    MPI = FALSE,
    std.lv = FALSE
  ),
  results = results,
  runtime_minutes = as.numeric(
    difftime(Sys.time(), started, units = "mins")
  ),
  session_info = capture.output(sessionInfo())
)

temporary_file <- tempfile(
  pattern = "userGWAS_",
  tmpdir = output_dir,
  fileext = ".rds"
)

saveRDS(
  output,
  temporary_file,
  compress = "gzip"
)

if (file.exists(output_file)) {
  unlink(temporary_file)
  stop("Output appeared during execution; refusing overwrite.")
}

if (!file.rename(temporary_file, output_file)) {
  unlink(temporary_file)
  stop("Failed to finalise output file.")
}

cat("\nPASS: GWAS completed and validated.\n")
cat("Saved: ", output_file, "\n", sep = "")
