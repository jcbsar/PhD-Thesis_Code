
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GenomicSEM)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3L) {
  stop("Usage: Rscript 02_run_userGWAS_chr.R <GGE|FE> <chromosome 1-22> <cores>")
}

analysis <- toupper(args[[1]])
chromosome <- as.integer(args[[2]])
cores <- as.integer(args[[3]])

if (!analysis %in% c("GGE", "FE")) {
  stop("Analysis must be GGE or FE.")
}

if (is.na(chromosome) || chromosome < 1L || chromosome > 22L) {
  stop("Chromosome must be an integer from 1 to 22.")
}

if (is.na(cores) || cores < 1L) {
  stop("Cores must be a positive integer.")
}

root_dir <- "/scratch/jcbsar002/gsem"
input_dir <- file.path(root_dir, "input")
split_root <- file.path(root_dir, "input_by_chr")
results_root <- file.path(root_dir, "results")

if (analysis == "GGE") {

  covstruc_file <- file.path(
    input_dir,
    "ldsc_GGE_brain_6region.rds"
  )

  snp_file <- file.path(
    split_root,
    "GGE_brain",
    sprintf("chr%02d_sumstats.rds", chromosome)
  )

  result_dir <- file.path(results_root, "GGE_brain")

  model <- '
    GGE =~ CAE + GTCS + JAE + JME
    Brain =~ accumbens + amygdala + caudate + pallidum + putamen + thalamus
    CAE ~~ 0*CAE
    GGE ~~ Brain
    GGE ~ SNP
    Brain ~ SNP
  '

  requested_paths <- c("GGE~SNP", "Brain~SNP")
  path_labels <- c("GGE", "Brain")

} else {

  covstruc_file <- file.path(
    input_dir,
    "ldsc_FE_brain_6region.rds"
  )

  snp_file <- file.path(
    split_root,
    "FE_brain",
    sprintf("chr%02d_sumstats.rds", chromosome)
  )

  result_dir <- file.path(results_root, "FE_brain")

  model <- '
    FE =~ F_HS + F_Neg + F_Oth
    Brain =~ accumbens + amygdala + caudate + pallidum + putamen + thalamus
    FE ~~ Brain
    FE ~ SNP
    Brain ~ SNP
  '

  requested_paths <- c("FE~SNP", "Brain~SNP")
  path_labels <- c("FE", "Brain")
}

dir.create(
  result_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

prefix <- sprintf("%s_chr%02d", analysis, chromosome)

done_file <- file.path(
  result_dir,
  paste0(prefix, ".done")
)

if (file.exists(done_file)) {
  cat("Completion marker already exists; skipping: ",
      done_file, "\n", sep = "")
  quit(save = "no", status = 0)
}

for (required_file in c(covstruc_file, snp_file)) {
  if (!file.exists(required_file)) {
    stop("Missing required input: ", required_file)
  }
}

cat(
  "========================================\n",
  "Analysis: ", analysis, "\n",
  "Chromosome: ", chromosome, "\n",
  "Cores: ", cores, "\n",
  "GenomicSEM version: ", as.character(packageVersion("GenomicSEM")), "\n",
  "R version: ", R.version.string, "\n",
  "========================================\n",
  sep = ""
)

covstruc <- readRDS(covstruc_file)
SNPs <- as.data.table(readRDS(snp_file))

required_snp_columns <- c(
  "SNP", "CHR", "BP", "A1", "A2", "MAF"
)

missing_snp_columns <- setdiff(
  required_snp_columns,
  names(SNPs)
)

if (length(missing_snp_columns) > 0L) {
  stop(
    "Chromosome summary statistics are missing: ",
    paste(missing_snp_columns, collapse = ", ")
  )
}

if (nrow(SNPs) == 0L) {
  stop("No SNPs were loaded.")
}

if (any(as.integer(SNPs$CHR) != chromosome)) {
  stop("The chromosome file contains rows from another chromosome.")
}

cat(
  "SNPs loaded: ",
  format(nrow(SNPs), big.mark = ","),
  "\n",
  sep = ""
)

start_time <- Sys.time()

result <- GenomicSEM::userGWAS(
  covstruc = covstruc,
  SNPs = SNPs,
  estimation = "DWLS",
  model = model,
  printwarn = TRUE,
  sub = requested_paths,
  cores = cores,
  toler = FALSE,
  SNPSE = FALSE,
  parallel = cores > 1L,
  GC = "standard",
  MPI = FALSE,
  smooth_check = TRUE,
  std.lv = FALSE,
  fix_measurement = TRUE,
  Q_SNP = TRUE
)

end_time <- Sys.time()

runtime_minutes <- as.numeric(
  difftime(end_time, start_time, units = "mins")
)

if (!is.list(result) || length(result) != length(requested_paths)) {
  stop("Unexpected userGWAS result structure.")
}

result_rds <- file.path(
  result_dir,
  paste0(prefix, "_all_paths.rds")
)

saveRDS(
  result,
  result_rds,
  compress = "gzip"
)

path_summaries <- vector(
  "list",
  length(result)
)

for (i in seq_along(result)) {

  result_dt <- as.data.table(result[[i]])

  out_file <- file.path(
    result_dir,
    paste0(prefix, "_", path_labels[[i]], ".csv.gz")
  )

  fwrite(
    result_dt,
    out_file,
    sep = ",",
    quote = FALSE,
    na = "NA",
    compress = "gzip"
  )

  error_count <- if ("error" %in% names(result_dt)) {
    sum(
      !is.na(result_dt$error) &
        !as.character(result_dt$error) %in% c("", "0", "FALSE")
    )
  } else {
    NA_integer_
  }

  max_z_smooth <- if (
    "Z_smooth" %in% names(result_dt) &&
    any(is.finite(result_dt$Z_smooth))
  ) {
    max(result_dt$Z_smooth, na.rm = TRUE)
  } else {
    NA_real_
  }

  path_summaries[[i]] <- data.table(
    analysis = analysis,
    chromosome = chromosome,
    path = requested_paths[[i]],
    label = path_labels[[i]],
    rows = nrow(result_dt),
    errors = error_count,
    maximum_Z_smooth = max_z_smooth,
    output_file = out_file
  )
}

summary_dt <- rbindlist(path_summaries)

summary_dt[, runtime_minutes := round(runtime_minutes, 2)]
summary_dt[, cores := cores]
summary_dt[, GenomicSEM_version := as.character(packageVersion("GenomicSEM"))]

fwrite(
  summary_dt,
  file.path(result_dir, paste0(prefix, "_summary.csv"))
)

writeLines(
  capture.output(sessionInfo()),
  file.path(result_dir, paste0(prefix, "_session_info.txt"))
)

writeLines(
  c(
    paste0("analysis=", analysis),
    paste0("chromosome=", chromosome),
    paste0(
      "completed=",
      format(Sys.time(), tz = "Africa/Johannesburg")
    ),
    paste0("runtime_minutes=", round(runtime_minutes, 2)),
    paste0("snps=", nrow(SNPs)),
    paste0(
      "GenomicSEM_version=",
      as.character(packageVersion("GenomicSEM"))
    )
  ),
  done_file
)

cat(
  "\nCompleted successfully.\n",
  "Runtime: ", round(runtime_minutes, 2),
  " minutes\nCompletion marker: ",
  done_file,
  "\n",
  sep = ""
)
