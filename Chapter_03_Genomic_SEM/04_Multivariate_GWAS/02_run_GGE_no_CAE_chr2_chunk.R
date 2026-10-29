
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GenomicSEM)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1L) {
  stop("Usage: Rscript 02_run_GGE_no_CAE_chr2_chunk.R <part 1|2>")
}

part <- as.integer(args[[1]])

if (is.na(part) || !part %in% c(1L, 2L)) {
  stop("Part must be 1 or 2.")
}

root_dir <- "/scratch/jcbsar002/gsem"

covstruc_file <- file.path(
  root_dir,
  "input",
  "ldsc_GGE_no_CAE_brain_6region.rds"
)

snp_file <- file.path(
  root_dir,
  "input_by_chr",
  "GGE_brain",
  "chr02_sumstats.rds"
)

result_dir <- file.path(
  root_dir,
  "results",
  "GGE_no_CAE_sensitivity",
  "chr02_parts"
)

dir.create(
  result_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(covstruc_file)) {
  stop("Missing covariance object: ", covstruc_file)
}

if (!file.exists(snp_file)) {
  stop("Missing chr2 SNP file: ", snp_file)
}

model <- '
  GGE =~ GTCS + JAE + JME

  Brain =~ accumbens + amygdala + caudate +
           pallidum + putamen + thalamus

  GGE ~~ Brain

  GGE ~ SNP
  Brain ~ SNP
'

requested_paths <- c(
  "GGE~SNP",
  "Brain~SNP"
)

path_labels <- c(
  "GGE",
  "Brain"
)

cat("========================================\n")
cat("GGE no-CAE chromosome 2 chunk analysis\n")
cat("Part: ", part, "\n", sep = "")
cat("GenomicSEM version: ",
    as.character(packageVersion("GenomicSEM")), "\n", sep = "")
cat("R version: ", R.version.string, "\n", sep = "")
cat("========================================\n")

covstruc <- readRDS(covstruc_file)

SNPs <- as.data.table(
  readRDS(snp_file)
)

# Remove CAE effects to match the 9-trait covariance structure
if (!all(c("beta.CAE", "se.CAE") %in% names(SNPs))) {
  stop("Expected CAE columns are missing.")
}

SNPs[, c("beta.CAE", "se.CAE") := NULL]

total_snps <- nrow(SNPs)

cat(
  "Total chr2 SNPs: ",
  format(total_snps, big.mark = ","),
  "\n",
  sep = ""
)

# Split chromosome 2 into two approximately equal halves
split_point <- ceiling(total_snps / 2)

if (part == 1L) {
  idx <- 1:split_point
} else {
  idx <- (split_point + 1L):total_snps
}

SNPs_part <- SNPs[idx]

cat(
  "Part ",
  part,
  " rows: ",
  format(nrow(SNPs_part), big.mark = ","),
  "\n",
  sep = ""
)

cat(
  "Original row range: ",
  min(idx),
  " - ",
  max(idx),
  "\n",
  sep = ""
)

rm(SNPs)
invisible(gc())

start_time <- Sys.time()

result <- GenomicSEM::userGWAS(
  covstruc = covstruc,
  SNPs = SNPs_part,
  estimation = "DWLS",
  model = model,
  printwarn = TRUE,
  sub = requested_paths,
  cores = 1L,
  toler = FALSE,
  SNPSE = FALSE,
  parallel = FALSE,
  GC = "standard",
  MPI = FALSE,
  smooth_check = TRUE,
  std.lv = FALSE,
  fix_measurement = TRUE,
  Q_SNP = TRUE
)

end_time <- Sys.time()

runtime_minutes <- as.numeric(
  difftime(
    end_time,
    start_time,
    units = "mins"
  )
)

if (!is.list(result) ||
    length(result) != length(requested_paths)) {
  stop("Unexpected userGWAS result structure.")
}

# Save each path separately.
# We deliberately avoid saving a redundant all_paths RDS here
# to reduce scratch-space use.
for (i in seq_along(result)) {

  result_dt <- as.data.table(
    result[[i]]
  )

  out_file <- file.path(
    result_dir,
    sprintf(
      "GGE_no_CAE_chr02_part%d_%s.csv.gz",
      part,
      path_labels[[i]]
    )
  )

  fwrite(
    result_dt,
    out_file,
    sep = ",",
    quote = FALSE,
    na = "NA",
    compress = "gzip"
  )
}

summary_dt <- data.table(
  part = part,
  total_chr2_snps = total_snps,
  first_original_row = min(idx),
  last_original_row = max(idx),
  rows_analyzed = length(idx),
  runtime_minutes = round(runtime_minutes, 2),
  GenomicSEM_version =
    as.character(packageVersion("GenomicSEM"))
)

fwrite(
  summary_dt,
  file.path(
    result_dir,
    sprintf(
      "GGE_no_CAE_chr02_part%d_summary.csv",
      part
    )
  )
)

writeLines(
  capture.output(sessionInfo()),
  file.path(
    result_dir,
    sprintf(
      "GGE_no_CAE_chr02_part%d_session_info.txt",
      part
    )
  )
)

done_file <- file.path(
  result_dir,
  sprintf(
    "GGE_no_CAE_chr02_part%d.done",
    part
  )
)

writeLines(
  c(
    paste0("part=", part),
    paste0("rows=", length(idx)),
    paste0(
      "runtime_minutes=",
      round(runtime_minutes, 2)
    ),
    paste0(
      "completed=",
      format(
        Sys.time(),
        tz = "Africa/Johannesburg"
      )
    )
  ),
  done_file
)

cat("\nCompleted successfully.\n")
cat("Part: ", part, "\n", sep = "")
cat(
  "Runtime: ",
  round(runtime_minutes, 2),
  " minutes\n",
  sep = ""
)
cat(
  "Completion marker: ",
  done_file,
  "\n",
  sep = ""
)
