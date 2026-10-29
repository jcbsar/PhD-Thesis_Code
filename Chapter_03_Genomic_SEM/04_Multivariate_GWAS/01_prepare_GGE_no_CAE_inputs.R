
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GenomicSEM)
  library(data.table)
})

root_gsem <- "/scratch/jcbsar002/gsem"
source_root <- "/scratch/jcbsar002/enigma_ilae_sumstats/gsem"

input_dir <- file.path(root_gsem, "input")

ld_dir <- file.path(
  source_root,
  "ldsc_refdata",
  "eur_w_ld_chr"
)

munged_dir <- file.path(source_root, "munged")

cat("========================================\n")
cat("Preparing GGE no-CAE sensitivity inputs\n")
cat("========================================\n\n")

############################################################
# PART 1: CREATE 9-TRAIT LDSC OBJECT
############################################################

traits <- c(
  file.path(munged_dir, "GTCS_cau_21_Jan_25.hm3.sumstats.gz"),
  file.path(munged_dir, "JAE_cau_21_Jan_25.hm3.sumstats.gz"),
  file.path(munged_dir, "JME_cau_21_Jan_25.hm3.sumstats.gz"),
  file.path(munged_dir, "accumbens_21_Jan_25.hm3.sumstats.gz"),
  file.path(munged_dir, "amygdala_21_Jan_25.hm3.sumstats.gz"),
  file.path(munged_dir, "caudate_21_Jan_25.hm3.sumstats.gz"),
  file.path(munged_dir, "pallidum_21_Jan_25.hm3.sumstats.gz"),
  file.path(munged_dir, "putamen_21_Jan_25.hm3.sumstats.gz"),
  file.path(munged_dir, "thalamus_21_Jan_25.hm3.sumstats.gz")
)

trait.names <- c(
  "GTCS", "JAE", "JME",
  "accumbens", "amygdala", "caudate",
  "pallidum", "putamen", "thalamus"
)

sample.prev <- c(
  0.50, 0.50, 0.50,
  NA, NA, NA, NA, NA, NA
)

population.prev <- c(
  0.00020, 0.00015, 0.00015,
  NA, NA, NA, NA, NA, NA
)

cat("Checking LDSC input files...\n")

missing_traits <- traits[!file.exists(traits)]

if (length(missing_traits) > 0L) {
  stop(
    "Missing munged trait files:\n",
    paste(missing_traits, collapse = "\n")
  )
}

if (!dir.exists(ld_dir)) {
  stop("Missing LD-score directory: ", ld_dir)
}

cat("All LDSC files found.\n\n")
cat("Running 9-trait LDSC...\n")

ldsc_no_CAE <- GenomicSEM::ldsc(
  traits = traits,
  sample.prev = sample.prev,
  population.prev = population.prev,
  ld = paste0(ld_dir, "/"),
  wld = paste0(ld_dir, "/"),
  trait.names = trait.names
)

ldsc_out <- file.path(
  input_dir,
  "ldsc_GGE_no_CAE_brain_6region.rds"
)

saveRDS(ldsc_no_CAE, ldsc_out)

cat("\nSaved LDSC object:\n")
cat(ldsc_out, "\n")

cat("\nLDSC dimensions:\n")
cat("S:", paste(dim(ldsc_no_CAE$S), collapse = " x "), "\n")
cat("V:", paste(dim(ldsc_no_CAE$V), collapse = " x "), "\n")

if (!identical(dim(ldsc_no_CAE$S), c(9L, 9L))) {
  stop("Unexpected S-matrix dimensions.")
}

if (!identical(dim(ldsc_no_CAE$V), c(45L, 45L))) {
  stop("Unexpected V-matrix dimensions.")
}

############################################################
# PART 2: CREATE MATCHED 9-TRAIT SNP OBJECT
############################################################

primary_sumstats_file <- file.path(
  input_dir,
  "GGE_brain_sumstats.rds"
)

if (!file.exists(primary_sumstats_file)) {
  stop(
    "Missing primary GGE sumstats object: ",
    primary_sumstats_file
  )
}

cat("\nLoading primary GGE SNP object...\n")

SNPs <- as.data.table(
  readRDS(primary_sumstats_file)
)

cat(
  "Primary object:",
  format(nrow(SNPs), big.mark = ","),
  "rows and",
  ncol(SNPs),
  "columns\n"
)

required_CAE_columns <- c(
  "beta.CAE",
  "se.CAE"
)

if (!all(required_CAE_columns %in% names(SNPs))) {
  stop(
    "Expected CAE columns not found: ",
    paste(required_CAE_columns, collapse = ", ")
  )
}

# Remove CAE SNP effects only.
# All remaining SNPs and harmonisation are identical
# to the primary analysis.
SNPs[, c("beta.CAE", "se.CAE") := NULL]

cat(
  "No-CAE object:",
  format(nrow(SNPs), big.mark = ","),
  "rows and",
  ncol(SNPs),
  "columns\n"
)

expected_columns <- c(
  "SNP",
  "CHR",
  "BP",
  "MAF",
  "A1",
  "A2",
  "beta.GTCS",
  "se.GTCS",
  "beta.JAE",
  "se.JAE",
  "beta.JME",
  "se.JME",
  "beta.accumbens",
  "se.accumbens",
  "beta.amygdala",
  "se.amygdala",
  "beta.caudate",
  "se.caudate",
  "beta.pallidum",
  "se.pallidum",
  "beta.putamen",
  "se.putamen",
  "beta.thalamus",
  "se.thalamus"
)

if (!identical(names(SNPs), expected_columns)) {
  stop(
    "Unexpected SNP column structure after removing CAE.\nObserved:\n",
    paste(names(SNPs), collapse = "\n")
  )
}

if (nrow(SNPs) != 4660537L) {
  stop(
    "Unexpected SNP count. Expected 4,660,537 but observed ",
    nrow(SNPs)
  )
}

if (ncol(SNPs) != 24L) {
  stop(
    "Unexpected number of columns. Expected 24 but observed ",
    ncol(SNPs)
  )
}

sumstats_out <- file.path(
  input_dir,
  "GGE_no_CAE_brain_sumstats.rds"
)

saveRDS(
  SNPs,
  sumstats_out,
  compress = "gzip"
)

cat("\nSaved no-CAE SNP object:\n")
cat(sumstats_out, "\n")

############################################################
# SESSION INFO
############################################################

session_file <- file.path(
  input_dir,
  "GGE_no_CAE_input_preparation_session_info.txt"
)

writeLines(
  capture.output(sessionInfo()),
  session_file
)

cat("\n========================================\n")
cat("PREPARATION COMPLETED SUCCESSFULLY\n")
cat("========================================\n")
cat("LDSC object: ", ldsc_out, "\n", sep = "")
cat("SNP object:  ", sumstats_out, "\n", sep = "")
cat("SNP rows:    ", format(nrow(SNPs), big.mark = ","), "\n", sep = "")
cat("SNP columns: ", ncol(SNPs), "\n", sep = "")
cat("S matrix:    ", paste(dim(ldsc_no_CAE$S), collapse = " x "), "\n", sep = "")
cat("V matrix:    ", paste(dim(ldsc_no_CAE$V), collapse = " x "), "\n", sep = "")
