#!/usr/bin/env Rscript
# Verify that retained six-region GenomicSEM covariance inputs agree with
# the matching phenotype subsets of archived July 2026 LDSC objects.
# This script does not run LDSC, smooth covariance matrices, or fit CFA models.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(paste(
    "Usage: Rscript 01_verify_six_region_objects.R",
    "GGE_July.rds GGE_6region.rds FE_July.rds FE_6region.rds"
  ))
}

read_ldsc <- function(path) {
  if (!file.exists(path)) stop("Input file not found: ", path)
  obj <- readRDS(path)
  needed <- c("S", "V", "I", "N")
  absent <- setdiff(needed, names(obj))
  if (length(absent)) {
    stop(basename(path), " lacks: ", paste(absent, collapse = ", "))
  }
  obj$S <- as.matrix(obj$S)
  obj$V <- as.matrix(obj$V)
  obj$I <- as.matrix(obj$I)
  n <- ncol(obj$S)
  tr <- colnames(obj$S)
  if (is.null(tr) || anyNA(tr) || anyDuplicated(tr)) {
    stop("S requires unique, non-missing phenotype column names: ", path)
  }
  if (!identical(dim(obj$S), c(n, n)) ||
      !identical(dim(obj$I), c(n, n))) {
    stop("S and I dimensions must be square and consistent: ", path)
  }
k <- as.integer(n * (n + 1L) / 2L)

if (!all(dim(obj$V) == c(k, k)) || length(obj$N) != k) {
  stop("V/N lengths do not match lower triangle of S: ", path)
}
  for (field in needed) {
    if (!is.numeric(obj[[field]]) || any(!is.finite(obj[[field]]))) {
      stop("Non-numeric or non-finite ", field, " in: ", path)
    }
  }
  obj
}

# GenomicSEM's V and N index the lower triangle of S, including the diagonal,
# using R column-major order. Translate the retained ordering into the source.

# Map retained covariance elements to their original LDSC positions.
vech_index <- function(all_traits, retained_traits) {

  n <- length(all_traits)
  k <- length(retained_traits)

  source_pairs <- which(
    lower.tri(matrix(FALSE, n, n), diag = TRUE),
    arr.ind = TRUE
  )

  retained_pairs <- which(
    lower.tri(matrix(FALSE, k, k), diag = TRUE),
    arr.ind = TRUE
  )

  source_row <- match(
    retained_traits[retained_pairs[, 1]],
    all_traits
  )

  source_col <- match(
    retained_traits[retained_pairs[, 2]],
    all_traits
  )

  if (anyNA(source_row) || anyNA(source_col)) {
    stop(
      "At least one retained phenotype is missing from ",
      "the source LDSC object"
    )
  }

  source_keys <- paste(
    source_pairs[, 1],
    source_pairs[, 2],
    sep = ":"
  )

  lookup <- paste(
    pmax(source_row, source_col),
    pmin(source_row, source_col),
    sep = ":"
  )

  idx <- match(lookup, source_keys)

  if (anyNA(idx) ||
      length(idx) != k * (k + 1L) / 2L) {
    stop(
      "Could not map all retained covariance-element indices"
    )
  }

  idx
}


compare_pair <- function(group, full_path, retained_path, expected_phenotypes) {
  full <- read_ldsc(full_path)
  retained <- read_ldsc(retained_path)
  original_traits <- colnames(full$S)
  kept_traits <- colnames(retained$S)

  if (!setequal(kept_traits, expected_phenotypes)) {
    stop(group, ": retained phenotype set does not match its specification")
  }
  if (!setequal(original_traits, c(expected_phenotypes, "brainstem"))) {
    stop(group, ": source phenotype set does not match its specification")
  }

  positions <- match(kept_traits, original_traits)
  vech_positions <- vech_index(original_traits, kept_traits)

  x <- list(
    S = full$S[positions, positions, drop = FALSE],
    I = full$I[positions, positions, drop = FALSE],
    V = full$V[vech_positions, vech_positions, drop = FALSE],
    N = as.numeric(full$N)[vech_positions]
  )
  y <- list(
    S = retained$S,
    I = retained$I,
    V = retained$V,
    N = as.numeric(retained$N)
  )

  output <- lapply(names(x), function(field) {
    a <- x[[field]]
    b <- y[[field]]
    if (!identical(dim(a), dim(b)) || length(a) != length(b)) {
      stop(group, ": incompatible dimensions for ", field)
    }
    delta <- max(abs(as.numeric(a) - as.numeric(b)))
    data.frame(
      model = group,
      field = field,
      max_absolute_difference = delta,
      exact_numeric_match = (delta == 0),
      tolerance_match_1e10 = (delta <= 1e-10),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, output)
}

results <- rbind(
  compare_pair(
    "GGE", args[[1]], args[[2]],
    c("CAE", "GTCS", "JAE", "JME",
      "accumbens", "amygdala", "caudate", "pallidum", "putamen", "thalamus")
  ),
  compare_pair(
    "FE", args[[3]], args[[4]],
    c("F_HS", "F_Neg", "F_Oth",
      "accumbens", "amygdala", "caudate", "pallidum", "putamen", "thalamus")
  )
)

print(results, row.names = FALSE)

command_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(command_file) != 1L) {
  stop("Run this script with Rscript so the output directory is defined")
}
script_file <- sub("^--file=", "", command_file)
section_dir <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
out_dir <- file.path(section_dir, "outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
report <- file.path(out_dir, "six_region_object_verification.csv")
write.csv(results, report, row.names = FALSE)
cat("\nVerification report: ", report, "\n", sep = "")

if (!all(results$tolerance_match_1e10)) {
  stop("The retained covariance inputs differ from the expected July LDSC subsets")
}
cat("PASS: all S, I, V and N values match the July-derived six-region subsets.\n")
