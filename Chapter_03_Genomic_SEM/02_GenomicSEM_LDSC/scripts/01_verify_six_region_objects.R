
#!/usr/bin/env Rscript

# Verify that the retained six-region Genomic SEM covariance inputs
# agree with the corresponding phenotype subsets of archived July
# 2026 LDSC objects.
#
# The script compares S, I, V and N, accounting for the ordering of
# unique covariance elements in the sampling covariance matrix.
#
# This script does not rerun LDSC, smooth matrices, modify the
# original RDS files or fit factor models.

# ==============================================================================
# 1. Command-line arguments
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4L) {
  stop(
    paste(
      "Usage: Rscript 01_verify_six_region_objects.R",
      "GGE_July.rds GGE_6region.rds FE_July.rds FE_6region.rds"
    )
  )
}


# ==============================================================================
# 2. Read and validate LDSC objects
# ==============================================================================

read_ldsc <- function(path) {

  if (!file.exists(path)) {
    stop("Input file not found: ", path)
  }

  obj <- readRDS(path)

  needed <- c("S", "V", "I", "N")
  absent <- setdiff(needed, names(obj))

  if (length(absent) > 0L) {
    stop(
      basename(path),
      " lacks required components: ",
      paste(absent, collapse = ", ")
    )
  }

  obj$S <- as.matrix(obj$S)
  obj$V <- as.matrix(obj$V)
  obj$I <- as.matrix(obj$I)

  n <- ncol(obj$S)
  tr <- colnames(obj$S)

  if (
    is.null(tr) ||
    anyNA(tr) ||
    anyDuplicated(tr) > 0L
  ) {
    stop(
      "S requires unique, non-missing phenotype column names: ",
      path
    )
  }

  if (
    !identical(dim(obj$S), c(n, n)) ||
    !identical(dim(obj$I), c(n, n))
  ) {
    stop(
      "S and I dimensions must be square and consistent: ",
      path
    )
  }

  # Number of unique covariance elements, including the diagonal.
  # Use integer type to match the type returned by dim().

  k <- as.integer(n * (n + 1L) / 2L)

  if (!identical(dim(obj$V), c(k, k))) {
    stop(
      "Unexpected V dimensions in: ", path,
      ". Expected ", k, " x ", k,
      "; observed ",
      paste(dim(obj$V), collapse = " x ")
    )
  }

  if (length(obj$N) != k) {
    stop(
      "Unexpected N length in: ", path,
      ". Expected ", k,
      "; observed ", length(obj$N)
    )
  }

  for (field in needed) {

    if (
      !is.numeric(obj[[field]]) ||
      any(!is.finite(obj[[field]]))
    ) {
      stop(
        "Non-numeric or non-finite values in ",
        field, ": ", path
      )
    }
  }

  obj
}


# ==============================================================================
# 3. Map lower-triangular covariance-element indices
# ==============================================================================

# Genomic SEM's V and N use the unique covariance elements of S,
# including its diagonal, ordered by the lower triangle in R's
# column-major order.
#
# This function maps the retained phenotype ordering to the
# corresponding positions in the original LDSC object.

vech_index <- function(all_traits, retained_traits) {

  n <- length(all_traits)
  k <- length(retained_traits)

  source_pairs <- which(
    lower.tri(
      matrix(FALSE, nrow = n, ncol = n),
      diag = TRUE
    ),
    arr.ind = TRUE
  )

  retained_pairs <- which(
    lower.tri(
      matrix(FALSE, nrow = k, ncol = k),
      diag = TRUE
    ),
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
      "At least one retained phenotype is missing ",
      "from the source LDSC object"
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

  expected_length <- as.integer(
    k * (k + 1L) / 2L
  )

  if (
    anyNA(idx) ||
    length(idx) != expected_length
  ) {
    stop(
      "Could not map all retained covariance-element indices"
    )
  }

  idx
}


# ==============================================================================
# 4. Compare original and retained covariance objects
# ==============================================================================

compare_pair <- function(
  group,
  full_path,
  retained_path,
  expected_phenotypes
) {

  cat("\n====================================================\n")
  cat("Verifying ", group, " covariance objects\n", sep = "")
  cat("====================================================\n")

  cat("Original: ", full_path, "\n", sep = "")
  cat("Retained: ", retained_path, "\n\n", sep = "")

  full <- read_ldsc(full_path)
  retained <- read_ldsc(retained_path)

  original_traits <- colnames(full$S)
  kept_traits <- colnames(retained$S)

  # Confirm that the retained object has the expected phenotypes.

  if (!setequal(kept_traits, expected_phenotypes)) {
    stop(
      group,
      ": retained phenotype set does not match ",
      "the expected six-region specification"
    )
  }

  # The original LDSC object must contain the retained phenotypes
  # plus brainstem.

  if (!setequal(
    original_traits,
    c(expected_phenotypes, "brainstem")
  )) {
    stop(
      group,
      ": original phenotype set does not match ",
      "the expected seven-region specification"
    )
  }

  cat("Original S dimensions: ",
      paste(dim(full$S), collapse = " x "), "\n", sep = "")

  cat("Retained S dimensions: ",
      paste(dim(retained$S), collapse = " x "), "\n", sep = "")

  # Match retained phenotypes to the original covariance object.

  positions <- match(
    kept_traits,
    original_traits
  )

  # Map the unique covariance-element indices used by V and N.

  vech_positions <- vech_index(
    original_traits,
    kept_traits
  )

  # Extract the corresponding original covariance elements.

  x <- list(
    S = full$S[
      positions,
      positions,
      drop = FALSE
    ],
    I = full$I[
      positions,
      positions,
      drop = FALSE
    ],
    V = full$V[
      vech_positions,
      vech_positions,
      drop = FALSE
    ],
    N = as.numeric(full$N)[vech_positions]
  )

  # Read the retained covariance elements.

  y <- list(
    S = retained$S,
    I = retained$I,
    V = retained$V,
    N = as.numeric(retained$N)
  )

  # Compare every numerical element of each component.

  output <- lapply(names(x), function(field) {

    a <- x[[field]]
    b <- y[[field]]

    if (
      !identical(dim(a), dim(b)) ||
      length(a) != length(b)
    ) {
      stop(
        group,
        ": incompatible dimensions for ",
        field
      )
    }

    delta <- max(
      abs(as.numeric(a) - as.numeric(b))
    )

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


# ==============================================================================
# 5. Verify GGE and FE covariance objects
# ==============================================================================

gge_traits <- c(
  "CAE",
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

fe_traits <- c(
  "F_HS",
  "F_Neg",
  "F_Oth",
  "accumbens",
  "amygdala",
  "caudate",
  "pallidum",
  "putamen",
  "thalamus"
)

results <- rbind(
  compare_pair(
    group = "GGE",
    full_path = args[[1]],
    retained_path = args[[2]],
    expected_phenotypes = gge_traits
  ),
  compare_pair(
    group = "FE",
    full_path = args[[3]],
    retained_path = args[[4]],
    expected_phenotypes = fe_traits
  )
)


# ==============================================================================
# 6. Display verification results
# ==============================================================================

cat("\n====================================================\n")
cat("Six-region covariance verification results\n")
cat("====================================================\n\n")

print(
  results,
  row.names = FALSE
)


# ==============================================================================
# 7. Save verification report
# ==============================================================================

command_file <- grep(
  "^--file=",
  commandArgs(trailingOnly = FALSE),
  value = TRUE
)

if (length(command_file) != 1L) {
  stop(
    "Run this script using Rscript ",
    "so the output directory can be identified"
  )
}

script_file <- sub(
  "^--file=",
  "",
  command_file[[1]]
)

section_dir <- normalizePath(
  file.path(
    dirname(script_file),
    ".."
  ),
  mustWork = TRUE
)

out_dir <- file.path(
  section_dir,
  "outputs"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

report <- file.path(
  out_dir,
  "six_region_object_verification.csv"
)

write.csv(
  results,
  file = report,
  row.names = FALSE
)

cat(
  "\nVerification report: ",
  report,
  "\n",
  sep = ""
)


# ==============================================================================
# 8. Final verification status
# ==============================================================================

if (!all(results$tolerance_match_1e10)) {

  stop(
    "\nFAIL: At least one retained covariance component ",
    "differs from the corresponding July LDSC subset. ",
    "Check that the original July objects were selected."
  )
}

cat(
  "\nPASS: All S, I, V and N values match ",
  "the July-derived six-region subsets.\n"
)
