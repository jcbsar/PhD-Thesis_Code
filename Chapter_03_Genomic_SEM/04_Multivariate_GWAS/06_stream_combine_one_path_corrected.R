# ============================================================
# 1. Command-line arguments
# ============================================================

args <- commandArgs(
  trailingOnly = TRUE
)

if (length(args) != 3L) {
  stop(
    paste0(
      "Usage:\n",
      "Rscript 06_stream_combine_one_path.R ",
      "<GGE|FE> <path_index> <path_label>"
    )
  )
}

analysis <- toupper(
  args[[1]]
)

path_index <- as.integer(
  args[[2]]
)

path_label <- args[[3]]


if (!analysis %in% c("GGE", "FE")) {
  stop(
    "Analysis must be either GGE or FE."
  )
}

if (
  is.na(path_index) ||
  !path_index %in% c(1L, 2L)
) {
  stop(
    "Path index must be 1 or 2."
  )
}


# ============================================================
# 2. Directories
# ============================================================

root_dir <- "/scratch/jcbsar002/gsem"

results_root <- file.path(
  root_dir,
  "results"
)

split_root <- file.path(
  root_dir,
  "input_by_chr"
)

output_dir <- file.path(
  root_dir,
  "final_streamed_corrected"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 3. Analysis-specific paths
# ============================================================

if (analysis == "GGE") {

  result_dir <- file.path(
    results_root,
    "GGE_brain_corrected"
  )

  split_summary_file <- file.path(
    split_root,
    "GGE_brain",
    "GGE_split_summary.csv"
  )

} else {

  result_dir <- file.path(
    results_root,
    "FE_brain"
  )

  split_summary_file <- file.path(
    split_root,
    "FE_brain",
    "FE_split_summary.csv"
  )
}


if (!dir.exists(result_dir)) {
  stop(
    paste0(
      "Cannot find chromosome result directory:\n",
      result_dir
    )
  )
}

if (!file.exists(split_summary_file)) {
  stop(
    paste0(
      "Cannot find chromosome split summary:\n",
      split_summary_file
    )
  )
}


# ============================================================
# 4. Read and validate chromosome split summary
# ============================================================

split_summary <- fread(
  split_summary_file
)

required_split_columns <- c(
  "chromosome",
  "rows"
)

missing_split_columns <- setdiff(
  required_split_columns,
  names(split_summary)
)

if (length(missing_split_columns) > 0L) {
  stop(
    paste0(
      "The split summary is missing these columns:\n",
      paste(
        missing_split_columns,
        collapse = ", "
      )
    )
  )
}

split_summary[
  ,
  chromosome := as.integer(
    chromosome
  )
]

split_summary[
  ,
  rows := as.integer(
    rows
  )
]

setorder(
  split_summary,
  chromosome
)

if (!identical(
  split_summary$chromosome,
  1:22
)) {
  stop(
    "The split summary does not contain chromosomes 1-22 in order."
  )
}

if (
  anyNA(split_summary$rows) ||
  any(split_summary$rows <= 0L)
) {
  stop(
    "The split summary contains invalid chromosome row counts."
  )
}

expected_total_rows <- sum(
  split_summary$rows
)


# ============================================================
# 5. Output files
# ============================================================

output_prefix <- paste0(
  analysis,
  "_brain_",
  path_label,
  "_userGWAS_full_streamed"
)

temporary_csv <- file.path(
  output_dir,
  paste0(
    output_prefix,
    ".csv.partial"
  )
)

final_csv <- file.path(
  output_dir,
  paste0(
    output_prefix,
    ".csv"
  )
)

validation_file <- file.path(
  output_dir,
  paste0(
    output_prefix,
    "_validation.csv"
  )
)

summary_file <- file.path(
  output_dir,
  paste0(
    output_prefix,
    "_summary.csv"
  )
)

session_file <- file.path(
  output_dir,
  paste0(
    output_prefix,
    "_session_info.txt"
  )
)


# Remove leftovers from an earlier failed or incomplete attempt.

for (old_file in c(
  temporary_csv,
  final_csv
)) {

  if (file.exists(old_file)) {
    file.remove(
      old_file
    )
  }
}


# ============================================================
# 6. Helper functions
# ============================================================

count_nonempty_errors <- function(
  result_dt
) {

  if (!"error" %in% names(result_dt)) {
    return(
      NA_integer_
    )
  }

  error_values <- as.character(
    result_dt$error
  )

  sum(
    !is.na(error_values) &
      !error_values %in% c(
        "",
        "0",
        "FALSE",
        "NA"
      )
  )
}


extract_warning_text <- function(
  result_dt
) {

  if (!"warning" %in% names(result_dt)) {
    return(
      ""
    )
  }

  warning_values <- unique(
    as.character(
      result_dt$warning
    )
  )

  warning_values <- warning_values[
    !is.na(warning_values) &
      !warning_values %in% c(
        "",
        "0",
        "FALSE",
        "NA"
      )
  ]

  if (length(warning_values) == 0L) {
    return(
      ""
    )
  }

  warning_values <- gsub(
    pattern = "[\r\n\t]+",
    replacement = " ",
    x = warning_values
  )

  paste(
    warning_values,
    collapse = " || "
  )
}


clean_character_columns <- function(
  result_dt
) {

  character_columns <- names(
    result_dt
  )[
    vapply(
      result_dt,
      is.character,
      logical(1)
    )
  ]

  for (column_name in character_columns) {

    cleaned_values <- gsub(
      pattern = "[\r\n\t]+",
      replacement = " ",
      x = result_dt[[column_name]]
    )

    set(
      result_dt,
      j = column_name,
      value = cleaned_values
    )
  }

  result_dt
}


# ============================================================
# 7. Initialise run-level diagnostics
# ============================================================

validation_list <- vector(
  mode = "list",
  length = 22L
)

total_rows <- 0L

total_significant <- 0L

overall_max_z_smooth <- NA_real_


cat("\n")
cat("====================================================\n")
cat("STREAMING USERGWAS COMBINATION\n")
cat("====================================================\n")
cat("Analysis: ", analysis, "\n", sep = "")
cat("Path index: ", path_index, "\n", sep = "")
cat("Path label: ", path_label, "\n", sep = "")
cat(
  "Expected total rows: ",
  format(
    expected_total_rows,
    big.mark = ","
  ),
  "\n",
  sep = ""
)
cat("Output directory: ", output_dir, "\n", sep = "")
cat("====================================================\n\n")


# ============================================================
# 8. Process one chromosome at a time
# ============================================================

for (chr in 1:22) {

  chromosome_file <- file.path(
    result_dir,
    sprintf(
      "%s_chr%02d_all_paths.rds",
      analysis,
      chr
    )
  )

  if (!file.exists(chromosome_file)) {
    stop(
      paste0(
        "Missing chromosome result file:\n",
        chromosome_file
      )
    )
  }


  # ----------------------------------------------------------
  # Read the chromosome-level userGWAS object
  # ----------------------------------------------------------

  result_object <- readRDS(
    chromosome_file
  )

  if (
    !is.list(result_object) ||
    length(result_object) < path_index
  ) {
    stop(
      paste0(
        "Unexpected result structure in:\n",
        chromosome_file
      )
    )
  }

  result_dt <- as.data.table(
    result_object[[path_index]]
  )


  # ----------------------------------------------------------
  # Expected chromosome row count
  # ----------------------------------------------------------

  expected_position <- match(
    chr,
    split_summary$chromosome
  )

  if (is.na(expected_position)) {
    stop(
      "Chromosome ",
      chr,
      " is absent from the split summary."
    )
  }

  expected_chr_rows <- split_summary$rows[
    expected_position
  ]


  # ----------------------------------------------------------
  # Validate row count
  # ----------------------------------------------------------

  if (nrow(result_dt) != expected_chr_rows) {
    stop(
      paste0(
        analysis,
        " ",
        path_label,
        " chromosome ",
        chr,
        " has ",
        format(
          nrow(result_dt),
          big.mark = ","
        ),
        " rows; expected ",
        format(
          expected_chr_rows,
          big.mark = ","
        ),
        "."
      )
    )
  }


  # ----------------------------------------------------------
  # Validate required columns
  # ----------------------------------------------------------

  required_result_columns <- c(
    "SNP",
    "CHR",
    "BP",
    "MAF",
    "A1",
    "A2",
    "lhs",
    "op",
    "rhs",
    "est",
    "SE",
    "Z_Estimate",
    "Pval_Estimate"
  )

  missing_result_columns <- setdiff(
    required_result_columns,
    names(result_dt)
  )

  if (length(missing_result_columns) > 0L) {
    stop(
      paste0(
        analysis,
        " ",
        path_label,
        " chromosome ",
        chr,
        " is missing these columns:\n",
        paste(
          missing_result_columns,
          collapse = ", "
        )
      )
    )
  }


  # ----------------------------------------------------------
  # Validate chromosome assignment
  # ----------------------------------------------------------

  observed_chr <- unique(
    as.integer(
      result_dt$CHR[
        !is.na(
          result_dt$CHR
        )
      ]
    )
  )

  if (
    length(observed_chr) != 1L ||
    observed_chr != chr
  ) {
    stop(
      paste0(
        "Invalid CHR values in ",
        analysis,
        " ",
        path_label,
        " chromosome ",
        chr,
        "."
      )
    )
  }


  # ----------------------------------------------------------
  # Validate SNP identifiers
  # ----------------------------------------------------------

  missing_snps <- sum(
    is.na(
      result_dt$SNP
    ) |
      result_dt$SNP == ""
  )

  if (missing_snps > 0L) {
    stop(
      paste0(
        analysis,
        " ",
        path_label,
        " chromosome ",
        chr,
        " contains ",
        missing_snps,
        " missing SNP identifiers."
      )
    )
  }

  duplicate_snps <- sum(
    duplicated(
      result_dt$SNP
    )
  )

  if (duplicate_snps > 0L) {
    stop(
      paste0(
        analysis,
        " ",
        path_label,
        " chromosome ",
        chr,
        " contains ",
        duplicate_snps,
        " duplicated SNP identifiers."
      )
    )
  }


  # ----------------------------------------------------------
  # Check model errors
  # ----------------------------------------------------------

  error_count <- count_nonempty_errors(
    result_dt
  )

  if (
    is.finite(error_count) &&
    error_count > 0L
  ) {
    stop(
      paste0(
        analysis,
        " ",
        path_label,
        " chromosome ",
        chr,
        " contains ",
        error_count,
        " nonempty model errors."
      )
    )
  }


  # ----------------------------------------------------------
  # Retain distinct warning text in validation only
  # ----------------------------------------------------------

  warning_text <- extract_warning_text(
    result_dt
  )


  # ----------------------------------------------------------
  # Smoothing diagnostic
  # ----------------------------------------------------------

  chromosome_max_z_smooth <- NA_real_

  if ("Z_smooth" %in% names(result_dt)) {

    z_smooth_values <- suppressWarnings(
      as.numeric(
        result_dt$Z_smooth
      )
    )

    if (any(is.finite(z_smooth_values))) {

      chromosome_max_z_smooth <- max(
        abs(
          z_smooth_values
        ),
        na.rm = TRUE
      )

      if (
        is.na(overall_max_z_smooth) ||
        chromosome_max_z_smooth >
          overall_max_z_smooth
      ) {

        overall_max_z_smooth <-
          chromosome_max_z_smooth
      }
    }
  }


  # ----------------------------------------------------------
  # Genome-wide significant SNP count
  # ----------------------------------------------------------

  p_values <- suppressWarnings(
    as.numeric(
      result_dt$Pval_Estimate
    )
  )

  significant_snps <- sum(
    is.finite(p_values) &
      p_values < 5e-8
  )

  total_significant <-
    total_significant +
    significant_snps


  # ----------------------------------------------------------
  # Remove repeated warning/error fields from analysis output
  # ----------------------------------------------------------

  columns_to_remove <- intersect(
    c(
      "warning",
      "error"
    ),
    names(result_dt)
  )

  if (length(columns_to_remove) > 0L) {

    result_dt[
      ,
      (columns_to_remove) := NULL
    ]
  }


  # ----------------------------------------------------------
  # Remove line breaks and tabs from remaining text fields
  # ----------------------------------------------------------

  result_dt <- clean_character_columns(
    result_dt
  )


  # ----------------------------------------------------------
  # Stream this chromosome directly to the output file
  # ----------------------------------------------------------

  fwrite(
    result_dt,
    file = temporary_csv,
    sep = ",",
    append = chr > 1L,
    col.names = chr == 1L,
    quote = "auto",
    na = "NA"
  )


  # ----------------------------------------------------------
  # Record chromosome validation
  # ----------------------------------------------------------

  validation_list[[chr]] <- data.table(

    analysis = analysis,

    path = path_label,

    chromosome = chr,

    expected_rows =
      expected_chr_rows,

    observed_rows =
      nrow(result_dt),

    unique_snps =
      uniqueN(
        result_dt$SNP
      ),

    duplicated_snps =
      duplicate_snps,

    missing_snps =
      missing_snps,

    model_errors =
      error_count,

    genome_wide_significant_snps =
      significant_snps,

    maximum_Z_smooth =
      chromosome_max_z_smooth,

    distinct_warning_text =
      warning_text
  )


  total_rows <-
    total_rows +
    nrow(result_dt)


  cat(
    analysis,
    " ",
    path_label,
    " chromosome ",
    chr,
    ": ",
    format(
      nrow(result_dt),
      big.mark = ","
    ),
    " rows written\n",
    sep = ""
  )


  # ----------------------------------------------------------
  # Release chromosome memory
  # ----------------------------------------------------------

  rm(
    result_object,
    result_dt,
    p_values
  )

  invisible(
    gc()
  )
}


# ============================================================
# 9. Validate the completed output
# ============================================================

if (total_rows != expected_total_rows) {
  stop(
    paste0(
      "Final row count was ",
      format(
        total_rows,
        big.mark = ","
      ),
      "; expected ",
      format(
        expected_total_rows,
        big.mark = ","
      ),
      "."
    )
  )
}

if (!file.exists(temporary_csv)) {
  stop(
    "The temporary output CSV was not created."
  )
}

rename_success <- file.rename(
  from = temporary_csv,
  to = final_csv
)

if (!rename_success) {
  stop(
    paste0(
      "Could not rename the completed file to:\n",
      final_csv
    )
  )
}


# ============================================================
# 10. Save validation and summary files
# ============================================================

validation_dt <- rbindlist(
  validation_list,
  use.names = TRUE,
  fill = TRUE
)

fwrite(
  validation_dt,
  file = validation_file,
  quote = "auto"
)

run_summary <- data.table(

  analysis = analysis,

  path = path_label,

  expected_rows =
    expected_total_rows,

  observed_rows =
    total_rows,

  chromosomes =
    22L,

  duplicated_snps =
    sum(
      validation_dt$duplicated_snps
    ),

  missing_snps =
    sum(
      validation_dt$missing_snps
    ),

  model_errors =
    sum(
      validation_dt$model_errors,
      na.rm = TRUE
    ),

  genome_wide_significant_snps =
    total_significant,

  maximum_Z_smooth =
    overall_max_z_smooth,

  output_file =
    final_csv,

  output_size_mb =
    round(
      file.info(final_csv)$size /
        1024^2,
      2
    )
)

fwrite(
  run_summary,
  file = summary_file
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  session_file
)


# ============================================================
# 11. Final output
# ============================================================

cat("\n")
cat("====================================================\n")
cat("STREAMING COMBINATION COMPLETE\n")
cat("====================================================\n")

print(
  run_summary
)

cat(
  "\nCSV created:\n",
  final_csv,
  "\n",
  sep = ""
)

cat(
  "\nValidation file:\n",
  validation_file,
  "\n",
  sep = ""
)
