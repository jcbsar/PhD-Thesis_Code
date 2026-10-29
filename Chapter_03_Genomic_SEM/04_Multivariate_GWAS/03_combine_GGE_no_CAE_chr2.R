
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

root_dir <- "/scratch/jcbsar002/gsem"

parts_dir <- file.path(
  root_dir,
  "results",
  "GGE_no_CAE_sensitivity",
  "chr02_parts"
)

final_dir <- file.path(
  root_dir,
  "results",
  "GGE_no_CAE_sensitivity"
)

# ---------------------------------------------------------
# Function to repair the multiline warning field
# ---------------------------------------------------------

read_repaired_chunk <- function(file) {

  cat("Reading and repairing:\n", file, "\n")

  con <- gzfile(file, open = "rt")
  lines <- readLines(con, warn = FALSE)
  close(con)

  header <- lines[1]
  body <- lines[-1]

  # Every true SNP record begins:
  # SNP,2,BP,...
  record_starts <- which(
    grepl("^[^,]+,2,[0-9]+,", body)
  )

  if (length(record_starts) == 0L) {
    stop("No SNP records detected in: ", file)
  }

  record_ends <- c(
    record_starts[-1] - 1L,
    length(body)
  )

  # Join continuation lines back onto the SNP record
  records <- vapply(
    seq_along(record_starts),
    function(i) {
      paste(
        body[record_starts[i]:record_ends[i]],
        collapse = " "
      )
    },
    character(1)
  )

  repaired_text <- paste(
    c(header, records),
    collapse = "\n"
  )

  dat <- fread(
    text = repaired_text,
    sep = ",",
    header = TRUE,
    quote = "",
    fill = FALSE
  )

  cat(
    "Recovered ",
    format(nrow(dat), big.mark = ","),
    " SNP rows and ",
    ncol(dat),
    " columns\n",
    sep = ""
  )

  dat
}

# ---------------------------------------------------------
# Repair and combine GGE + Brain outputs
# ---------------------------------------------------------

paths <- c("GGE", "Brain")

combined_results <- vector(
  "list",
  length(paths)
)

names(combined_results) <- paths

for (path in paths) {

  cat("\n========================================\n")
  cat("Processing ", path, "\n", sep = "")
  cat("========================================\n")

  part1_file <- file.path(
    parts_dir,
    paste0(
      "GGE_no_CAE_chr02_part1_",
      path,
      ".csv.gz"
    )
  )

  part2_file <- file.path(
    parts_dir,
    paste0(
      "GGE_no_CAE_chr02_part2_",
      path,
      ".csv.gz"
    )
  )

  if (!file.exists(part1_file)) {
    stop("Missing: ", part1_file)
  }

  if (!file.exists(part2_file)) {
    stop("Missing: ", part2_file)
  }

  part1 <- read_repaired_chunk(part1_file)
  part2 <- read_repaired_chunk(part2_file)

  if (nrow(part1) != 205347L) {
    stop(
      path,
      " part 1 has ",
      nrow(part1),
      " rows; expected 205,347."
    )
  }

  if (nrow(part2) != 205347L) {
    stop(
      path,
      " part 2 has ",
      nrow(part2),
      " rows; expected 205,347."
    )
  }

  combined <- rbindlist(
    list(part1, part2),
    use.names = TRUE
  )

  if (nrow(combined) != 410694L) {
    stop(
      path,
      " combined row count is ",
      nrow(combined),
      "; expected 410,694."
    )
  }

  if (uniqueN(combined$SNP) != 410694L) {
    stop(
      path,
      " does not contain 410,694 unique SNPs."
    )
  }

  combined_results[[path]] <- combined

  # Write a CLEAN CSV.
  # quote=TRUE protects warning text properly.
  out_file <- file.path(
    final_dir,
    paste0(
      "GGE_no_CAE_chr02_",
      path,
      ".csv.gz"
    )
  )

  fwrite(
    combined,
    out_file,
    sep = ",",
    quote = TRUE,
    na = "NA",
    compress = "gzip"
  )

  cat(
    path,
    " successfully combined: ",
    format(nrow(combined), big.mark = ","),
    " SNPs\n",
    sep = ""
  )
}

# ---------------------------------------------------------
# Save an RDS matching the other chromosome outputs
# ---------------------------------------------------------

all_paths_file <- file.path(
  final_dir,
  "GGE_no_CAE_chr02_all_paths.rds"
)

saveRDS(
  unname(combined_results),
  all_paths_file,
  compress = "gzip"
)

# ---------------------------------------------------------
# Completion marker
# ---------------------------------------------------------

done_file <- file.path(
  final_dir,
  "GGE_no_CAE_chr02.done"
)

writeLines(
  c(
    "analysis=GGE_no_CAE",
    "chromosome=2",
    "method=combined_from_two_single_core_chunks",
    "snps=410694",
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

writeLines(
  capture.output(sessionInfo()),
  file.path(
    final_dir,
    "GGE_no_CAE_chr02_session_info.txt"
  )
)

cat("\n========================================\n")
cat("CHROMOSOME 2 COMBINED SUCCESSFULLY\n")
cat("========================================\n")
cat("SNPs: 410,694\n")
cat("Completion marker: ", done_file, "\n", sep = "")
