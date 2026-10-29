
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

root_dir <- "/scratch/jcbsar002/gsem"

result_dir <- file.path(
  root_dir,
  "results",
  "GGE_no_CAE_sensitivity"
)

final_dir <- file.path(
  root_dir,
  "final_streamed"
)

dir.create(
  final_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

paths <- c("GGE", "Brain")

for (path in paths) {

  cat("\n========================================\n")
  cat("Combining all chromosomes for ", path, "\n", sep = "")
  cat("========================================\n")

  chr_results <- vector("list", 22L)

  for (chr in 1:22) {

    chr_file <- file.path(
      result_dir,
      sprintf(
        "GGE_no_CAE_chr%02d_all_paths.rds",
        chr
      )
    )

    if (!file.exists(chr_file)) {
      stop("Missing chromosome RDS: ", chr_file)
    }

    x <- readRDS(chr_file)

    if (!is.list(x) || length(x) < 2L) {
      stop(
        "Unexpected RDS structure for chromosome ",
        chr
      )
    }

    # userGWAS output order:
    # [[1]] = GGE~SNP
    # [[2]] = Brain~SNP
    idx <- if (path == "GGE") 1L else 2L

    chr_dt <- as.data.table(
      x[[idx]]
    )

    chr_results[[chr]] <- chr_dt

    cat(
      "chr",
      sprintf("%02d", chr),
      ": ",
      format(nrow(chr_dt), big.mark = ","),
      " rows\n",
      sep = ""
    )

    rm(x, chr_dt)
    invisible(gc())
  }

  combined <- rbindlist(
    chr_results,
    use.names = TRUE,
    fill = TRUE
  )

  setorder(
    combined,
    CHR,
    BP
  )

  cat(
    "\nCombined ",
    path,
    " rows: ",
    format(nrow(combined), big.mark = ","),
    "\n",
    sep = ""
  )

  cat(
    "Unique SNPs: ",
    format(uniqueN(combined$SNP), big.mark = ","),
    "\n",
    sep = ""
  )

  out_file <- file.path(
    final_dir,
    paste0(
      "GGE_no_CAE_",
      path,
      "_userGWAS_full_streamed.csv.gz"
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
    "Saved: ",
    out_file,
    "\n",
    sep = ""
  )

  rm(combined, chr_results)
  invisible(gc())
}

writeLines(
  capture.output(sessionInfo()),
  file.path(
    final_dir,
    "GGE_no_CAE_userGWAS_full_streamed_session_info.txt"
  )
)

cat("\n========================================\n")
cat("ALL CHROMOSOMES COMBINED SUCCESSFULLY\n")
cat("========================================\n")
