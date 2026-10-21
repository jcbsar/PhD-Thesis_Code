
cat("Select the new GGE CFA script\n")
script <- file.choose()

cat("Select ldsc_GGE_brain_6region.rds\n")
gge_ldsc <- file.choose()

invisible(parse(file = script))

rscript <- file.path(R.home("bin"), "Rscript.exe")

result <- suppressWarnings(
  system2(
    rscript,
    args = shQuote(c(script, gge_ldsc)),
    stdout = TRUE,
    stderr = TRUE
  )
)

cat(result, sep = "\n")

status <- attr(result, "status")
if (is.null(status)) status <- 0L

cat("\nExit status:", status, "\n")
