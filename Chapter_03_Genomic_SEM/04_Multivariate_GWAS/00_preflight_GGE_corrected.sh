
#!/bin/bash
#SBATCH --job-name=GGE_preflight
#SBATCH --account=psychiatry
#SBATCH --partition=ada
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --output=/scratch/jcbsar002/gsem/logs/GGE_preflight_%j.out

set -euo pipefail
export R_LIBS_USER="$HOME/R/4.6"

Rscript - <<'RSCRIPT'
root <- "/scratch/jcbsar002/gsem"

file <- file.path(
  root, "input",
  "ldsc_GGE_brain_6region_corrected.rds"
)

x <- readRDS(file)

stopifnot(
  identical(dim(x$S), c(10L, 10L)),
  identical(dim(x$I), c(10L, 10L)),
  identical(dim(x$V), c(55L, 55L)),
  abs(x$S[1, 1] - 0.2404654) < 1e-5,
  abs(x$S[4, 4] - 0.2534618) < 1e-5
)

cat("PASS: Corrected GGE LDSC covariance object\n")
cat("CAE variance:", x$S[1, 1], "\n")
cat("JME variance:", x$S[4, 4], "\n")

snp_files <- file.path(
  root, "input_by_chr", "GGE_brain",
  sprintf("chr%02d_sumstats.rds", 1:22)
)

missing <- snp_files[!file.exists(snp_files)]

if (length(missing) > 0) {
  stop("Missing chromosome inputs: ",
       paste(missing, collapse = ", "))
}

cat("PASS: All 22 chromosome SNP inputs exist\n")

out <- file.path(root, "results", "GGE_brain_corrected")

if (dir.exists(out) && length(list.files(out)) > 0) {
  stop("Corrected output directory is not empty.")
}

cat("PASS: Corrected output directory is available\n")
cat("\nALL PREFLIGHT CHECKS PASSED\n")
RSCRIPT
