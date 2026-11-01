
#!/usr/bin/env Rscript

# ============================================================
# Chapter 3: FE CFA model comparisons
# ============================================================
#
# Reproduces Appendix Table A14:
#
# 1. One-factor FE/brain model
# 2. Two-factor FE/brain model (retained)
# 3. Three-factor FE/striatal-thalamic/limbic model
#
# All CFA comparisons exclude brainstem.
#
# The three-factor model fixes the F_HS residual variance
# to zero following a marginal Heywood case in the
# unconstrained model.
#
# The three-factor specification was reconstructed and
# independently verified against Appendix Table A14.
# Its exact original script was not recovered.
#
# Input:
# FE_subtypes_brain_ldsc_July2026.rds
#
# Usage:
# Rscript 06_compare_FE_CFA_models.R <path_to_ldsc_rds>
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1L) {
  stop(
    "Usage: Rscript 06_compare_FE_CFA_models.R ",
    "<path_to_FE_subtypes_brain_ldsc_July2026.rds>"
  )
}

if (!requireNamespace("GenomicSEM", quietly = TRUE)) {
  stop("GenomicSEM is required.")
}

# ------------------------------------------------------------
# 1. Load and validate input
# ------------------------------------------------------------

input_file <- normalizePath(
  args[[1]],
  mustWork = TRUE
)

fe_ldsc <- readRDS(input_file)

if (is.null(fe_ldsc$S) || is.null(fe_ldsc$V)) {
  stop("Input must contain LDSC S and V matrices.")
}

expected_traits <- c(
  "F_HS",
  "F_Neg",
  "F_Oth",
  "accumbens",
  "amygdala",
  "brainstem",
  "caudate",
  "pallidum",
  "putamen",
  "thalamus"
)

observed_traits <- colnames(fe_ldsc$S)

if (
  length(observed_traits) != length(expected_traits) ||
  !setequal(observed_traits, expected_traits)
) {
  stop("Unexpected phenotype set in FE LDSC input.")
}

cat("LDSC input:", input_file, "\n")
cat(
  "Phenotypes:",
  paste(observed_traits, collapse = ", "),
  "\n"
)

# ------------------------------------------------------------
# 2. Define CFA models
# ------------------------------------------------------------

# One-factor: all FE and retained brain indicators
# load onto a common factor.

model_FE_1F <- '
  Common =~ F_HS + F_Neg + F_Oth +
            accumbens + amygdala + caudate +
            pallidum + putamen + thalamus
'

# Two-factor: retained FE model.
# FE subtypes and six brain regions form separate,
# correlated latent factors.

model_FE_2F <- '
  FE =~ F_HS + F_Neg + F_Oth

  Brain =~ accumbens + amygdala + caudate +
           pallidum + putamen + thalamus

  FE ~~ Brain
'

# Three-factor: FE, striatal/thalamic and limbic.
#
# In the unconstrained model, the F_HS residual
# variance was marginally negative (-0.00277).
#
# Fixing F_HS residual variance to zero resolves
# that boundary estimate and reproduces A14.

model_FE_3F <- '
  FE =~ F_HS + F_Neg + F_Oth

  StrThal =~ caudate + pallidum +
             putamen + thalamus

  Limbic =~ amygdala + accumbens

  F_HS ~~ 0*F_HS
'

# ------------------------------------------------------------
# 3. Fit models
# ------------------------------------------------------------

fit_cfa <- function(model) {
  GenomicSEM::usermodel(
    covstruc = fe_ldsc,
    model = model,
    estimation = "DWLS",
    CFIcalc = TRUE,
    std.lv = TRUE,
    imp_cov = FALSE
  )
}

cat("\n===== FITTING ONE-FACTOR MODEL =====\n")
fit_1F <- fit_cfa(model_FE_1F)

cat("\n===== FITTING TWO-FACTOR MODEL =====\n")
fit_2F <- fit_cfa(model_FE_2F)

cat("\n===== FITTING THREE-FACTOR MODEL =====\n")
fit_3F <- fit_cfa(model_FE_3F)

fits <- list(
  One_factor = fit_1F,
  Two_factor_retained = fit_2F,
  Three_factor_FHS_boundary = fit_3F
)

if (
  any(vapply(
    fits,
    function(x) is.null(x$modelfit),
    logical(1)
  ))
) {
  stop("At least one CFA model failed to return fit statistics.")
}

# ------------------------------------------------------------
# 4. Assemble model comparison
# ------------------------------------------------------------

comparison <- do.call(
  rbind,
  lapply(names(fits), function(name) {

    data.frame(
      Model = name,
      as.data.frame(fits[[name]]$modelfit),
      row.names = NULL,
      check.names = FALSE
    )
  })
)

rownames(comparison) <- NULL

cat("\n===== APPENDIX TABLE A14: MODEL FIT =====\n")
print(comparison, row.names = FALSE)

# ------------------------------------------------------------
# 5. Validate against independently reproduced results
# ------------------------------------------------------------

reference <- data.frame(
  Model = names(fits),
  chisq = c(
    96.62329,
    70.30092,
    67.88547
  ),
  df = c(27, 26, 25),
  AIC = c(
    132.6233,
    108.3009,
    107.8855
  ),
  CFI = c(
    0.7851359,
    0.8632832,
    0.8676514
  ),
  SRMR = c(
    0.1794701,
    0.07659914,
    0.07363168
  )
)

tolerances <- c(
  chisq = 0.02,
  df = 0,
  AIC = 0.02,
  CFI = 0.002,
  SRMR = 0.002
)

for (metric in names(tolerances)) {

  observed <- comparison[[metric]]
  expected <- reference[[metric]]

  if (
    any(!is.finite(observed)) ||
    any(abs(observed - expected) > tolerances[[metric]])
  ) {
    stop(
      "Reproduction check failed for: ",
      metric
    )
  }
}

cat("\nAppendix Table A14 fit checks: PASSED\n")

# ------------------------------------------------------------
# 6. Check the three-factor boundary constraint
# ------------------------------------------------------------

three_factor_results <- fit_3F$results

boundary_row <- subset(
  three_factor_results,
  lhs == "F_HS" &
    op == "~~" &
    rhs == "F_HS"
)

if (
  nrow(boundary_row) != 1L ||
  !is.finite(boundary_row$Unstand_Est) ||
  abs(boundary_row$Unstand_Est) > 1e-6
) {
  stop("F_HS residual variance is not fixed to zero.")
}

cat("\nThree-factor F_HS boundary check: PASSED\n")

# ------------------------------------------------------------
# 7. Prepare output directory
# ------------------------------------------------------------

script_arg <- grep(
  "^--file=",
  commandArgs(trailingOnly = FALSE),
  value = TRUE
)

if (length(script_arg) != 1L) {
  stop("Unable to determine script location.")
}

script_file <- sub(
  "^--file=",
  "",
  script_arg[[1]]
)

section_dir <- normalizePath(
  file.path(dirname(script_file), ".."),
  mustWork = TRUE
)

output_dir <- file.path(
  section_dir,
  "outputs",
  "FE_CFA_comparison"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 8. Save results and model definitions
# ------------------------------------------------------------

write.csv(
  comparison,
  file.path(
    output_dir,
    "FE_CFA_model_comparison.csv"
  ),
  row.names = FALSE
)

saveRDS(
  fits,
  file.path(
    output_dir,
    "FE_CFA_comparison_models.rds"
  )
)

for (name in names(fits)) {
  write.csv(
    fits[[name]]$results,
    file.path(
      output_dir,
      paste0(name, "_parameters.csv")
    ),
    row.names = FALSE
  )
}

model_syntax <- c(
  "ONE-FACTOR MODEL",
  model_FE_1F,
  "",
  "TWO-FACTOR RETAINED MODEL",
  model_FE_2F,
  "",
  "THREE-FACTOR BOUNDARY-CONSTRAINED MODEL",
  model_FE_3F
)

writeLines(
  model_syntax,
  file.path(
    output_dir,
    "FE_CFA_model_syntax.txt"
  )
)

# ------------------------------------------------------------
# 9. Save provenance and session information
# ------------------------------------------------------------

provenance <- c(
  paste("Input:", input_file),
  paste(
    "Input MD5:",
    unname(tools::md5sum(input_file))
  ),
  "",
  "Appendix: Chapter 3, Table A14",
  "",
  "One-factor:",
  "Nine indicators; brainstem excluded.",
  "",
  "Two-factor:",
  "FE and six-region brain factor;",
  "brainstem excluded.",
  "",
  "Three-factor:",
  "FE, striatal/thalamic and limbic factors;",
  "brainstem excluded.",
  "F_HS residual variance fixed to zero.",
  "",
  "The unconstrained three-factor model produced",
  "a marginal negative F_HS residual variance.",
  "The boundary-constrained specification",
  "reproduces the A14 fit statistics.",
  "",
  "The exact historical three-factor source script",
  "was not recovered. This specification was",
  "reconstructed and independently verified.",
  "",
  "Estimation: DWLS",
  "Latent variances fixed to one.",
  "GenomicSEM usermodel(); CFIcalc = TRUE;",
  "imp_cov = FALSE."
)

writeLines(
  provenance,
  file.path(
    output_dir,
    "FE_CFA_comparison_provenance.txt"
  )
)

writeLines(
  capture.output(sessionInfo()),
  file.path(
    output_dir,
    "FE_CFA_comparison_session_info.txt"
  )
)

cat("\nResults saved to:", output_dir, "\n")
cat("FE CFA model comparison completed.\n")
