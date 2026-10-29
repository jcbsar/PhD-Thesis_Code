
# Chapter 5: Conditional instrument-strength sensitivity
# Educational attainment and type 2 diabetes
#
# R 4.4.2 | MVMR 0.4.6
#
# Uses the final jointly clumped, harmonised MVMR datasets.
# Validates against thesis Appendix Tables A9 and A10.

library(MVMR)

# 1. Input directory ---------------------------------------

data_dir <- Sys.getenv("MVMR_DATA_DIR", unset = "")

if (!nzchar(data_dir) || !dir.exists(data_dir)) {
  stop("Set MVMR_DATA_DIR to the JOINT_CLUMPED archive.")
}

outcomes <- c(
  All_Epilepsy = "All epilepsy",
  Focal_Epilepsy = "FE",
  GGE_Epilepsy = "GGE",
  FocalHS_Epilepsy = "Focal HS epilepsy",
  FocalLesionNegative_Epilepsy = "Focal lesion negative epilepsy",
  JME_Epilepsy = "JME"
)

rho_values <- c(
  -0.30, -0.20, -0.10, -0.05,
   0.00,
   0.05,  0.10,  0.20,  0.30
)

# 2. Calculate conditional F-statistics -------------------

calculate_F <- function(formatted, exposure_se, rho) {

  if (rho == 0) {
    covariance <- 0
  } else {
    covariance <- lapply(
      seq_len(nrow(exposure_se)),
      function(i) {
        s1 <- exposure_se[i, 1]
        s2 <- exposure_se[i, 2]

        matrix(
          c(
            s1^2, rho * s1 * s2,
            rho * s1 * s2, s2^2
          ),
          nrow = 2
        )
      }
    )
  }

  result <- NULL

  invisible(capture.output({
    result <- suppressWarnings(
      MVMR::strength_mvmr(
        r_input = formatted,
        gencov = covariance
      )
    )
  }))

  as.numeric(unlist(
    result[1, c("exposure1", "exposure2")],
    use.names = FALSE
  ))
}

# 3. Run six outcomes and nine sensitivity settings -------

results <- list()
k <- 1L

for (outcome in names(outcomes)) {

  dat <- readRDS(file.path(
    data_dir,
    paste0("harmonised_MVMR_data_", outcome, ".rds")
  ))

  bx <- as.matrix(dat$exposure_beta)
  bxse <- as.matrix(dat$exposure_se)

  stopifnot(
    ncol(bx) == 2L,
    identical(dim(bx), dim(bxse)),
    all(is.finite(bx)),
    all(is.finite(bxse)),
    all(bxse > 0)
  )

  snps <- rownames(bx)

  if (is.null(snps)) {
    snps <- as.character(seq_len(nrow(bx)))
  }

  formatted <- MVMR::format_mvmr(
    BXGs = bx,
    BYG = as.numeric(dat$outcome_beta),
    seBXGs = bxse,
    seBYG = as.numeric(dat$outcome_se),
    RSID = snps
  )

  for (rho in rho_values) {

    F_values <- calculate_F(formatted, bxse, rho)

    results[[k]] <- data.frame(
      outcome_label = outcome,
      outcome_display = unname(outcomes[[outcome]]),
      exposure = c("Educational attainment", "Type 2 diabetes"),
      nsnp = nrow(bx),
      rho = rho,
      conditional_F = F_values
    )

    k <- k + 1L
  }
}

results <- do.call(rbind, results)

# 4. Validate zero-covariance F-statistics ----------------

reference_zero <- read.csv(file.path(
  data_dir, "17_conditional_F_lambda_0_validation.csv"
))

zero_results <- subset(results, rho == 0)

for (i in seq_len(nrow(reference_zero))) {

  ref <- reference_zero[i, ]

  calculated <- subset(
    zero_results,
    outcome_label == ref$outcome_label &
      exposure == ref$exposure_label
  )

  stopifnot(
    nrow(calculated) == 1L,
    abs(calculated$conditional_F -
          ref$baseline_conditional_F) < 1e-8
  )
}

cat("Appendix A9 zero-covariance validation: PASS\n")

# 5. Validate nine-point covariance sensitivity ----------

reference_sensitivity <- read.csv(
  file.path(
    data_dir,
    "16_conditional_F_covariance_sensitivity_thesis_table.csv"
  ),
  sep = ";",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

stopifnot(
  nrow(reference_sensitivity) == 12L,
  ncol(reference_sensitivity) == 12L
)

for (i in seq_len(nrow(reference_sensitivity))) {

  ref <- reference_sensitivity[i, ]

  exposure_label <- if (
    ref[[2]] == "T2D"
  ) "Type 2 diabetes" else ref[[2]]

  calculated <- subset(
    results,
    outcome_display == ref[[1]] &
      exposure == exposure_label
  )

  calculated <- calculated[
    match(rho_values, calculated$rho), ,
    drop = FALSE
  ]

  expected_F <- as.numeric(unlist(
    ref[1, 4:12],
    use.names = FALSE
  ))

  stopifnot(
    nrow(calculated) == 9L,
    !anyNA(calculated$rho),
    all(abs(calculated$conditional_F -
              expected_F) <= 0.00051)
  )
}

cat("Appendix A10 covariance sensitivity validation: PASS\n")

# 6. Export results ---------------------------------------

output_dir <- "outputs"
dir.create(output_dir, showWarnings = FALSE)

write.csv(
  results,
  file.path(output_dir, "MVMR_conditional_F_sensitivity.csv"),
  row.names = FALSE
)

cat("\nAll 108 sensitivity estimates validated.\n")
cat("Minimum conditional F:",
    round(min(results$conditional_F), 3), "\n")
cat("All conditional F-statistics > 10:",
    all(results$conditional_F > 10), "\n")
