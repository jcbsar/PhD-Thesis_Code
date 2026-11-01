
# Chapter 5: MVMR instrument-flow validation
# Educational attainment and type 2 diabetes
#
# Validates the final archived instrument-selection records
# against Appendix Table A6.
#
# This script audits the recorded workflow.
# It does not independently rerun LD clumping.

# 1. Input directory ---------------------------------------

data_dir <- Sys.getenv("MVMR_DATA_DIR", unset = "")

if (!nzchar(data_dir) || !dir.exists(data_dir)) {
  stop("Set MVMR_DATA_DIR to the JOINT_CLUMPED archive.")
}

read_input <- function(filename) {
  read.csv(
    file.path(data_dir, filename),
    stringsAsFactors = FALSE
  )
}

# 2. Load final archived records ---------------------------

candidates <- read_input(
  "01_candidate_union_before_joint_clumping.csv"
)

clumping_audit <- read_input(
  "04_joint_clumping_audit_all_candidates.csv"
)

clumped <- read_input(
  "03_jointly_clumped_instrument_index.csv"
)

harmonisation <- read_input(
  "06_cross_exposure_harmonisation_audit.csv"
)

final_exposures <- read_input(
  "07_MVMR_exposure_dataset_jointly_clumped.csv"
)

reference <- read_input(
  "08_global_instrument_flow_counts.csv"
)

# 3. Validate instrument identities -----------------------

stopifnot(
  nrow(candidates) == 877L,
  !anyDuplicated(candidates$SNP),

  setequal(
    candidates$SNP,
    clumping_audit$SNP
  ),

  !anyDuplicated(clumped$SNP),

  setequal(
    clumped$SNP,
    clumping_audit$SNP[
      clumping_audit$retained_after_joint_clumping
    ]
  ),

  setequal(
    clumped$SNP,
    harmonisation$SNP
  )
)

# 4. Validate the final harmonised SNP set ----------------

final_snps <- unique(final_exposures$SNP)

expected_snps <- harmonisation$SNP[
  harmonisation$retained_after_cross_exposure_harmonisation
]

stopifnot(
  length(final_snps) == 503L,
  setequal(final_snps, expected_snps)
)

# Verify one record per SNP for each exposure
pair_counts <- table(
  final_exposures$SNP,
  final_exposures$exposure
)

stopifnot(
  all(c(
    "Educational attainment",
    "Type 2 diabetes"
  ) %in% colnames(pair_counts)),
  all(pair_counts == 1L)
)

# Confirm alleles agree between exposures
ea <- subset(
  final_exposures,
  exposure == "Educational attainment"
)

t2d <- subset(
  final_exposures,
  exposure == "Type 2 diabetes"
)

ea <- ea[order(ea$SNP), ]
t2d <- t2d[order(t2d$SNP), ]

stopifnot(
  identical(ea$SNP, t2d$SNP),

  identical(
    ea$effect_allele.exposure,
    t2d$effect_allele.exposure
  ),

  identical(
    ea$other_allele.exposure,
    t2d$other_allele.exposure
  )
)

# 5. Recalculate all Appendix A6 counts -------------------

observed <- c(
  "Exposure-specific educational-attainment instruments" =
    sum(grepl(
      "Educational attainment",
      candidates$selected_for,
      fixed = TRUE
    )),

  "Exposure-specific T2D instruments" =
    sum(grepl(
      "Type 2 diabetes",
      candidates$selected_for,
      fixed = TRUE
    )),

  "Unique variants in combined candidate union" =
    nrow(candidates),

  "Candidate variants selected for both exposures" =
    sum(candidates$n_selecting_exposures == 2),

  "Candidate variants present in EUR LD reference" =
    sum(clumping_audit$present_in_EUR_LD_reference),

  "Variants retained after joint LD clumping" =
    nrow(clumped),

  "Jointly clumped variants valid in educational-attainment GWAS" =
    sum(harmonisation$education_row_valid),

  "Jointly clumped variants valid in T2D GWAS" =
    sum(harmonisation$t2d_row_valid),

  "Jointly clumped variants valid in both exposure GWAS" =
    sum(harmonisation$retained_in_both_exposure_GWAS),

  "Variants retained after cross-exposure harmonisation" =
    length(final_snps)
)

# Match recorded counts by stage name
stopifnot(setequal(reference$stage, names(observed)))

validation <- data.frame(
  stage = reference$stage,
  archived = reference$n_snps,
  recalculated = unname(observed[reference$stage]),
  row.names = NULL
)

validation$match <- (
  validation$archived == validation$recalculated
)

print(validation, row.names = FALSE)

stopifnot(all(validation$match))

# 6. Validate exclusion categories ------------------------

stopifnot(
  sum(
    !clumping_audit$retained_after_joint_clumping
  ) == 355L,

  sum(
    !harmonisation$retained_in_both_exposure_GWAS
  ) == 3L,

  sum(
    harmonisation$retained_in_both_exposure_GWAS &
    !harmonisation$retained_after_cross_exposure_harmonisation
  ) == 16L
)

cat("\nAll instrument-flow validation checks passed.\n")
cat("Final unique SNPs:", length(final_snps), "\n")
cat("Total exposure records:", nrow(final_exposures), "\n")
cat("Appendix A6: PASS\n")
