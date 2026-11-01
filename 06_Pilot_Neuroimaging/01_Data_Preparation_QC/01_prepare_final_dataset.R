# Chapter 6 — final ABHGP imaging dataset preparation (2026)
# Restricted participant-level data MUST remain outside Git.
# Supply a locally authorised, completed-study master RDS via ABHGP_INPUT_RDS.
# This is a reconstructed script; numerical equivalence has not yet been established.

suppressPackageStartupMessages({
  library(dplyr)
})

input_path <- Sys.getenv("ABHGP_INPUT_RDS")
output_dir <- Sys.getenv("ABHGP_PRIVATE_OUTPUT_DIR")
if (!nzchar(input_path) || !file.exists(input_path)) {
  stop("Set ABHGP_INPUT_RDS to the authorised local final-study RDS file")
}
if (!nzchar(output_dir)) {
  stop("Set ABHGP_PRIVATE_OUTPUT_DIR to a private location outside the repository")
}
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

dat <- readRDS(input_path)
if (!is.data.frame(dat)) stop("Input must be a data frame")
required <- c("SubjID", "diagnosis", "age", "sex", "ICV", "Euler_LH", "Euler_RH",
              "Lthal", "Rthal", "Lcaud", "Rcaud", "Lput", "Rput", "Lpal", "Rpal",
              "Lhippo", "Rhippo", "Lamyg", "Ramyg", "Laccumb", "Raccumb",
              "LLatVent", "RLatVent", "Whole_brainstem")
missing_cols <- setdiff(required, names(dat))
if (length(missing_cols)) stop("Required fields missing: ", paste(missing_cols, collapse = ", "))
if (anyDuplicated(dat$SubjID)) stop("Duplicate participant keys detected")
if (!all(na.omit(unique(dat$diagnosis)) %in% c("Control", "GGE", "FE", "FDS"))) {
  stop("Unexpected diagnosis values")
}

# Strict conversion prevents silent coercion of malformed numeric data.
convert_numeric <- function(x, field) {
  if (is.numeric(x)) return(x)
  original <- as.character(x)
  out <- suppressWarnings(as.numeric(original))
  bad <- !is.na(original) & trimws(original) != "" & is.na(out)
  if (any(bad)) stop("Non-numeric entries in ", field, " (", sum(bad), " rows)")
  out
}
raw_volume <- c("Lthal", "Rthal", "Lcaud", "Rcaud", "Lput", "Rput", "Lpal", "Rpal",
                "Lhippo", "Rhippo", "Lamyg", "Ramyg", "Laccumb", "Raccumb",
                "LLatVent", "RLatVent", "Whole_brainstem", "ICV")
# Not every prefixed field is necessarily a measurement; select recognised types only.
measure_cols <- unique(c(raw_volume, "age", "Euler_LH", "Euler_RH",
                         grep("_thickavg$|_surfavg$|^Left-|^Right-", names(dat), value = TRUE)))
for (nm in measure_cols) dat[[nm]] <- convert_numeric(dat[[nm]], nm)

dat <- dat %>% mutate(
  diagnosis4 = factor(diagnosis, levels = c("Control", "GGE", "FE", "FDS")),
  diagnosis3 = factor(case_when(
    diagnosis == "Control" ~ "Control",
    diagnosis %in% c("GGE", "FE") ~ "Epilepsy",
    diagnosis == "FDS" ~ "FDS",
    TRUE ~ NA_character_
  ), levels = c("Control", "Epilepsy", "FDS")),
  Thalamus = Lthal + Rthal,
  Caudate = Lcaud + Rcaud,
  Putamen = Lput + Rput,
  Pallidum = Lpal + Rpal,
  Hippocampus = Lhippo + Rhippo,
  Amygdala = Lamyg + Ramyg,
  Accumbens = Laccumb + Raccumb,
  LateralVentricles = LLatVent + RLatVent,
  Brainstem = Whole_brainstem,
  Euler_Mean = (Euler_LH + Euler_RH) / 2
)

# Explicitly validate final sample. Do not silently filter to desired counts.
expected <- c(Control = 61L, GGE = 21L, FE = 49L, FDS = 19L)
observed <- table(dat$diagnosis4)
if (!identical(as.integer(observed), as.integer(expected))) {
  stop("Input is not the reported final analytical cohort; observed counts: ",
       paste(names(observed), observed, collapse = ", "))
}

# Nuclei/subfields: use an explicit local dictionary rather than a broad prefix
# match that can unintentionally include unrelated FreeSurfer measurements.
# The dictionary must be checked against the final 25/18/9 outcome sets.
# See README; subregional validation will be incorporated in the next module.

saveRDS(dat, file.path(output_dir, "final_prepared_private.rds"))
message("Private prepared dataset written to: ", output_dir)
