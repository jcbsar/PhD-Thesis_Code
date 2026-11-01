
# ============================================================
# PILOT NEUROIMAGING STUDY
# 01 - FINAL DATA PREPARATION AND QUALITY CONTROL
#
# This script prepares the completed-study analytical
# dataset used for primary subcortical, subregional,
# cortical thickness and surface-area analyses.
#
# Participant-level data, QC decisions and mappings
# are restricted and are not included in GitHub.
#
# Reconstructed from the original 2026 workflow.
# Numerical equivalence has not yet been verified.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
})

# ------------------------------------------------------------
# 1. PRIVATE INPUTS AND OUTPUTS
# ------------------------------------------------------------

input_path <- Sys.getenv("ABHGP_INPUT_RDS")
output_dir <- Sys.getenv("ABHGP_PRIVATE_OUTPUT_DIR")
dictionary_path <- Sys.getenv("ABHGP_SUBREGION_DICTIONARY")
qc_path <- Sys.getenv("ABHGP_QC_MANIFEST_RDS")

private_files <- c(
  input_path,
  dictionary_path,
  qc_path
)

if (any(!nzchar(private_files)) ||
    !all(file.exists(private_files))) {
  stop(
    "Set ABHGP_INPUT_RDS, ABHGP_SUBREGION_DICTIONARY ",
    "and ABHGP_QC_MANIFEST_RDS to authorised local files."
  )
}

if (!nzchar(output_dir)) {
  stop("Set ABHGP_PRIVATE_OUTPUT_DIR.")
}


# Resolve the actual Git repository root, independently
# of the current R working directory.

git_executable <- Sys.which("git")

if (!nzchar(git_executable)) {
  stop("Git is required to verify the private output location.")
}

repo_root_raw <- suppressWarnings(
  system2(
    git_executable,
    c("rev-parse", "--show-toplevel"),
    stdout = TRUE,
    stderr = FALSE
  )
)

if (
  length(repo_root_raw) != 1L ||
  !is.null(attr(repo_root_raw, "status"))
) {
  stop(
    "Cannot determine Git repository root. ",
    "Run from within the repository."
  )
}

repo_root <- normalizePath(
  trimws(repo_root_raw),
  winslash = "/",
  mustWork = TRUE
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

output_dir <- normalizePath(
  output_dir,
  winslash = "/",
  mustWork = TRUE
)

# Protect against both the repository itself and
# directories nested inside it.

is_within <- function(path, parent) {
  path <- tolower(
    sub("/+$", "", normalizePath(
      path, winslash = "/", mustWork = TRUE
    ))
  )
  parent <- tolower(
    sub("/+$", "", normalizePath(
      parent, winslash = "/", mustWork = TRUE
    ))
  )

  identical(path, parent) ||
    startsWith(path, paste0(parent, "/"))
}

if (is_within(output_dir, repo_root)) {
  stop(
    "Private outputs must not be written ",
    "inside the Git repository."
  )
}

# Do not overwrite the original input files.

private_input_paths <- normalizePath(
  private_files,
  winslash = "/",
  mustWork = TRUE
)

if (any(vapply(
  private_input_paths,
  is_within,
  logical(1),
  parent = repo_root
))) {
  stop(
    "Restricted input data and QC files must be ",
    "stored outside the Git repository."
  )
}

# ------------------------------------------------------------
# LOAD RESTRICTED SOURCE DATA
# ------------------------------------------------------------

dat <- readRDS(input_path)

if (!is.data.frame(dat)) {
  stop("Input must be a data frame.")
}

# The original master workbook contains incidental
# whitespace in some anatomical column headers.

names(dat) <- trimws(names(dat))

if (anyDuplicated(names(dat))) {
  stop(
    "Duplicate column names after whitespace normalisation."
  )
}

# ------------------------------------------------------------
# 2. VALIDATE SOURCE VARIABLES
# ------------------------------------------------------------

primary_raw <- c(
  "Lthal", "Rthal",
  "Lcaud", "Rcaud",
  "Lput", "Rput",
  "Lpal", "Rpal",
  "Lhippo", "Rhippo",
  "Lamyg", "Ramyg",
  "Laccumb", "Raccumb",
  "LLatVent", "RLatVent",
  "Whole_brainstem"
)

required <- c(
  "SubjID", "diagnosis",
  "age", "sex", "ICV",
  "Euler_LH", "Euler_RH",
  primary_raw
)

missing <- setdiff(required, names(dat))

if (length(missing) > 0) {
  stop(
    "Missing source variables: ",
    paste(missing, collapse = ", ")
  )
}

if (anyNA(dat$SubjID) ||
    anyDuplicated(dat$SubjID)) {
  stop("Missing or duplicated participant keys.")
}

# ------------------------------------------------------------
# 3. VALIDATE FINAL COHORT
# ------------------------------------------------------------

dat$diagnosis4 <- factor(
  as.character(dat$diagnosis),
  levels = c("Control", "GGE", "FE", "FDS")
)

dat$diagnosis3 <- factor(
  case_when(
    dat$diagnosis == "Control" ~ "Control",
    dat$diagnosis %in% c("GGE", "FE") ~ "Epilepsy",
    dat$diagnosis == "FDS" ~ "FDS",
    TRUE ~ NA_character_
  ),
  levels = c("Control", "Epilepsy", "FDS")
)

expected_counts <- c(
  Control = 61L,
  GGE = 21L,
  FE = 49L,
  FDS = 19L
)

observed_counts <- table(dat$diagnosis4)

stopifnot(
  nrow(dat) == 150L,
  !anyNA(dat$diagnosis4),
  identical(
    as.integer(observed_counts),
    as.integer(expected_counts)
  )
)

dat$sex <- factor(dat$sex)

# ------------------------------------------------------------
# 4. LOAD PRIVATE SUBREGION DICTIONARY
# ------------------------------------------------------------

# Dictionary columns:
#
# family        - Thalamus, Hippocampus or Amygdala
# output_name   - Final bilateral variable name
# left_column   - Original left measurement
# right_column  - Original right measurement
#
# This dictionary is maintained privately.

dictionary <- read.csv(
  dictionary_path,
  stringsAsFactors = FALSE
)

dictionary_fields <- c(
  "family",
  "output_name",
  "left_column",
  "right_column"
)

stopifnot(
  all(dictionary_fields %in% names(dictionary))
)

expected_families <- c(
  Thalamus = 25L,
  Hippocampus = 18L,
  Amygdala = 9L
)

observed_families <- table(
  factor(
    dictionary$family,
    levels = names(expected_families)
  )
)

if (
  nrow(dictionary) != 52L ||
  anyNA(dictionary[dictionary_fields]) ||
  any(!nzchar(trimws(
    as.matrix(dictionary[dictionary_fields])
  ))) ||
  !identical(
    as.integer(observed_families),
    as.integer(expected_families)
  ) ||
  anyDuplicated(dictionary$output_name)
) {
  stop("Invalid subregion dictionary.")
}

prefixes <- c(
  Thalamus = "Thal_",
  Hippocampus = "Hip_",
  Amygdala = "Amyg_"
)

correct_prefix <- mapply(
  startsWith,
  dictionary$output_name,
  unname(prefixes[dictionary$family])
)

if (!all(correct_prefix)) {
  stop("Subregion output naming mismatch.")
}

source_cols <- unique(c(
  dictionary$left_column,
  dictionary$right_column
))

if (!all(source_cols %in% names(dat))) {
  stop("Subregion source columns are missing.")
}

if (any(
  dictionary$left_column == dictionary$right_column
)) {
  stop("Left and right columns must differ.")
}

# Check against the original FreeSurfer naming
# conventions used in the August 2026 analysis.

expected_source_prefix <- function(family, side) {
  if (family == "Thalamus") {
    return(if (side == "L") "^Left-" else "^Right-")
  }
  if (side == "L") "^L_" else "^R_"
}

for (i in seq_len(nrow(dictionary))) {
  fam <- dictionary$family[i]

  if (!grepl(
    expected_source_prefix(fam, "L"),
    dictionary$left_column[i]
  ) || !grepl(
    expected_source_prefix(fam, "R"),
    dictionary$right_column[i]
  )) {
    stop("Unexpected source naming at dictionary row ", i)
  }
}

# ------------------------------------------------------------
# 5. CONVERT MEASUREMENTS TO NUMERIC
# ------------------------------------------------------------

convert_numeric <- function(x, field) {
  if (is.numeric(x)) return(x)

  original <- as.character(x)
  converted <- suppressWarnings(
    as.numeric(original)
  )

  invalid <- (
    !is.na(original) &
    nzchar(trimws(original)) &
    is.na(converted)
  )

  if (any(invalid)) {
    stop("Non-numeric values in ", field)
  }

  converted
}

thickness_cols <- grep(
  "^[LR]_.*_thickavg$",
  names(dat),
  value = TRUE
)

surface_cols <- grep(
  "^[LR]_.*_surfavg$",
  names(dat),
  value = TRUE
)

for (hemi in c("L", "R")) {
  if (
    sum(startsWith(thickness_cols, paste0(hemi, "_"))) != 34L ||
    sum(startsWith(surface_cols, paste0(hemi, "_"))) != 34L
  ) {
    stop("Expected 34 cortical regions per hemisphere.")
  }
}

numeric_cols <- unique(c(
  primary_raw,
  "age", "ICV",
  "Euler_LH", "Euler_RH",
  source_cols,
  thickness_cols,
  surface_cols
))

for (nm in numeric_cols) {
  dat[[nm]] <- convert_numeric(
    dat[[nm]],
    nm
  )
}

# ------------------------------------------------------------
# 6. APPLY PRIVATE QUALITY-CONTROL DECISIONS
# ------------------------------------------------------------

# Private QC manifest structure:
#
# region_missing:
#   data.frame(SubjID, variable)
#
# right_euler_exclusion_id:
#   the participant key for the original sensitivity
#   exclusion, stored only in the private manifest.
#
# The manifest must reproduce the original 2026 QC
# decisions. It is not generated automatically.

qc <- readRDS(qc_path)

if (
  !is.list(qc) ||
  !all(c(
    "region_missing",
    "right_euler_exclusion_id"
  ) %in% names(qc))
) {
  stop("Invalid private QC manifest.")
}

qc_regions <- qc$region_missing

if (
  !is.data.frame(qc_regions) ||
  !all(c("SubjID", "variable") %in% names(qc_regions))
) {
  stop("Invalid regional QC entries.")
}

if (
  anyNA(qc_regions$SubjID) ||
  anyNA(qc_regions$variable) ||
  !all(qc_regions$SubjID %in% dat$SubjID) ||
  !all(qc_regions$variable %in%
       c(thickness_cols, surface_cols))
) {
  stop("Regional QC entries do not match source data.")
}

           
# ------------------------------------------------------------
# 6A. VERIFY ORIGINAL CORTICAL QC DECISIONS
# ------------------------------------------------------------

# These conditions were recovered from the original
# August 2026 analytical script.
#
# Participant identifiers remain in the private QC
# manifest and are not included in this repository.

qc_expected <- c(
  "L_entorhinal_thickavg",
  "L_entorhinal_surfavg",
  "L_temporalpole_surfavg"
)

# Exactly three region-specific corrections were
# documented in the original analytical source.

if (
  nrow(qc_regions) != 3L ||
  anyDuplicated(qc_regions$variable) ||
  anyDuplicated(qc_regions[c("SubjID", "variable")]) ||
  !setequal(qc_regions$variable, qc_expected)
) {
  stop(
    "Private QC manifest does not match the ",
    "three original cortical QC decisions."
  )
}

# The three corrections concern the same participant.

if (length(unique(qc_regions$SubjID)) != 1L) {
  stop("Expected one participant for regional QC corrections.")
}

# Check the source measurements before replacing them.

for (i in seq_len(nrow(qc_regions))) {

  participant <- qc_regions$SubjID[i]
  variable <- qc_regions$variable[i]

  value <- dat[
    dat$SubjID == participant,
    variable
  ][[1]]

  if (length(value) != 1L || !is.finite(value)) {
    stop("Missing or invalid original QC measurement.")
  }

  if (
    variable %in% c(
      "L_entorhinal_thickavg",
      "L_entorhinal_surfavg"
    ) && value != 0
  ) {
    stop(
      "Original entorhinal QC rule requires ",
      "a zero measurement."
    )
  }

  if (
    variable == "L_temporalpole_surfavg" &&
    value > 10
  ) {
    stop(
      "Original temporal pole QC rule requires ",
      "a measurement <= 10 mm2."
    )
  }
}

message(
  "Original regional cortical QC conditions verified."
)

# The existing region-specific NA assignment follows.

# Apply region-specific missing-value corrections.

for (i in seq_len(nrow(qc_regions))) {
  participant <- qc_regions$SubjID[i]
  variable <- qc_regions$variable[i]

  dat[
    dat$SubjID == participant,
    variable
  ] <- NA_real_
}

# Define the original right-hemisphere
# sensitivity exclusion without embedding IDs.

sensitivity_id <- qc$right_euler_exclusion_id

if (
  length(sensitivity_id) != 1L ||
  is.na(sensitivity_id) ||
  sum(dat$SubjID == sensitivity_id) != 1L
) {
  stop("Invalid sensitivity exclusion key.")
}

dat$exclude_right_euler_sensitivity <- (
  dat$SubjID == sensitivity_id
)

stopifnot(
  is.logical(dat$exclude_right_euler_sensitivity),
  sum(dat$exclude_right_euler_sensitivity) == 1L
)

# ------------------------------------------------------------
# 7. PRIMARY BILATERAL VOLUMES
# ------------------------------------------------------------

dat <- dat %>%
  mutate(
    Thalamus = Lthal + Rthal,
    Caudate = Lcaud + Rcaud,
    Putamen = Lput + Rput,
    Pallidum = Lpal + Rpal,
    Hippocampus = Lhippo + Rhippo,
    Amygdala = Lamyg + Ramyg,
    Accumbens = Laccumb + Raccumb,
    LateralVentricles = LLatVent + RLatVent,
    Brainstem = Whole_brainstem
  )

# ------------------------------------------------------------
# 8. BILATERAL SUBREGIONAL VOLUMES
# ------------------------------------------------------------


# ------------------------------------------------------------
# 8A. VERIFY ANATOMICAL MAPPINGS
# ------------------------------------------------------------

# Verify that left and right source columns refer to
# the same anatomical structure and that output names
# reproduce the August 2026 naming convention.

source_suffix <- function(x, family, side) {

  prefix <- if (family == "Thalamus") {
    if (side == "left") "Left-" else "Right-"
  } else {
    if (side == "left") "L_" else "R_"
  }

  if (!startsWith(x, prefix)) {
    stop("Unexpected source prefix: ", x)
  }

  substring(x, nchar(prefix) + 1L)
}

expected_outputs <- character(nrow(dictionary))

for (i in seq_len(nrow(dictionary))) {

  family <- dictionary$family[i]

  left_name <- source_suffix(
    dictionary$left_column[i],
    family,
    "left"
  )

  right_name <- source_suffix(
    dictionary$right_column[i],
    family,
    "right"
  )

  if (!identical(left_name, right_name)) {
    stop(
      "Left/right anatomical mismatch at row ", i,
      ": ", left_name, " vs ", right_name
    )
  }

  expected_outputs[i] <- paste0(
    prefixes[[family]],
    gsub("[^A-Za-z0-9]", "_", left_name)
  )
}

if (!identical(
  dictionary$output_name,
  expected_outputs
)) {
  stop(
    "Dictionary output names differ from the ",
    "original August 2026 naming convention."
  )
}
           
# ------------------------------------------------------------
# 8B. VERIFY HISTORICAL SUBREGIONAL OUTCOME SET
# ------------------------------------------------------------

# Exact ROI names recovered from the archived
# August 2026 full-results CSV files.
#
# This reference contains anatomical names only.
# It does not include participant-level information.

archived_rois <- list(

  Thalamus = c(
    "Thal_AV", "Thal_CL", "Thal_CM", "Thal_CeM",
    "Thal_LD", "Thal_LGN", "Thal_LP", "Thal_L_Sg",
    "Thal_MDl", "Thal_MDm", "Thal_MGN",
    "Thal_MV_Re_", "Thal_Pc", "Thal_Pf",
    "Thal_Pt", "Thal_PuA", "Thal_PuI",
    "Thal_PuL", "Thal_PuM", "Thal_VA",
    "Thal_VAmc", "Thal_VLa", "Thal_VLp",
    "Thal_VM", "Thal_VPL"
  ),

  Hippocampus = c(
    "Hip_CA1_body", "Hip_CA1_head",
    "Hip_CA3_body", "Hip_CA3_head",
    "Hip_CA4_body", "Hip_CA4_head",
    "Hip_GC_ML_DG_body", "Hip_GC_ML_DG_head",
    "Hip_HATA", "Hip_Hippocampal_tail",
    "Hip_fimbria",
    "Hip_molecular_layer_HP_body",
    "Hip_molecular_layer_HP_head",
    "Hip_parasubiculum",
    "Hip_presubiculum_body",
    "Hip_presubiculum_head",
    "Hip_subiculum_body",
    "Hip_subiculum_head"
  ),

  Amygdala = c(
    "Amyg_Accessory_Basal_nucleus",
    "Amyg_Anterior_amygdaloid_area_AAA",
    "Amyg_Basal_nucleus",
    "Amyg_Central_nucleus",
    "Amyg_Cortical_nucleus",
    "Amyg_Corticoamygdaloid_transitio",
    "Amyg_Lateral_nucleus",
    "Amyg_Medial_nucleus",
    "Amyg_Paralaminar_nucleus"
  )
)

for (family in names(archived_rois)) {

  observed <- dictionary$output_name[
    dictionary$family == family
  ]

  expected <- archived_rois[[family]]

  if (!setequal(observed, expected)) {
    stop(
      "Anatomical outcomes do not match the ",
      "archived August 2026 results for: ",
      family,
      "\nMissing: ",
      paste(setdiff(expected, observed), collapse = ", "),
      "\nUnexpected: ",
      paste(setdiff(observed, expected), collapse = ", ")
    )
  }
}

message(
  "All 52 subregional outcome names match ",
  "the archived August 2026 result sets."
)


if (
  anyDuplicated(dictionary$left_column) ||
  anyDuplicated(dictionary$right_column) ||
  anyDuplicated(expected_outputs)
) {
  stop(
    "Duplicate source measurements or bilateral ",
    "output names detected."
  )
}

message(
  "Subregional dictionary structural checks passed: ",
  "25 thalamic, 18 hippocampal and 9 amygdala regions."
)

# These checks validate naming, pairing and counts.
# They do not independently establish that the
# anatomical selections match the final thesis.


for (i in seq_len(nrow(dictionary))) {

  output <- dictionary$output_name[i]
  left <- dictionary$left_column[i]
  right <- dictionary$right_column[i]

  if (output %in% names(dat)) {
    stop("Bilateral variable already exists: ", output)
  }

  dat[[output]] <- (
    dat[[left]] + dat[[right]]
  )
}

# Confirm all 52 expected bilateral measures.

stopifnot(
  sum(grepl("^Thal_", names(dat))) == 25L,
  sum(grepl("^Hip_", names(dat))) == 18L,
  sum(grepl("^Amyg_", names(dat))) == 9L
)

# ------------------------------------------------------------
# 9. FINAL CORTICAL VARIABLES
# ------------------------------------------------------------

dat$Euler_Mean <- (
  dat$Euler_LH + dat$Euler_RH
) / 2

# Cortical thickness and surface area remain
# separate by hemisphere.

# Thickness models use age, sex and hemisphere-
# specific Euler number.
#
# Surface-area models additionally include ICV.

# ------------------------------------------------------------
# 10. SAVE RESTRICTED PREPARED DATASET
# ------------------------------------------------------------

output_file <- file.path(
  output_dir,
  "final_prepared_private.rds"
)

saveRDS(
  dat,
  output_file
)

message(
  "Final analytical dataset prepared. ",
  "Restricted output saved locally."
)

# ============================================================
# END
# ============================================================
