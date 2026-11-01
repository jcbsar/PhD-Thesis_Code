# Final 2026 cortical models — reconstruction of the original model specifications.
# Restricted dataset and participant identifiers remain outside the repository.
suppressPackageStartupMessages({
  library(dplyr)
  library(emmeans)
})
private_dir <- Sys.getenv("ABHGP_PRIVATE_OUTPUT_DIR")
if (!nzchar(private_dir) || !file.exists(file.path(private_dir, "final_prepared_private.rds"))) {
  stop("Run the preparation script first; set ABHGP_PRIVATE_OUTPUT_DIR")
}
dat <- readRDS(file.path(private_dir, "final_prepared_private.rds"))

spec <- list(
  diagnosis4 = list("GGE vs Control" = c(-1,1,0,0), "FE vs Control" = c(-1,0,1,0),
                    "FDS vs Control" = c(-1,0,0,1)),
  diagnosis3 = list("Epilepsy vs Control" = c(-1,1,0), "FDS vs Control" = c(-1,0,1),
                    "Epilepsy vs FDS" = c(0,1,-1))
)
target_regions <- c("caudalmiddlefrontal", "lateralorbitofrontal", "medialorbitofrontal",
                    "parsopercularis", "parsorbitalis", "parstriangularis",
                    "rostralmiddlefrontal", "superiorfrontal", "frontalpole",
                    "precentral", "postcentral", "paracentral")

run_cortex <- function(data, suffix, group, target_only = FALSE, area = FALSE) {
  stopifnot(group %in% names(spec))
  if (target_only && area) stop("Target ROI set was only specified for thickness")
  rlist <- list()
  for (hemi in c("L", "R")) {
    roi_cols <- grep(paste0("^", hemi, "_.*", suffix, "$"), names(data), value = TRUE)
    if (target_only) roi_cols <- paste0(hemi, "_", target_regions, suffix)
    if (!all(roi_cols %in% names(data))) stop("Missing cortical ROI columns")
    if (length(roi_cols) != if (target_only) 12 else 34) {
      stop("Unexpected number of cortical regions in ", hemi, ": ", length(roi_cols))
    }
    euler <- if (hemi == "L") "Euler_LH" else "Euler_RH"
    covars <- if (area) paste("age + sex + ICV +", euler) else paste("age + sex +", euler)
    for (roi in roi_cols) {
      mod <- lm(stats::as.formula(paste(roi, "~", group, "+", covars)), data = data)
      em <- emmeans(mod, specs = group)
      ans <- as.data.frame(summary(contrast(em, method = spec[[group]]), infer = c(TRUE, TRUE)))
      ans$ROI <- roi
      ans$Hemisphere <- if (hemi == "L") "Left" else "Right"
      rlist[[roi]] <- ans
    }
  }
  bind_rows(rlist) %>%
    group_by(Hemisphere, contrast) %>%
    mutate(q_BH = p.adjust(p.value, method = "BH")) %>%
    ungroup() %>%
    select(Hemisphere, ROI, contrast, estimate, SE, df, lower.CL, upper.CL,
           t.ratio, p.value, q_BH)
}

results <- list(
  targeted_thickness_4group = run_cortex(dat, "_thickavg", "diagnosis4", target_only = TRUE),
  targeted_thickness_3group = run_cortex(dat, "_thickavg", "diagnosis3", target_only = TRUE),
  whole_thickness_4group = run_cortex(dat, "_thickavg", "diagnosis4"),
  whole_thickness_3group = run_cortex(dat, "_thickavg", "diagnosis3"),
  surface_area_4group = run_cortex(dat, "_surfavg", "diagnosis4", area = TRUE),
  surface_area_3group = run_cortex(dat, "_surfavg", "diagnosis3", area = TRUE)
)

# Right-hemisphere sensitivity analyses require an institutionally maintained
# private exclusion flag. No participant key is embedded in the public code.
# Apply identical QC and study eligibility rules to the input before running.
if ("exclude_right_euler_sensitivity" %in% names(dat)) {
  if (!is.logical(dat$exclude_right_euler_sensitivity)) stop("QC flag must be logical")
  keep <- !is.na(dat$exclude_right_euler_sensitivity) & !dat$exclude_right_euler_sensitivity
  sensitivity_dat <- dat[keep, , drop = FALSE]
  if (sum(!keep) != 1L) stop("Expected one excluded sensitivity participant")
  results$targeted_thickness_4group_right_sensitivity <- run_cortex(sensitivity_dat, "_thickavg", "diagnosis4", TRUE) %>% filter(Hemisphere == "Right")
  results$targeted_thickness_3group_right_sensitivity <- run_cortex(sensitivity_dat, "_thickavg", "diagnosis3", TRUE) %>% filter(Hemisphere == "Right")
  results$whole_thickness_4group_right_sensitivity <- run_cortex(sensitivity_dat, "_thickavg", "diagnosis4") %>% filter(Hemisphere == "Right")
  results$whole_thickness_3group_right_sensitivity <- run_cortex(sensitivity_dat, "_thickavg", "diagnosis3") %>% filter(Hemisphere == "Right")
  results$surface_area_4group_right_sensitivity <- run_cortex(sensitivity_dat, "_surfavg", "diagnosis4", area = TRUE) %>% filter(Hemisphere == "Right")
  results$surface_area_3group_right_sensitivity <- run_cortex(sensitivity_dat, "_surfavg", "diagnosis3", area = TRUE) %>% filter(Hemisphere == "Right")
} else {
  message("Sensitivity analyses skipped: private exclude_right_euler_sensitivity flag missing")
}

saveRDS(results, file.path(private_dir, "cortical_aggregate_results_private.rds"))
message("Cortical aggregate results saved to private output directory")
