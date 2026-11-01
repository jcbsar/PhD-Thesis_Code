# Final 2026 primary subcortical models — reconstructed from August source.
# Requires private, locally prepared final dataset. No participant-level output is written.
suppressPackageStartupMessages({
  library(dplyr)
  library(emmeans)
})
private_dir <- Sys.getenv("ABHGP_PRIVATE_OUTPUT_DIR")
if (!nzchar(private_dir) || !file.exists(file.path(private_dir, "final_prepared_private.rds"))) {
  stop("Run the preparation script first; set ABHGP_PRIVATE_OUTPUT_DIR")
}
dat <- readRDS(file.path(private_dir, "final_prepared_private.rds"))

primary_rois <- c("Thalamus", "Caudate", "Putamen", "Pallidum", "Hippocampus",
                  "Amygdala", "Accumbens", "LateralVentricles", "Brainstem")
stopifnot(length(primary_rois) == 9L, all(primary_rois %in% names(dat)))

contrast_spec <- list(
  diagnosis4 = list(
    "GGE vs Control" = c(-1, 1, 0, 0),
    "FE vs Control" = c(-1, 0, 1, 0),
    "FDS vs Control" = c(-1, 0, 0, 1)
  ),
  diagnosis3 = list(
    "Epilepsy vs Control" = c(-1, 1, 0),
    "FDS vs Control" = c(-1, 0, 1),
    "Epilepsy vs FDS" = c(0, 1, -1)
  )
)

fit_primary <- function(data, group_var, log_ventricles = FALSE) {
  stopifnot(group_var %in% c("diagnosis4", "diagnosis3"))
  rois <- primary_rois
  out <- lapply(rois, function(roi) {
    response <- if (log_ventricles && roi == "LateralVentricles") {
      if (any(data[[roi]] <= 0, na.rm = TRUE)) stop("Log ventricular volume requires positive values")
      "log(LateralVentricles)"
    } else roi
    fit <- lm(stats::as.formula(paste(response, "~", group_var, "+ age + sex + ICV")), data = data)
    em <- emmeans(fit, specs = group_var)
    ans <- as.data.frame(summary(contrast(em, method = contrast_spec[[group_var]]),
                                 infer = c(TRUE, TRUE)))
    ans$ROI <- roi
    ans$scale <- if (log_ventricles && roi == "LateralVentricles") "log" else "original"
    ans
  })
  bind_rows(out) %>%
    group_by(contrast) %>%
    mutate(q_BH = p.adjust(p.value, method = "BH")) %>%
    ungroup() %>%
    select(ROI, contrast, scale, estimate, SE, df, lower.CL, upper.CL, t.ratio, p.value, q_BH)
}

primary_4group <- fit_primary(dat, "diagnosis4")
primary_3group <- fit_primary(dat, "diagnosis3")

# Nine-ROI sensitivity family: same models except log transformation of
# lateral ventricular volume. This follows the appendix-described FDR family.
# NUMERICAL EQUIVALENCE IS NOT YET VERIFIED; do not treat derived q as final.
primary_3group_log_vent_sensitivity <- fit_primary(dat, "diagnosis3", log_ventricles = TRUE)

# Output only aggregate, region-wise contrast results to private directory.
saveRDS(list(primary_4group = primary_4group,
             primary_3group = primary_3group,
             primary_3group_log_vent_sensitivity = primary_3group_log_vent_sensitivity),
        file.path(private_dir, "primary_aggregate_results_private.rds"))
message("Primary aggregate results saved to private output directory")
