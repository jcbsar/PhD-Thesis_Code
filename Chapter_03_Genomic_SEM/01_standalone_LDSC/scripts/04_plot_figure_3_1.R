
#!/usr/bin/env Rscript

# ============================================================
# Figure 3.1: SNP-based heritability and genetic correlations
# ============================================================
#
# Generates Figure 3.1A and Figure 3.1B using summary-level
# results from standalone linkage disequilibrium score
# regression (LDSC).
#
# Figure 3.1A: SNP-based heritability estimates with 95% CIs.
# Figure 3.1B: Pairwise genetic correlation heatmap.
#
# This script generates figures from saved LDSC results.
# It does not perform LDSC or process GWAS summary statistics.
#
# Requirements: R, ggplot2, pheatmap
#
# Run:
# Rscript scripts/04_plot_figure_3_1.R
# ============================================================


# ------------------------------------------------------------
# 1. Directories and dependencies
# ------------------------------------------------------------

command_file <- grep(
  "^--file=",
  commandArgs(trailingOnly = FALSE),
  value = TRUE
)

if (length(command_file) != 1L) {
  stop(
    "Run using: Rscript scripts/04_plot_figure_3_1.R"
  )
}

script_file <- sub("^--file=", "", command_file)

root <- normalizePath(
  file.path(dirname(script_file), ".."),
  mustWork = TRUE
)

data_dir <- file.path(root, "data")
out_dir <- file.path(root, "outputs", "figures")

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

required_packages <- c("ggplot2", "pheatmap")

for (package in required_packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      "Missing R package: ", package,
      ". Install it before generating figures."
    )
  }
}


# ------------------------------------------------------------
# 2. Figure colours
# ------------------------------------------------------------

blue <- "#2166AC"
orange <- "#E66101"


# ============================================================
# FIGURE 3.1A: SNP-BASED HERITABILITY
# ============================================================

# ------------------------------------------------------------
# 3. Import and validate heritability estimates
# ------------------------------------------------------------

h2 <- read.csv(
  file.path(data_dir, "figure_3_1A_heritability.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns <- c(
  "Phenotype",
  "h2",
  "SE",
  "Group",
  "CI_lower",
  "CI_upper"
)

stopifnot(
  nrow(h2) == 16L,
  all(required_columns %in% names(h2)),
  !anyDuplicated(h2$Phenotype),
  all(is.finite(h2$h2)),
  all(is.finite(h2$SE)),
  all(is.finite(h2$CI_lower)),
  all(is.finite(h2$CI_upper)),
  all(h2$SE >= 0)
)

# Verify confidence intervals using estimate +/- 1.96 SE.
stopifnot(
  all(
    abs(
      h2$CI_lower - (h2$h2 - 1.96 * h2$SE)
    ) < 1e-6
  ),
  all(
    abs(
      h2$CI_upper - (h2$h2 + 1.96 * h2$SE)
    ) < 1e-6
  )
)

# Use the publication label GTCSA.
h2$Phenotype[h2$Phenotype == "GTCS"] <- "GTCSA"

# Specify the phenotype order.
plot_order <- c(
  "GGE",
  "FE",
  "GTCSA",
  "CAE",
  "JME",
  "JAE",
  "F-HS",
  "F-Neg",
  "F-Oth",
  "ACC",
  "AMY",
  "BRS",
  "CAU",
  "PAL",
  "PUT",
  "THA"
)

stopifnot(
  setequal(h2$Phenotype, plot_order),
  setequal(
    unique(h2$Group),
    c("Epilepsy", "Subcortical brain volumes")
  )
)

h2$Phenotype <- factor(
  h2$Phenotype,
  levels = rev(plot_order)
)


# ------------------------------------------------------------
# 4. Generate heritability plot
# ------------------------------------------------------------

pA <- ggplot2::ggplot(
  h2,
  ggplot2::aes(
    x = h2,
    y = Phenotype,
    colour = Group
  )
) +

  # 95% confidence intervals
  ggplot2::geom_segment(
    ggplot2::aes(
      x = CI_lower,
      xend = CI_upper,
      y = Phenotype,
      yend = Phenotype
    ),
    linewidth = 0.55
  ) +

  # Point estimates
  ggplot2::geom_point(
    ggplot2::aes(shape = Group),
    size = 2.9
  ) +

  # Reference line at zero
  ggplot2::geom_vline(
    xintercept = 0,
    colour = "grey65",
    linetype = "dashed"
  ) +

  # Phenotype group colours
  ggplot2::scale_colour_manual(
    values = c(
      "Epilepsy" = blue,
      "Subcortical brain volumes" = orange
    )
  ) +

  # Phenotype group symbols
  ggplot2::scale_shape_manual(
    values = c(
      "Epilepsy" = 16,
      "Subcortical brain volumes" = 17
    )
  ) +

  ggplot2::labs(
    x = expression("SNP-based heritability (" * h^2 * ")"),
    y = "Phenotype",
    colour = NULL,
    shape = NULL
  ) +

  ggplot2::theme_classic(
    base_size = 11
  ) +

  ggplot2::theme(
    legend.position = "bottom",
    axis.title.y = ggplot2::element_text(size = 11),
    axis.text.y = ggplot2::element_text(size = 10)
  )


# ============================================================
# FIGURE 3.1B: GENETIC CORRELATIONS
# ============================================================

# ------------------------------------------------------------
# 5. Import genetic correlation matrix
# ------------------------------------------------------------

heatmap_order <- c(
  "ACC",
  "AMY",
  "BRS",
  "CAU",
  "PAL",
  "PUT",
  "THA",
  "CAE",
  "GTCS",
  "JAE",
  "JME",
  "F-HS",
  "F-Neg",
  "F-Oth"
)

M <- as.matrix(
  read.csv(
    file.path(data_dir, "figure_3_1B_rg_matrix.csv"),
    row.names = 1,
    check.names = FALSE
  )
)

storage.mode(M) <- "double"

# Confirm matrix dimensions, phenotypes and symmetry.
stopifnot(
  identical(dim(M), c(14L, 14L)),
  setequal(rownames(M), heatmap_order),
  setequal(colnames(M), heatmap_order),
  all(is.finite(M)),
  max(abs(M - t(M))) < 1e-9,
  all(abs(diag(M) - 1) < 1e-9)
)

# Arrange phenotypes in the specified order.
M <- M[
  heatmap_order,
  heatmap_order,
  drop = FALSE
]

# Apply publication phenotype labels.
plot_labels <- heatmap_order

plot_labels[
  plot_labels == "GTCS"
] <- "GTCSA"

dimnames(M) <- list(
  plot_labels,
  plot_labels
)


# ------------------------------------------------------------
# 6. Prepare lower-triangular correlation matrix
# ------------------------------------------------------------

# Limit the colour scale to [-1, 1].
# Original correlation estimates remain unchanged in M.
M_colour <- pmax(
  -1,
  pmin(1, M)
)

# Display the lower triangle and diagonal only.
M_colour[
  upper.tri(M_colour)
] <- NA_real_

# Display original correlation estimates as text.
number_labels <- matrix(
  "",
  nrow = 14,
  ncol = 14,
  dimnames = dimnames(M)
)

number_labels[
  lower.tri(M, diag = TRUE)
] <- sprintf(
  "%.2f",
  M[lower.tri(M, diag = TRUE)]
)

diag(number_labels) <- "1"


# ------------------------------------------------------------
# 7. Generate genetic correlation heatmap
# ------------------------------------------------------------

heat_colors <- grDevices::colorRampPalette(
  c(blue, "white", orange)
)(201)

heat_breaks <- seq(
  -1,
  1,
  length.out = length(heat_colors) + 1L
)

heat <- pheatmap::pheatmap(
  M_colour,

  cluster_rows = FALSE,
  cluster_cols = FALSE,

  color = heat_colors,
  breaks = heat_breaks,

  display_numbers = number_labels,
  number_color = "#313131",

  fontsize_number = 8,
  fontsize_row = 9,
  fontsize_col = 9,

  angle_col = "0",

  na_col = "white",
  border_color = NA,

  legend_breaks = c(
    -1,
    -0.5,
    0,
    0.5,
    1
  ),

  silent = TRUE
)


# ============================================================
# EXPORT FIGURES
# ============================================================

# ------------------------------------------------------------
# 8. Save Figure 3.1A
# ------------------------------------------------------------

ggplot2::ggsave(
  file.path(
    out_dir,
    "Figure_3_1A_heritability.png"
  ),
  plot = pA,
  width = 7.5,
  height = 6.5,
  dpi = 300,
  bg = "white"
)

ggplot2::ggsave(
  file.path(
    out_dir,
    "Figure_3_1A_heritability.pdf"
  ),
  plot = pA,
  width = 7.5,
  height = 6.5,
  bg = "white"
)


# ------------------------------------------------------------
# 9. Save Figure 3.1B
# ------------------------------------------------------------

save_heat <- function(file, type) {

  if (type == "png") {
    grDevices::png(
      file,
      width = 1700,
      height = 1500,
      res = 200,
      bg = "white"
    )
  } else {
    grDevices::pdf(
      file,
      width = 8.5,
      height = 7.5,
      bg = "white"
    )
  }

  grid::grid.newpage()
  grid::grid.draw(heat$gtable)

  grDevices::dev.off()
}

save_heat(
  file.path(
    out_dir,
    "Figure_3_1B_genetic_correlations.png"
  ),
  "png"
)

save_heat(
  file.path(
    out_dir,
    "Figure_3_1B_genetic_correlations.pdf"
  ),
  "pdf"
)


# ------------------------------------------------------------
# 10. Save combined Figure 3.1
# ------------------------------------------------------------

save_combined <- function(file, type) {

  if (type == "png") {
    grDevices::png(
      file,
      width = 3200,
      height = 1400,
      res = 200,
      bg = "white"
    )
  } else {
    grDevices::pdf(
      file,
      width = 16,
      height = 7,
      bg = "white"
    )
  }

  grid::grid.newpage()

  grid::pushViewport(
    grid::viewport(
      layout = grid::grid.layout(
        1,
        2,
        widths = grid::unit(
          c(1.05, 1.0),
          "null"
        )
      )
    )
  )

  # Panel A: Heritability
  grid::pushViewport(
    grid::viewport(
      layout.pos.row = 1,
      layout.pos.col = 1
    )
  )

  grid::grid.draw(
    ggplot2::ggplotGrob(pA)
  )

  grid::grid.text(
    "A",
    x = grid::unit(0.015, "npc"),
    y = grid::unit(0.985, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(
      fontface = "bold",
      fontsize = 15
    )
  )

  grid::popViewport()

  # Panel B: Genetic correlations
  grid::pushViewport(
    grid::viewport(
      layout.pos.row = 1,
      layout.pos.col = 2
    )
  )

  grid::grid.draw(
    heat$gtable
  )

  grid::grid.text(
    "B",
    x = grid::unit(0.015, "npc"),
    y = grid::unit(0.985, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(
      fontface = "bold",
      fontsize = 15
    )
  )

  grid::popViewport(2)

  grDevices::dev.off()
}

save_combined(
  file.path(
    out_dir,
    "Figure_3_1_combined.png"
  ),
  "png"
)

save_combined(
  file.path(
    out_dir,
    "Figure_3_1_combined.pdf"
  ),
  "pdf"
)


# ------------------------------------------------------------
# 11. Completion
# ------------------------------------------------------------

message(
  "Figure 3.1 exports saved to: ",
  out_dir
)
