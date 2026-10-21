# Reproduce thesis Figure 3.1 from archived, independently verified numerical
# LDSC exports. This script does NOT rerun LD-score regression.
# The August 2026 original plotting source has been adapted/reconstructed.
# Requires R with ggplot2 and pheatmap. Run: Rscript scripts/04_plot_figure_3_1.R

command_file <- grep('^--file=', commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(command_file) != 1L) {
  stop('Run from the command line using Rscript scripts/04_plot_figure_3_1.R')
}
script_file <- sub('^--file=', '', command_file)
root <- normalizePath(file.path(dirname(script_file), '..'), mustWork = TRUE)
data_dir <- file.path(root, 'data')
out_dir <- file.path(root, 'outputs', 'figures')
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (package in c('ggplot2', 'pheatmap')) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop('Missing R package ', package, '. Install it before plotting.')
  }
}

blue <- '#2166AC'
orange <- '#E66101'

# ---- Panel A: univariate h², liability scale for epilepsy ------------------
h2 <- read.csv(file.path(data_dir, 'figure_3_1A_heritability.csv'),
               stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(nrow(h2) == 16L,
          all(c('Phenotype', 'h2', 'SE', 'Group', 'CI_lower', 'CI_upper') %in% names(h2)))
stopifnot(all(abs(h2$CI_lower - (h2$h2 - 1.96 * h2$SE)) < 1e-6))
stopifnot(all(abs(h2$CI_upper - (h2$h2 + 1.96 * h2$SE)) < 1e-6))

h2$Phenotype[h2$Phenotype == 'GTCS'] <- 'GTCSA'  # publication label
plot_order <- c('GGE','FE','GTCSA','CAE','JME','JAE',
                'F-HS','F-Neg','F-Oth','ACC','AMY','BRS','CAU','PAL','PUT','THA')
stopifnot(setequal(h2$Phenotype, plot_order))
h2$Phenotype <- factor(h2$Phenotype, levels = rev(plot_order))

pA <- ggplot2::ggplot(h2, ggplot2::aes(x = h2, y = Phenotype, colour = Group)) +
  ggplot2::geom_segment(ggplot2::aes(x = CI_lower, xend = CI_upper,
                                    y = Phenotype, yend = Phenotype), linewidth = 0.55) +
  ggplot2::geom_point(ggplot2::aes(shape = Group), size = 2.9) +
  ggplot2::geom_vline(xintercept = 0, colour = 'grey65', linetype = 'dashed') +
  ggplot2::scale_colour_manual(values = c('Epilepsy' = blue,
                                         'Subcortical brain volumes' = orange)) +
  ggplot2::scale_shape_manual(values = c('Epilepsy' = 16,
                                        'Subcortical brain volumes' = 17)) +
  ggplot2::labs(x = expression('SNP-based heritability ('*h^2*')'),
                y = 'Phenotype', colour = NULL, shape = NULL) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(legend.position = 'bottom',
                 axis.title.y = ggplot2::element_text(size = 11),
                 axis.text.y = ggplot2::element_text(size = 10))

# ---- Panel B: heatmap of 91 LDSC rg estimates ------------------------------
heatmap_order <- c('ACC','AMY','BRS','CAU','PAL','PUT','THA',
                   'CAE','GTCS','JAE','JME','F-HS','F-Neg','F-Oth')
M <- as.matrix(read.csv(file.path(data_dir, 'figure_3_1B_rg_matrix.csv'),
                        row.names = 1, check.names = FALSE))
storage.mode(M) <- 'double'
stopifnot(identical(dim(M), c(14L, 14L)),
          setequal(rownames(M), heatmap_order),
          setequal(colnames(M), heatmap_order),
          max(abs(M - t(M))) < 1e-9)
M <- M[heatmap_order, heatmap_order, drop = FALSE]

# Publish GTCSA, not the source GWAS shorthand GTCS.
plot_labels <- heatmap_order
plot_labels[plot_labels == 'GTCS'] <- 'GTCSA'
dimnames(M) <- list(plot_labels, plot_labels)

M_colour <- pmax(-1, pmin(1, M))  # cap COLOURS only; preserve actual rg estimates
M_colour[upper.tri(M_colour)] <- NA_real_
number_labels <- matrix('', nrow = 14, ncol = 14, dimnames = dimnames(M))
number_labels[lower.tri(M, diag = TRUE)] <- sprintf('%.2f', M[lower.tri(M, diag = TRUE)])
diag(number_labels) <- '1'

heat_colors <- grDevices::colorRampPalette(c(blue, 'white', orange))(201)
heat_breaks <- seq(-1, 1, length.out = length(heat_colors) + 1L)
heat <- pheatmap::pheatmap(
  M_colour,
  cluster_rows = FALSE, cluster_cols = FALSE,
  color = heat_colors, breaks = heat_breaks,
  display_numbers = number_labels,
  number_color = '#313131', fontsize_number = 8,
  fontsize_row = 9, fontsize_col = 9,
  angle_col = '0', na_col = 'white', border_color = NA,
  legend_breaks = c(-1, -0.5, 0, 0.5, 1),
  silent = TRUE
)

# ---- Save both panels plus a combined figure -------------------------------
ggplot2::ggsave(file.path(out_dir, 'Figure_3_1A_heritability.png'),
                pA, width = 7.5, height = 6.5, dpi = 300, bg = 'white')
ggplot2::ggsave(file.path(out_dir, 'Figure_3_1A_heritability.pdf'),
                pA, width = 7.5, height = 6.5, bg = 'white')

save_heat <- function(file, type) {
  if (type == 'png') {
    grDevices::png(file, width = 1700, height = 1500, res = 200, bg = 'white')
  } else {
    grDevices::pdf(file, width = 8.5, height = 7.5, bg = 'white')
  }
  grid::grid.newpage()
  grid::grid.draw(heat$gtable)
  grDevices::dev.off()
}
save_heat(file.path(out_dir, 'Figure_3_1B_genetic_correlations.png'), 'png')
save_heat(file.path(out_dir, 'Figure_3_1B_genetic_correlations.pdf'), 'pdf')

save_combined <- function(file, type) {
  if (type == 'png') {
    grDevices::png(file, width = 3200, height = 1400, res = 200, bg = 'white')
  } else {
    grDevices::pdf(file, width = 16, height = 7, bg = 'white')
  }
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(
    1, 2, widths = grid::unit(c(1.05, 1.0), 'null'))))
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid::grid.draw(ggplot2::ggplotGrob(pA))
  grid::grid.text('A', x = grid::unit(0.015, 'npc'), y = grid::unit(0.985, 'npc'),
                  just = c('left', 'top'), gp = grid::gpar(fontface = 'bold', fontsize = 15))
  grid::popViewport()
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  grid::grid.draw(heat$gtable)
  grid::grid.text('B', x = grid::unit(0.015, 'npc'), y = grid::unit(0.985, 'npc'),
                  just = c('left', 'top'), gp = grid::gpar(fontface = 'bold', fontsize = 15))
  grid::popViewport(2)
  grDevices::dev.off()
}
save_combined(file.path(out_dir, 'Figure_3_1_combined.png'), 'png')
save_combined(file.path(out_dir, 'Figure_3_1_combined.pdf'), 'pdf')
message('Figure 3.1 exports saved to: ', out_dir)
