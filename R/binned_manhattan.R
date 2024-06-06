#' Preprocess GWAS Result for Binned Manhattan Plot
#' 
#' @export
binned_manhattan_plot <- function(x, ...) UseMethod("binned_manhattan_plot")

#' Preprocess GWAS Result for Binned Manhattan Plot
#' 
#' @export
binned_manhattan_plot.MPdataBinned <- function(
  x, ...
) {
  ggplot(x$data) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = log10(n))) +
    scale_y_continuous(
      expand = c(0.02, 0.01)
    ) +
    scale_x_continuous(
      name = "Chromosome",
      breaks = x$center_pos,
      labels = x$chr.labels,
      expand = c(0.01, 0.01),
      limits = c(min(x$start_pos), max(x$end_pos))
    ) +
    geom_hline(
      yintercept = -log10(x$signif),
      linetype = "dashed",
      color = x$signif.col
    ) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    ylab(expression(-log[10](p))) +
    ggtitle(label = waiver())
}
