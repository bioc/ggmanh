#' Preprocess GWAS Result for Binned Manhattan Plot
#' 
#' @export
binned_manhattan_plot <- function(x, ...) UseMethod("binned_manhattan_plot")

#' Preprocess GWAS Result for Binned Manhattan Plot
#' 
#' @export
binned_manhattan_plot.MPdataBinned <- function(
  x, outfn = NULL, bin.palette = "viridis::plasma", signif.lwd = 1,
  bin.alpha = 0.9, bin.outline = TRUE, bin.outline.alpha = 0.2,
  highlight.counts = TRUE,
  highlight.col = NULL,
  show.legend = TRUE,
  background.col = c("grey90", "white"),
  plot.title = ggplot2::waiver(), plot.subtitle = ggplot2::waiver(),
  plot.width = 10, plot.height = 5, 
  ...
) {
  chr_pos_info <- x$chr.pos.info
  if (bin.outline) {
    bin.outline <- bin.outline
  } else {
    bin.outline <- NULL
  }
  
  if (isFALSE(show.legend)) {
    show.legend <- "none"
  } else if (show.legend %in% c("left", "right", "bottom", "top")) {
    show.legend <- show.legend
  } else {
    show.legend <- "right"
  }

  if (highlight.counts) {
    fill_color <- ".n_points"
    fill_scale_definition <- list(
      paletteer::scale_fill_paletteer_c(
        name = "# Points",
        trans = "log10",
        breaks = scales::breaks_log(n = 5),
        palette = bin.palette
      )
    )
    manh_geom <- ggplot2::geom_rect(aes(
      xmin = .xmin,
      xmax = .xmax,
      ymin = .ymin,
      ymax = .ymax,
      fill = .data[[fill_color]],
      color = bin.outline
    ), alpha = bin.alpha)
  } else if (!is.null(highlight.col)) {
    fill_color <- highlight.col
    fill_scale_definition <- list(
      paletteer::scale_fill_paletteer_c(
        palette = bin.palette
      )
    )
    manh_geom <- ggplot2::geom_rect(aes(
      xmin = .xmin,
      xmax = .xmax,
      ymin = .ymin,
      ymax = .ymax,
      fill = .data[[fill_color]],
      color = bin.outline
    ), alpha = bin.alpha)
  } else {
    manh_geom <- ggplot2::geom_rect(aes(
      xmin = .xmin,
      xmax = .xmax,
      ymin = .ymin,
      ymax = .ymax,
      color = bin.outline
    ), alpha = bin.alpha)
    fill_scale_definition <- NULL
  }
  
  # alternating background panel color
  if (!is.null(background.col)) {
    panel_pos_df <- get_background_panel_df(chr_pos_info)
    background_panel_geom <- ggplot2::geom_rect(
      data = panel_pos_df,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = panel_pos_df$panel_col,
      alpha = 0.7
    )
  } else {
    background_panel_geom <- NULL
  }
  
  p <- ggplot2::ggplot(x$data) +
    background_panel_geom +
    manh_geom +
    fill_scale_definition +
    ggplot2::scale_y_continuous(
      expand = c(0.02, 0.01)
    ) +
    ggplot2::scale_x_continuous(
      name = "Chromosome",
      breaks = chr_pos_info$center_pos,
      labels = x$chr.labels,
      expand = expansion(0.003)
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(x$signif),
      linetype = "dashed",
      color = x$signif.col,
      linewidth = signif.lwd
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = show.legend 
    ) +
    ggplot2::ylab(expression(-log[10](p))) +
    ggplot2::ggtitle(label = plot.title, subtitle = plot.subtitle) +
    scale_color_manual(values = alpha("black", bin.outline.alpha), guide = "none")
  
  if (!is.null(outfn)) {
    ggplot2::ggsave(outfn, plot=p, width=plot.width, height=plot.height, units = "in")
    invisible()
  } else {
    return(p)
  }
}
