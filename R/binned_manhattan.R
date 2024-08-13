#' Preprocess GWAS Result for Binned Manhattan Plot
#' 
#' @export
binned_manhattan_plot <- function(x, ...) UseMethod("binned_manhattan_plot")

#' Preprocess GWAS Result for Binned Manhattan Plot
#' 
#' @export
binned_manhattan_plot.MPdataBinned <- function(
  x, outfn = NULL, signif.lwd = 1,
  bin.outline = FALSE, bin.outline.alpha = 0.2,
  highlight.colname = NULL,
  highlight.counts = TRUE,
  bin.palette = "viridis::plasma",
  bin.alpha = 0.9,
  palette.direction = 1,
  nonsignif_default = NULL,
  show.legend = TRUE,
  background.col = c("grey90", "white"),
  background.alpha = 0.7,
  plot.title = ggplot2::waiver(), plot.subtitle = ggplot2::waiver(),
  plot.width = 10, plot.height = 5, 
  ...
) {
  
  # check if each block should be outlined
  if (bin.outline) {
    bin.outline <- bin.outline
  } else {
    bin.outline <- NULL
  }
  
  # check if legend is displayed
  if (isFALSE(show.legend)) {
    show.legend <- "none"
  } else if (show.legend %in% c("left", "right", "bottom", "top")) {
    show.legend <- show.legend
  } else {
    show.legend <- "right"
  }
  
  # coloring blocks
  # use scale_fill_paletteer to allow user to choose palette
  fill_scale_definition <- NULL
  if (!is.null(highlight.colname)) {
    # if highlight.colname is specified, use the 
    # column to fill block
    
    # different paletteer function is used depending on if the 
    # column type is categorical or numeric
    if (is.factor(x$data[[highlight.colname]]) | is.character(x$data[[highlight.colname]])) {
      fill_scale_definition <- list(
        paletteer::scale_fill_paletteer_d(
          palette = bin.palette,
          direction = palette.direction
        )
      )
    } else {
      fill_scale_definition <- list(
        paletteer::scale_fill_paletteer_c(
          palette = bin.palette,
          direction = palette.direction
        )
      )
    }
    
    if (!is.null(nonsignif_default)) {
      x$data[[highlight.colname]] <- ifelse(
        x$data$.ymax <= -log10(x$signif[1]),
        nonsignif_default,
        x$data[[highlight.colname]]
      )
    }
  } else if (highlight.counts) {
    highlight.colname <- ".n_points"
    fill_scale_definition <- list(
      paletteer::scale_fill_paletteer_c(
        name = "# Points",
        trans = "log10",
        breaks = scales::breaks_log(n = 5),
        palette = bin.palette
      )
    )
    
  }
  
  # ggplot block layer
  if (!is.null(highlight.colname) | highlight.counts) {
    manh_geom <- ggplot2::geom_rect(aes(
      xmin = .xmin,
      xmax = .xmax,
      ymin = .ymin,
      ymax = .ymax,
      fill = .data[[highlight.colname]],
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
  }
  
  # alternating background panel color
  if (!is.null(background.col)) {
    panel_pos_df <- get_background_panel_df(x$chr.pos.info, bg_colors = background.col)
    background_panel_geom <- ggplot2::geom_rect(
      data = panel_pos_df,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = panel_pos_df$panel_col,
      alpha = background.alpha
    )
  } else {
    background_panel_geom <- NULL
  }
  
  # create final ggplot object
  p <- ggplot2::ggplot(x$data) +
    background_panel_geom +
    manh_geom +
    fill_scale_definition +
    ggplot2::scale_y_continuous(
      expand = c(0.02, 0.01)
    ) +
    ggplot2::scale_x_continuous(
      name = "Chromosome",
      breaks = x$chr.pos.info$center_pos,
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
