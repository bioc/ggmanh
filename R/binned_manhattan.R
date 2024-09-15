#' Plot Binned Manhattan Plot
#' 
#' Contrary to the traditional manhattan plot, which plots all the points, the 
#' binned manhattan plot vertically and horizontally bins the variants into blocks.
#' This speeds up plotting and reduces clutter in the case of high number of variants.
#' The colors of the blocks can also be used to summarise the variants within each
#' block and highlight certain features.
#' 
#' @details
#' Similar to \code{manhattan_plot}, this function accepts summary statistics from GWAS and plots manhattan plot.
#' The difference is that the variants are binned. 
#' The number of blocks can be controlled by \code{bins.x} and \code{bins.y}.
#' The blocks can be colored based on a column in the data frame.
#' 
#' Palette for coloring the bins can be chosen from the package \href{https://emilhvitfeldt.github.io/paletteer/}{\code{paletteer}}.
#' Only palettes available in \code{paletteer} are supported. Furthermore, what palette you can use depends on what kind of 
#' variable you are using to fill the bins. Use discrete palette for categorical variable and continuous palette for continuous variable.
#' 
#' @param x a \code{data.frame} or any other extension of a data frame. It can also be a \code{MPdata} object.
#' @param ... Ignored
#' @inheritParams manhattan_plot
#' @inheritParams binned_manhattan_preprocess
#' 
#' @rdname binned_manhattan_plot
#' @export
binned_manhattan_plot <- function(x, ...) UseMethod("binned_manhattan_plot")

#' @rdname binned_manhattan_plot
#' @method binned_manhattan_plot default
#' @export
binned_manhattan_plot.default <- function(x, ...) stop("Provide a valid a valid data.frame, MPdataBinned object, or GRanges object to plot.")

#' Preprocess GWAS Result for Binned Manhattan Plot
#' 
#' @param signif.lwd a number. Line width of the significance threshold line.
#' @param bin.outline a logical. Outline each bin. The bins are colored black.
#' @param bin.outline.alpha a number. Alpha value of the bin outline.
#' @param highlight.counts a logical. If logical, the bins are colored based on the number of points in each block.
#' @param bin.palette a character. Palette to color the bins. Only palettes supported by the package \code{paletteer} are supported. More in details.
#' @param bin.alpha a number. Alpha value of the bins.
#' @param palette.direction a number. Direction of the palette. 1 for increasing and -1 for decreasing.
#' @param nonsignif.default a character. Default color for bins that are not significant.
#' @param show.legend a logical. Show legend if bins are colored based on a variable.
#' @param legend.title a character. Title of the legend.
#' @param background.col a character. Color of the background panels. Set to \code{NULL} for no color.
#' @param background.alpha a number. Alpha value of the background panels.
#' 
#' @rdname binned_manhattan_plot
#' @method binned_manhattan_plot MPdataBinned
#' @export
binned_manhattan_plot.MPdataBinned <- function(
  # data
  x, 
  # output file name
  outfn = NULL, 
  # significance threshold line width
  signif.lwd = 1,
  # bin outline
  bin.outline = FALSE,
  bin.outline.alpha = 0.2,
  # bin fill
  highlight.colname = NULL,
  highlight.counts = TRUE,
  bin.palette = "viridis::plasma",
  bin.alpha = 0.9,
  palette.direction = 1,
  nonsignif.default = NULL,
  show.legend = TRUE,
  legend.title = NULL,
  # background
  background.col = c("grey90", "white"),
  background.alpha = 0.7,
  # plot title
  plot.title = ggplot2::waiver(),
  plot.subtitle = ggplot2::waiver(),
  plot.width = 10, plot.height = 5, 
  plot.scale = 1,
  ...
) {
  
  if (!is.numeric(signif.lwd) | length(signif.lwd) != 1) {
    warning("Invalid signif.lwd value. Setting to 1.")
    signif.lwd <- 1
  }
  
  # check if each block should be outlined
  if (bin.outline) {
    bin.outline <- bin.outline
  } else {
    bin.outline <- NULL
  }
  
  if (!is.numeric(bin.outline.alpha) | length(bin.outline.alpha) != 1) {
    warning("Invalid bin.outline.alpha value. Setting to 0.2.")
    bin.outline.alpha <- 0.2
  } else if (bin.outline.alpha < 0 | bin.outline.alpha > 1) {
    warning("Invalid bin.outline.alpha value. Setting to 0.2.")
    bin.outline.alpha <- 0.2
  }
  
  if (!is.numeric(palette.direction) | length(palette.direction) != 1) {
    warning("Invalid palette.direction value. Setting to 1.")
    palette.direction <- 1
  } else if (palette.direction <= 0) {
    palette.direction <- -1
  } else {
    palette.direction <- 1
  }
  
  if (!is.null(nonsignif.default)) {
    if (length(nonsignif.default) != 1) {
      warning("nonsignif.default should be a single value. Ignoring this argument")
      nonsignif.default <- NULL
    }
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
    
    if (any(is.na(x$data[[highlight.colname]]))) {
      warning("NA values found in highlight column. These bins may not show in the plot.")
    }
    
    # set legend title if specified
    if (is.null(legend.title) || is.na(legend.title)) {
      legend.title <- waiver()
    }
    
    # different paletteer function is used depending on if the 
    # column type is categorical or numeric
    if (is.factor(x$data[[highlight.colname]]) | is.character(x$data[[highlight.colname]])) {
      fill_scale_definition <- list(
        paletteer::scale_fill_paletteer_d(
          name = legend.title,
          palette = bin.palette,
          direction = palette.direction
        )
      )
    } else {
      fill_scale_definition <- list(
        paletteer::scale_fill_paletteer_c(
          name = legend.title,
          palette = bin.palette,
          direction = palette.direction
        )
      )
    }
    
    if (!is.null(nonsignif.default)) {
      x$data[[highlight.colname]] <- ifelse(
        x$data$.ymax <= -log10(x$signif[1]),
        nonsignif.default,
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
        palette = bin.palette,
        direction = palette.direction
      )
    )
    
  }
  
  # define following variables to avoid "NOTE" in R CMD check
  .xmin <- .xmax <- .ymin <- .ymax <- NULL
  xmin <- xmax <- ymin <- ymax <- NULL
  
  # ggplot block layer
  if (!is.null(highlight.colname) | highlight.counts) {
    manh_geom <- ggplot2::geom_rect(ggplot2::aes(
      xmin = .xmin,
      xmax = .xmax,
      ymin = .ymin,
      ymax = .ymax,
      fill = .data[[highlight.colname]],
      color = bin.outline
    ), alpha = bin.alpha)
  } else {
    manh_geom <- ggplot2::geom_rect(ggplot2::aes(
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
    
    # build background
    background_panel_geom +
    
    # binned points
    manh_geom +
    
    # fill definition
    fill_scale_definition +
    
    ggplot2::scale_y_continuous(
      expand = c(0.02, 0.01)
    ) +
    ggplot2::scale_x_continuous(
      name = "Chromosome",
      breaks = x$chr.pos.info$center_pos,
      labels = x$chr.labels,
      expand = ggplot2::expansion(0.003)
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
    ggplot2::scale_color_manual(values = scales::alpha("black", bin.outline.alpha), guide = "none")
  
  if (!is.null(outfn)) {
    ggplot2::ggsave(outfn, plot=p, width=plot.width, height=plot.height, units = "in", scale=plot.scale)
    invisible()
  } else {
    return(p)
  }
}

#' @rdname binned_manhattan_plot
#' @method binned_manhattan_plot data.frame
#' @export
binned_manhattan_plot.data.frame <- function(
  x, 
  
  # preprocess arguments
  bins.x = 10, bins.y = 100, chr.gap.scaling = 0.4, signif = c(5e-8, 1e-5), pval.colname = "pval",
  chr.colname = "chr", pos.colname = "pos", chr.order = NULL,
  signif.col = NULL, preserve.position = TRUE,
  pval.log.transform = TRUE, summarise.expression.list = NULL,
  
  # plotting arguments
  outfn = NULL, signif.lwd = 1, bin.outline = FALSE, bin.outline.alpha = 0.2,
  highlight.colname = NULL, highlight.counts = TRUE, bin.palette = "viridis::plasma",
  bin.alpha = 0.9, palette.direction = 1,
  nonsignif.default = NULL, show.legend = TRUE, legend.title = NULL,
  background.col = c("grey90", "white"), background.alpha = 0.7,
  plot.title = ggplot2::waiver(), plot.subtitle = ggplot2::waiver(),
  plot.width = 10, plot.height = 5, plot.scale = 1, ...
) {
  
  mpdat <- binned_manhattan_preprocess(
    x = x, bins.x = bins.x , bins.y = bins.y, chr.gap.scaling = chr.gap.scaling,
    signif = signif, pval.colname = pval.colname, chr.colname = chr.colname,
    pos.colname = pos.colname, chr.order = chr.order, signif.col = signif.col,
    preserve.position = preserve.position, pval.log.transform = pval.log.transform,
    summarise.expression.list = summarise.expression.list, show.message = FALSE
  )
  
  p <- binned_manhattan_plot.MPdataBinned(
    x = mpdat, outfn = outfn, signif.lwd = signif.lwd, bin.outline = bin.outline,
    bin.outline.alpha = bin.outline.alpha, highlight.colname = highlight.colname,
    highlight.counts = highlight.counts, bin.palette = bin.palette, bin.alpha = bin.alpha,
    palette.direction = palette.direction, nonsignif.default = nonsignif.default,
    show.legend = show.legend, legend.title = legend.title, background.col = background.col,
    background.alpha = background.alpha, plot.title = plot.title, plot.subtitle = plot.subtitle,
    plot.width = plot.width, plot.height = plot.height, plot.scale = plot.scale, ...
  )
  
  return(p)
}

#' @rdname binned_manhattan_plot
#' @method binned_manhattan_plot GRanges
#' @export
setMethod(
  "binned_manhattan_plot", signature = "GRanges",
  function(
    x, 
    
    # preprocess arguments
    bins.x = 10, bins.y = 100, chr.gap.scaling = 0.4, signif = c(5e-8, 1e-5), pval.colname = "pval",
    chr.order = NULL,
    signif.col = NULL, preserve.position = TRUE,
    pval.log.transform = TRUE, summarise.expression.list = NULL,
    
    # plotting arguments
    outfn = NULL, signif.lwd = 1, bin.outline = FALSE, bin.outline.alpha = 0.2,
    highlight.colname = NULL, highlight.counts = TRUE, bin.palette = "viridis::plasma",
    bin.alpha = 0.9, palette.direction = 1,
    nonsignif.default = NULL, show.legend = TRUE, legend.title = NULL,
    background.col = c("grey90", "white"), background.alpha = 0.7,
    plot.title = ggplot2::waiver(), plot.subtitle = ggplot2::waiver(),
    plot.width = 10, plot.height = 5, plot.scale = 1, ...
    
  ) {
    mpdat <- binned_manhattan_preprocess(
      x, bins.x = bins.x, bins.y = bins.y, chr.gap.scaling = chr.gap.scaling, 
      signif = signif, pval.colname = pval.colname, chr.order = chr.order,
      signif.col = signif.col, preserve.position = preserve.position,
      pval.log.transform = pval.log.transform, summarise.expression.list = summarise.expression.list,
      show.message = FALSE
    )
    
    p <- binned_manhattan_plot.MPdataBinned(
      x = mpdat, outfn = outfn, signif.lwd = signif.lwd, bin.outline = bin.outline,
      bin.outline.alpha = bin.outline.alpha, highlight.colname = highlight.colname,
      highlight.counts = highlight.counts, bin.palette = bin.palette, bin.alpha = bin.alpha,
      palette.direction = palette.direction, nonsignif.default = nonsignif.default,
      show.legend = show.legend, legend.title = legend.title, background.col = background.col,
      background.alpha = background.alpha, plot.title = plot.title, plot.subtitle = plot.subtitle,
      plot.width = plot.width, plot.height = plot.height, plot.scale = plot.scale, ...
    )
  }
)

