#' Preprocess GWAS Result for Binned Manhattan Plot
#' 
#' @export
binned_manhattan_preprocess <- function(x, ...) UseMethod("binned_manhattan_preprocess")

#' @rdname binned_manhattan_preprocess
#' @method binned_manhattan_preprocess default
#' @export
binned_manhattan_preprocess.default <- function(x, ...) stop("Provide a valid data.frame, MPdata object, or GRanges object to preprocess.")

#' @rdname binned_manhattan_preprocess
#' @method binned_manhattan_preprocess MPdata
#' @export
binned_manhattan_preprocess.MPdata <- function(
  x, bins_x = 20, bins_y = 100, chr.gap.scaling = 1, 
  summarise_expression_list = NULL
) {
  
  if (!is.null(summarise_expression_list)) {
    if (!is.list(summarise_expression_list)) {
      stop("summarise_expression_list should be a list.")
    }
    check_formula_list(colnames(x$data), summarise_expression_list)
  }
  
  # calculate the lengths of chromosome
  chr_width <- x$chr.pos.info$chr_width
  longest_chr <- which.max(chr_width)
  
  # number of bins to horizontally span each chromosome
  n_blocks <- ceiling(chr_width / chr_width[longest_chr] * bins_x)
  
  # calculate horizontal length of a block
  h_length <- unname(chr_width[longest_chr] / n_blocks[longest_chr])
  
  # gets breaks for y and vertical bin length
  y_breaks <- seq(from = 0, to = ceiling(max(x$data[[x$pval.colname]])), length.out = bins_y)
  v_length <- ceiling(max(x$data[[x$pval.colname]])) / (bins_y - 1)
  
  # get chromosome position info
  new_chr_pos_info <- get_chr_pos_info(
    chr_width = n_blocks * h_length,
    chr_gap_scaling = chr.gap.scaling
    )

  # construct data frame with new positions
  x$data$new_pos_unscaled <- calc_new_pos_(
    x$data$new_pos_unscaled, x$data[[x$chr.colname]], 
    new_chr_pos_info, reposition = FALSE)
  
  new_data <- x$data |>
    # determine which bin each observation falls under
    dplyr::mutate(
      .nth_bin_x = cut(
        new_pos_unscaled, breaks = seq(0, max(new_pos_unscaled), by = h_length),
        labels = FALSE, include.lowest = TRUE
      ),
      .nth_bin_y = cut(
        .data[[x$pval.colname]], breaks = y_breaks,
        labels = FALSE, include.lowest = TRUE
      )
    ) |>
    # group by chromosome, x bin, and y bin
    dplyr::group_by(.data[[x$chr.colname]], .nth_bin_x, .nth_bin_y)
  
  # calculate summary for each bin 
  if (is.null(summarise_expression_list)) {
    summarise_expression_list <- list()
  }
  
  summarise_expression_list <- append(
    list(.n_points ~ n()),
    summarise_expression_list
  )
  
  new_data <- summarise_with_list(new_data, summarise_expression_list)
  
  bin_info <- list(
    n_blocks = n_blocks, # number of blocks (horizontally) for each chromosome
    h_length = h_length, # width of each block
    y_breaks = y_breaks, # breaks for y
    v_length = v_length # height of each block
  )
  
  # calculate xmin, xmax, ymin, ymax for each block for geom_rect
  new_data <- new_data |>
    dplyr::mutate(
      .xmin = new_chr_pos_info$start_pos[as.character(.data[[x$chr.colname]])] + 
        bin_info$h_length * (.nth_bin_x-1),
      .xmax = .xmin + bin_info$h_length,
      .ymin = bin_info$y_breaks[.nth_bin_y],
      .ymax = .ymin + bin_info$v_length
    )
  
  # update MPdata object to MPdataBinned
  mpdat_binned <- list(
    data = new_data,
    chr.pos.info = new_chr_pos_info,
    bin.info = bin_info,
    signif = x$signif,
    signif.col = x$signif.col,
    chr.labels = x$chr.labels,
    chr.colname = x$chr.colname
  )
  class(mpdat_binned) <- c("MPdataBinned")
  
  return(mpdat_binned)
}

#' @rdname binned_manhattan_preprocess
#' @method binned_manhattan_preprocess data.frame
#' @export
binned_manhattan_preprocess.data.frame <- function(
  x, bins_x = 20, bins_y = 100, chr.gap.scaling = 1, signif = c(5e-8, 1e-5), pval.colname = "pval",
  chr.colname = "chr", pos.colname = "pos", chr.order = NULL,
  signif.col = NULL, preserve.position = TRUE,
  pval.log.transform = TRUE, summarise_expression_list = NULL, ...
) {
  
  mpdat <- manhattan_data_preprocess(
    x = x, chromosome = NULL, signif = signif, pval.colname = pval.colname,
    chr.colname = chr.colname, pos.colname = pos.colname, highlight.colname = NULL,
    chr.order = chr.order, signif.col = signif.col, chr.col = NULL, highlight.col = NULL, 
    preserve.position = preserve.position, pval.log.transform = pval.log.transform, 
    thin = FALSE, ...
  )
  
  return(binned_manhattan_preprocess(
    mpdat, bins_x = bins_x, bins_y = bins_y, chr.gap.scaling = chr.gap.scaling,
    summarise_expression_list = summarise_expression_list)
  )
}
