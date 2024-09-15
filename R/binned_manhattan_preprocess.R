#' Preprocess GWAS Result for Binned Manhattan Plot
#' 
#' Preprocess a result from Genome Wide Association Study before creating a
#' binned manhattan plot. Works similar to \code{\link{manhattan_data_preprocess}}.
#' Returns a \code{MPdataBinned} object. It can be created using a \code{data.frame}
#' or a \code{MPdata} object. Go to details to read how to use \code{summarise.expression.list}.
#' 
#' @param x a \code{data.frame} or any other extension of a data frame. It can also be a \code{MPdata} object.
#' @param ... Ignored
#' @inheritParams manhattan_data_preprocess
#' 
#' @details
#' 
#' If \code{x} is a data frame or something alike, then it creates a \code{MPdata} object first
#' and then builds \code{MPdataBinned} S3 object.
#' 
#' \code{x} can also be a \code{MPdata} object. Be sure to check if \code{thin} has been applied because this can
#' affect what's being aggregated such as number of variables in each bin.
#'
#' Positions of each point relative to the plot are first calculated 
#' via \code{\link{manhattan_data_preprocess}}.
#' Then the data is binned into blocks. \code{bins.x} indicates number of blocks
#' allocated to the chromsome with the widest width. The number of blocks 
#' for other chromosomes is proportional to the widest chromosome.
#' \code{bins.y} indicates the number of blocks allocated to the y-axis.
#' The number may be slightly adjusted to have the block height end
#' exactly at the significance threshold.
#'
#' Since points are aggregated into bins, users have the choice
#' to freely specify expressions to summarise the data for each bin
#' through \code{summarise.expression.list} argument. This argument takes a list of 
#' two-sided formulas, where the left side is the name of the new column and 
#' the right side is the expression to calculate the column. This expression is 
#' then passed to \code{\link[dplyr]{summarise}}.
#' For example, to calculate the mean, min, max of a column named \code{beta} in each bin,
#' \code{summarise.expression.list} arument would be 
#' \preformatted{
#' # inside binned_manhattan_preprocess function
#' summarise.expression.list = list(
#'   mean_beta ~ mean(beta),
#'   min_beta ~ min(beta),
#'   max_beta ~ max(beta)
#' )
#' }
#' 
#' @returns a \code{MPdataBinned} object. This object contains necessary components 
#' for creating a binned manhattan plot.
#' 
#' @examples
#' gwasdat <- data.frame(
#'   "chromosome" = rep(1:5, each = 1500),
#'   "position" = c(replicate(5, sample(1:15000, 30))),
#'   "pvalue" = rbeta(7500, 1, 1)^5,
#'   "beta" = rnorm(7500)
#' )
#' 
#' tmp <- binned_manhattan_preprocess(
#'   gwasdat, pval.colname = "pvalue", chr.colname = "chromosome",
#'   pos.colname = "position", chr.order = as.character(1:5),
#'   bins.x = 10, bins.y = 50,
#'   summarise.expression.list = list(
#'     mean_beta ~ mean(beta, na.rm = TRUE),
#'     max_abs_beta ~ max(abs(beta), na.rm = TRUE)
#'   )
#' )
#' 
#' print(tmp)
#' 
#' @rdname binned_manhattan_preprocess
#' @export
binned_manhattan_preprocess <- function(x, ...) UseMethod("binned_manhattan_preprocess")

#' @rdname binned_manhattan_preprocess
#' @method binned_manhattan_preprocess default
#' @export
binned_manhattan_preprocess.default <- function(x, ...) stop("Provide a valid data.frame, MPdata object, or GRanges object to preprocess.")


#' @rdname binned_manhattan_preprocess
#' @method binned_manhattan_preprocess MPdata
#' 
#' @param bins.x an integer. number of blocks to horizontally span the longest chromosome
#' @param bins.y an integer. number of blocks to vertically span the plot
#' @param chr.gap.scaling a number. scaling factor for the gap between chromosomes
#' @param summarise.expression.list a list of formulas to summarise data for each bin. Check details for more information.
#' @param show.message a logical. Show warning if \code{MPdata} directly used. Set to FALSE to suppress warning.
#' 
#' @export
#' @importFrom rlang .data
binned_manhattan_preprocess.MPdata <- function(
    x, bins.x = 10, bins.y = 100, chr.gap.scaling = 0.4, 
    summarise.expression.list = NULL, show.message = TRUE,
    ...
) {
  
  # check arguments
  if (!is.numeric(bins.x)) {
    stop("bins.x should be a number.")
  }
  if (!is.numeric(bins.y)) {
    stop("bins.y should be a number.")
  }
  if (!is.numeric(chr.gap.scaling)) {
    stop("chr.gap.scaling should be a number.")
  }
  
  if (length(bins.x) != 1) {
    stop("bins.x should be a single number.")
  }
  if (length(bins.y) != 1) {
    stop("bins.y should be a single number.")
  }
  if (length(chr.gap.scaling) != 1) {
    stop("chr.gap.scaling should be a single number.")
  }
  
  # verify summarise.expression.list integrity
  if (!is.null(summarise.expression.list)) {
    if (!is.list(summarise.expression.list)) {
      stop("summarise.expression.list should be a list.")
    }
    check_formula_list(colnames(x$data), summarise.expression.list)
  }
  
  if (bins.x > 30 || bins.y > 100) {
    message("Large number of bins may affect visibility of some variants in the plot.")
  }
  
  # show warning if MPdata directly used. can be suppressed with show.message
  if (show.message) {
    message("MPdata object is directly used. Please make sure that \"thin\" has not been applied to 
            the object for accurate calcuation.")
  }
  
  # calculate the lengths of chromosome
  chr_width <- x$chr.pos.info$chr_width
  longest_chr <- which.max(chr_width)
  
  # number of bins to horizontally span each chromosome
  n_blocks <- ceiling(chr_width / chr_width[longest_chr] * bins.x)
  
  # calculate horizontal length of a block
  h_length <- unname(chr_width[longest_chr] / n_blocks[longest_chr])
  
  # gets breaks for y and vertical bin length
  y_breaks <- seq(from = 0, to = ceiling(max(x$data[[x$pval.colname]])), length.out = bins.y)
  v_length <- ceiling(max(x$data[[x$pval.colname]])) / (bins.y - 1)
  
  # set block height to end exactly at the significance threshold
  tmp <- round(-log10(x$signif[1]) / v_length)
  if (tmp < 1) {
    tmp <- 1
  }
  v_length <- -log10(x$signif[1]) / tmp
  y_breaks <- seq(from = 0, to = ceiling(max(x$data[[x$pval.colname]])), by = v_length)
  if (tail(y_breaks, 1) < ceiling(max(x$data[[x$pval.colname]]))) {
    y_breaks <- c(y_breaks, tail(y_breaks, 1) + v_length)
  }
  
  # get chromosome position info
  new_chr_pos_info <- get_chr_pos_info(
    chr_width = n_blocks * h_length,
    chr_gap_scaling = chr.gap.scaling
  )
  
  # construct data frame with new positions
  x$data$new_pos_unscaled <- calc_new_pos_(
    x$data$new_pos_unscaled, x$data[[x$chr.colname]], 
    new_chr_pos_info, reposition = FALSE)
  
  .nth_bin_x <- .nth_bin_y <- NULL
  
  new_data <- x$data %>%
    # determine which bin each observation falls under
    dplyr::mutate(
      .nth_bin_x = cut(
        .data$new_pos_unscaled, breaks = seq(0, max(.data$new_pos_unscaled), by = h_length),
        labels = FALSE, include.lowest = TRUE
      ),
      .nth_bin_y = cut(
        .data[[x$pval.colname]], breaks = y_breaks,
        labels = FALSE, include.lowest = TRUE
      )
    ) %>%
    # group by chromosome, x bin, and y bin
    dplyr::group_by(.data[[x$chr.colname]], .nth_bin_x, .nth_bin_y)
  
  # calculate summary for each bin 
  if (is.null(summarise.expression.list)) {
    summarise.expression.list <- list()
  }
  
  summarise.expression.list <- append(
    list(.n_points ~ dplyr::n()),
    summarise.expression.list
  )
  
  new_data <- summarise_with_list(new_data, summarise.expression.list)
  
  bin_info <- list(
    n_blocks = n_blocks, # number of blocks (horizontally) for each chromosome
    h_length = h_length, # width of each block
    y_breaks = y_breaks, # breaks for y
    v_length = v_length # height of each block
  )
  
  # need to set following variables to NULL to avoid "NOTE" in R CMD check
  .xmin <- .xmax <- .ymin <- .ymax <- NULL
  
  # calculate xmin, xmax, ymin, ymax for each block for geom_rect
  new_data <- new_data %>%
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
    x, bins.x = 10, bins.y = 100, chr.gap.scaling = 0.4, signif = c(5e-8, 1e-5), pval.colname = "pval",
    chr.colname = "chr", pos.colname = "pos", chr.order = NULL,
    signif.col = NULL, preserve.position = TRUE,
    pval.log.transform = TRUE, summarise.expression.list = NULL, ...
) {
  
  # verify summarise.expression.list integrity
  if (!is.null(summarise.expression.list)) {
    if (!is.list(summarise.expression.list)) {
      stop("summarise.expression.list should be a list.")
    }
    check_formula_list(colnames(x), summarise.expression.list)
  }
  
  mpdat <- manhattan_data_preprocess(
    x = x, chromosome = NULL, signif = signif, pval.colname = pval.colname,
    chr.colname = chr.colname, pos.colname = pos.colname, highlight.colname = NULL,
    chr.order = chr.order, signif.col = signif.col, chr.col = NULL, highlight.col = NULL, 
    preserve.position = preserve.position, pval.log.transform = pval.log.transform, 
    thin = FALSE, ...
  )
  
  return(
    binned_manhattan_preprocess.MPdata(
      mpdat, bins.x = bins.x, bins.y = bins.y, chr.gap.scaling = chr.gap.scaling,
      summarise.expression.list = summarise.expression.list, show.message = FALSE
    )
  )
}

#' @rdname binned_manhattan_preprocess
#' @method binned_manhattan_preprocess GRanges
#' @export
setMethod(
  "binned_manhattan_preprocess", signature = "GRanges",
  function(
    x, bins.x = 10, bins.y = 100, chr.gap.scaling = 0.4, 
    signif = c(5e-8, 1e-5), pval.colname = "pval", chr.order = NULL,
    signif.col = NULL, preserve.position = TRUE,
    pval.log.transform = TRUE, summarise.expression.list = NULL,
    ...
  ) {
    grdat <- as.data.frame(x)
    grdat$pos <- (grdat$start + grdat$end) %/% 2
    
    chr.colname <- "seqnames"
    pos.colname <- "pos"
    
    binned_manhattan_preprocess.data.frame(
      grdat, bins.x = bins.x, bins.y = bins.y, signif = signif, pval.colname = pval.colname,
      chr.colname = chr.colname, pos.colname = pos.colname, chr.order = chr.order,
      signif.col = signif.col, preserve.position = preserve.position,
      pval.log.transform = pval.log.transform, thin = FALSE,
      chr.gap.scaling = chr.gap.scaling, summarise.expression.list = summarise.expression.list,
      ...
    )
  }
)
