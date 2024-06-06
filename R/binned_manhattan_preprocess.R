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
  x, bins_x = 20, bins_y = 100
) {
  # calculate the lengths of chromosome
  chr_lengths <- x$end_pos - x$start_pos
  longest_chr <- which.max(chr_lengths)
  
  # number of bins to horizontally span each chromosome
  n_blocks <- round(chr_lengths / chr_lengths[longest_chr] * bins_x)
  
  # calculate horizontal length of a block
  h_length <- chr_lengths[longest_chr] / n_blocks[longest_chr]
  
  # gets breaks for y and vertical bin length
  y_breaks <- seq(from = 0, to = ceiling(max(x$data[[x$pval.colname]])), length.out = bins_y)
  v_length <- ceiling(max(x$data[[x$pval.colname]])) / (bins_y - 1)
  
  # calculate new position info
  nchr <- length(chr_lengths)
  new_chr_width <- n_blocks * h_length
  
  lg <- 0.15 / 26 * nchr
  new_start_pos <- c(0, cumsum(new_chr_width)[-nchr]) + ((1:nchr - 1) * lg)
  names(new_start_pos) <- x$chr.labels # starting x-coordinate for each chr
  new_end_pos <- new_start_pos + new_chr_width # ending x-coordinate for each chr
  new_center_pos <- (new_start_pos + new_end_pos) / 2 # middle x-coordinate for each chr... used for x axis labelling
  
  # construct data frame with new positions
  new_data <- x$data |>
    
    # reposition and rescale all points to fit the new blocks
    group_by(.data[[x$chr.colname]]) |>
    mutate(
      new_pos = new_pos - x$start_pos[as.character(.data[[x$chr.colname]])],
      new_pos = new_pos / max(new_pos) * new_chr_width[as.character(.data[[x$chr.colname]])],
    ) |>
    ungroup() |>
    
    # determine which bin each observation falls under
    mutate(
      nth_bin_x = cut(
        new_pos, breaks = seq(0, max(new_pos), by = h_length),
        labels = FALSE, include.lowest = TRUE
      ),
      nth_bin_y = cut(
        .data[[x$pval.colname]], breaks = y_breaks,
        labels = FALSE, include.lowest = TRUE
      )
    ) |>
    
    # count the number of observations in each bin
    group_by(.data[[x$chr.colname]], nth_bin_x, nth_bin_y) |>
    summarise(
      n = n(),
      .groups = "drop"
    ) |>
    
    # calculate xmin, xmax, ymin, ymax for each block for geom_rect
    mutate(
      xmin = new_start_pos[as.character(.data[[x$chr.colname]])] + 
        h_length * (nth_bin_x-1),
      xmax = xmin + h_length,
      ymin = y_breaks[nth_bin_y],
      ymax = ymin + v_length
    ) |>
    select(-nth_bin_x, -nth_bin_y)
  
  # update MPdata object to MPdataBinned
  mpdat_binned <- list(
    data = new_data,
    chr_width = new_chr_width,
    start_pos = new_start_pos,
    center_pos = new_center_pos,
    end_pos = new_end_pos,
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
  x, bins_x = 20, bins_y = 100, signif = c(5e-8, 1e-5), pval.colname = "pval",
  chr.colname = "chr", pos.colname = "pos", chr.order = NULL,
  signif.col = NULL, preserve.position = FALSE,
  pval.log.transform = TRUE, ...
) {
  
  mpdat <- manhattan_data_preprocess(
    x = x, chromosome = NULL, signif = signif, pval.colname = pval.colname,
    chr.colname = chr.colname, pos.colname = pos.colname, highlight.colname = NULL,
    chr.order = chr.order, signif.col = signif.col, chr.col = NULL, highlight.col = NULL, 
    preserve.position = preserve.position, pval.log.transform = pval.log.transform, 
    thin = FALSE, ...
  )
  
  return(binned_manhattan_preprocess(mpdat, bins_x = bins_x, bins_y = bins_y))
}

