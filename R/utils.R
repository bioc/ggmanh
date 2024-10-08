# check preprocessing arguments
preprocess_arg_check <- function(
    x, 
    chromosome = NULL, 
    signif = NULL, 
    signif.col = NULL,
    pval.colname = NULL,
    chr.colname = NULL,
    pos.colname = NULL, 
    preserve.position = NULL, 
    pval.log.transform = NULL,
    chr.gap.scaling = NULL
) {
  
  preprocess_checklist <- list()
  
  # signif input check
  if (!is.null(signif)) {
    preprocess_checklist$signif.col <- signif.col
    # check significance cutoff exists
    if (length(signif) < 1) stop("At least one significance threshold should be provided.")
    
    # check that significance cutoff is numeric
    if (!is.numeric(signif)) stop("signif should be a numeric vector.")
    
    # check signif.col
    if (is.null(signif.col)) {
      preprocess_checklist$signif.col <- rep("grey", length(signif))
      preprocess_checklist$signif.col[which.max(-log10(signif))] <- "black"
    } else if (!all(valid_colors(signif.col))) {
      warning("invalid signif.col colors... using default colors.")
      preprocess_checklist$signif.col <- rep("grey", length(signif))
      preprocess_checklist$signif.col[which.max(-log10(signif))] <- "black"
    } else if (length(signif) != length(signif.col)) {
      warning("Length of signif and signif.col do not match.")
      if (length(signif.col) > length(signif)) {
        signif.col <- signif.col[1:length(signif)]
      } else {
        signif.col <- rep(signif.col, length.out = length(signif))
      }
      preprocess_checklist$signif.col <- signif.col
    }
  }
  
  # check that the supplied column names exist
  if (any(c(!is.null(pval.colname), !is.null(chr.colname), !is.null(pos.colname)))) {
    if (!all(c(is.character(pval.colname), is.character(chr.colname), is.character(pos.colname)))) {
      stop("Column names should be characters.")
    }
    if (!all(c(pval.colname, chr.colname, pos.colname) %in% colnames(x))) {
      tmp <- c(pval.colname, chr.colname, pos.colname)
      stop(sprintf("Column name(s) not in data: %s.", paste0(tmp[!(tmp %in% colnames(x))], collapse = ", ")))
    }
    if (pval.colname == "log10pval") {
      stop("Choose a different name for pvalue column name.")
    }
    
    if (!is.numeric(x[[pval.colname]])) stop(pval.colname, " should be a numeric column.")
    if (!is.numeric(x[[pos.colname]])) stop(pos.colname, " should be a numeric column.")
  }
  
  # check that the supplied chromosomes exist
  if (!is.null(chromosome)) {
    valid_chr(x, chromosome, chr.colname)
  }
  
  # check that the values in p-value column are valid
  if (!is.null(pval.log.transform)) {
    if (pval.log.transform) {
      if (any(x[[pval.colname]] < 0, na.rm = TRUE) | any(x[[pval.colname]] > 1, na.rm = TRUE)) {
        stop("p.value is a probability between 0 and 1.")
      }
    }
  }
  
  if (!is.null(preserve.position)) {
    if (length(preserve.position) != 1 | !is.logical(preserve.position)) {
      stop("preserve.position should be TRUE or FALSE.")
    }
  }
  
  if (!is.null(chr.gap.scaling)) {
    if (!is.numeric(chr.gap.scaling)) {
      warning("chr.gap.scaling should be numeric. Setting value to 1.")
      chr.gap.scaling <- 1
    } else if (length(chr.gap.scaling) > 1) {
      warning("Multiple values for chr.gap.scaling provided. Only using the first value.")
      chr.gap.scaling <- chr.gap.scaling[1]
    }
    if (chr.gap.scaling < 0) {
      warning("chr.gap.scaling should be a positive number. Setting value to 1.")
      chr.gap.scaling <- 1
    }
    preprocess_checklist$chr.gap.scaling <- chr.gap.scaling
  }
  
  return(preprocess_checklist)
}

# check valid chromosome argument
valid_chr <- function(x, chromosome, chr.colname) {
  if (length(chromosome) != 1) stop("Only 1 chromosome should be specified")
  if (!(chromosome %in% x[[chr.colname]])) stop("The supplied chromosome does not exist.")
  invisible()
}

# remove entries where position, chromosome, or pvalue is missing
remove_na <- function(x, chr.colname = NULL, pos.colname = NULL, pval.colname = NULL) {
  n_before <- nrow(x)
  x <- tidyr::drop_na(x, tidyr::all_of(c(chr.colname, pos.colname, pval.colname)))
  n_after <- nrow(x)
  
  present_col <- paste(
    c("chromosome", "position", "p-value")[
      c(!is.null(chr.colname), !is.null(pos.colname), !is.null(pval.colname))
      ], collapse = "/")
  
  if (n_before != n_after) {
    warning("Removed ", n_before - n_after, " rows due to missing ", present_col, ".\n")
  }
  if (n_after < 1) {
    stop("Empty rows after omitting missing ", present_col, ".\n")
  }
  return(x)
}

set_thin_logical <- function(thin, chromosome) {
  if (is.null(thin)) {
    if (is.null(chromosome)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  return(thin)
}

# TEMPORARY: replace entries where p-value is zero with the minimum
# used by manhattan data preprocess andqqunif
replace_0_pval <- function(x) {
  zero_pval <- which(x == 0)
  if (length(zero_pval) > 0) {
    warning("Replacing p-value of 0 with the minimum.")
    x[zero_pval] <- min(x[-zero_pval], na.rm = TRUE)
  }
  return(x)
}

# remove entries where p-value is zero
remove_0_pval <- function(x) {
  zero_pval <- which(x == 0)
  if (length(zero_pval) > 0) {
    warning("Removing p-value of 0.")
    x <- x[-zero_pval]
  }
  return(x)
}

set_chr_col <- function(chr.col, nchr, chr.order) {
  if (is.null(chr.col)) {
    chr.col <- stats::setNames(rep_len(RColorBrewer::brewer.pal(8, "Dark2"), nchr), chr.order)
  } else {
    if (!all(valid_colors(chr.col))) {
      warning("Invalid chr.col colors. Using default colors")
      chr.col <- stats::setNames(rep_len(RColorBrewer::brewer.pal(8, "Dark2"), nchr), chr.order)
    }
    if (!is.null(names(chr.col))) {
      if (!all(chr.order %in% names(chr.col))) {
        stop("names(chr.col) is missing values from chr.colname.")
      }
    } else {
      if (nchr > length(chr.col)) {
        warning("chr.col is recycled to match chr.order length")
        chr.col <- rep(chr.col, length.out = nchr)
        names(chr.col) <- chr.order
      } else {
        warning(paste0("Using first ", nchr, " colors for chr.col."))
        chr.col <- chr.col[1:nchr]
        names(chr.col) <- chr.order
      }
    }
  }
  return(chr.col)
}

set_highlight_col <- function(x, highlight.colname, highlight.col) {
  if (!(highlight.colname %in% colnames(x))) stop(paste0(highlight.colname, " not in data."))
  highlight.levels <- unique(x[[highlight.colname]])
  if (is.null(highlight.col)) {
    highlight.col <- RColorBrewer::brewer.pal(length(highlight.levels), "Dark2")
    names(highlight.col) <- as.character(highlight.levels)
  } else {
    if (!all(valid_colors(highlight.col))) stop("Please provide valid colors.")
    if (!is.null(names(highlight.col))) {
      if (!all(highlight.levels %in% names(highlight.col))) {
        stop("names(highlight.col) is missing values from column ", highlight.colname, ".")
      }
    } else {
      if (length(highlight.levels) > length(highlight.col)) {
        warning("highlight.col is recycled to match unique values of ", highlight.colname, ".")
        highlight.col <- rep(highlight.col, length.out = length(highlight.levels))
        names(highlight.col) <- highlight.levels
      } else if (length(highlight.levels) < length(highlight.col)) {
        warning("Using first ", length(highlight.levels), " colors.")
        highlight.col <- highlight.col[1:length(highlight.levels)]
        names(highlight.col) <- highlight.levels
      } else {
        names(highlight.col) <- highlight.levels
      }
    }
  }
  
  return(highlight.col)
}

# check that a character string is a valid color
valid_colors_ <- function(clr) {
  tryCatch(is.matrix(grDevices::col2rgb(clr)), error = function(clr) FALSE)
}

# vectorized version of valid_colors_
valid_colors <- function(clr) {
  vapply(clr, valid_colors_, logical(1))
}

# create spaced points of length(pos)
sequence_along_chr_scaled <- function(pos) {
  if (any(pos < 0)) {
    pos <- pos - min(pos)
  }
  if (max(pos) != 0) {
    return(pos / max(pos))
  } else {
    return(pos)
  }
}

# create an equally spaced points of length(pos)
sequence_along_chr_unscaled <- function(pos) {
  if (length(pos) == 1) {
    return(1/2)
  } else if (length(pos) > 1) {
    return(seq(from = 0, to = 1, length.out = length(pos)))
  } else {
    stop("Invalid pos")
  }
}

# concatenate elements across the list
concat_list <- function(dflist, concat_char = "/") {
  
  if (!all(unlist(lapply(dflist, is.vector)) | unlist(lapply(dflist, is.factor)))) {
    stop("All elements in the list should be a vector.")
  }
  
  check_lengths <- lapply(dflist, length)
  if (length(unique(unlist(check_lengths))) != 1) {
    stop("Length of all list elements should be equal.")
  }
  
  if (!is.character(concat_char)) stop("concat_char should be of character type.")
  
  if (length(concat_char) > 1) {
    warning("concat_char should be a character vector of length 1. Using first element.")
    concat_char <- concat_char[1]
  } else if (length(concat_char) < 1) {
    warning("concat_char should be a character of length 1. Using \"/\" by default.")
  }
  
  if (nchar(concat_char) < 1) {
    warning("concat_char should be a character of length 1. Using \"/\" by default.")
  }
  
  dflist <- lapply(dflist, function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    return(x)
  })
  dflist <- unname(dflist)
  dflist$sep <- concat_char
  dflist <- do.call(paste, dflist)
  dflist <- gsub(paste0("(", concat_char, ")", "+$"), "", dflist)
  dflist <- gsub(paste0("^", "(", concat_char, ")", "+"), "", dflist)
  return(dflist)
}

# concatenate columns of data.frame and produce a character vector
concat_df_cols <- function(df, concat_char = "/") {
  if (!is.data.frame(df)) stop("df should be a data.frame.")
  if (nrow(df) == 0) {
    return("")
  }
  if (length(concat_char) > 1) {
    warning("concat_char should be a character vector of length 1. Using first element.")
    concat_char <- concat_char[1]
  } else if (length(concat_char) < 1) {
    warning("concat_char should be a character of length 1. Using \"/\" by default.")
  }
  if (nchar(concat_char) < 1) {
    warning("concat_char should be a character of length 1. Using \"/\" by default.")
  }
  
  return(concat_list(as.list(df), concat_char))
}

# check that gds node exists
gds_node_exists <- function(gds, nodes) {
  all(nodes %in% gdsfmt::ls.gdsn(gds, recursive = TRUE))
}

# sample from a vector iff the number to sample from is below the length of x
sample_vec <- function(x, n) {
  if (length(x) == 1) {
    return(x)
  } else if (length(x) <= n) {
    return(x)
  } else {
    return(sample(x, size = n, replace = FALSE))
  }
}

data_is_unsorted <- function(x, chr.colname, pos.colname) {
  chr_unsorted <- is.unsorted(x[[chr.colname]])
  pos_unsorted <- tapply(x[[pos.colname]], x[[chr.colname]], is.unsorted, default = FALSE)
  return(chr_unsorted | any(pos_unsorted))
}

get_chr_pos_info <- function(chr_width, chr_gap_scaling = 1) {
  nchr <- length(chr_width)
 
  # gap between chromosome (should be robust with different lengths of chromosme)
  chr_gap <- 0.15 / 24 * nchr * chr_gap_scaling
  
  # starting x-coordinate for each chr
  start_pos <- c(0, cumsum(chr_width)[-nchr]) + ((1:nchr - 1) * chr_gap)
  names(start_pos) <- names(chr_width)
  
  # ending x-coordinate for each chr
  end_pos <- start_pos + chr_width
  
  # middle x-coordinate for each chr... used for x axis labelling
  center_pos <- (start_pos + end_pos) / 2
  
  pos_info_list <- list(
    chr_gap = chr_gap,
    chr_width = chr_width,
    start_pos = start_pos,
    center_pos = center_pos,
    end_pos = end_pos
  )
  
  return(pos_info_list)
}

calc_new_pos_ <- function(new_pos_unscaled, chr, chr_pos_info, rescale = TRUE, reposition = TRUE) {
  if (rescale) {
    new_pos_unscaled <- new_pos_unscaled * 
      unname(chr_pos_info$chr_width[as.character(chr)])
  }
  if (reposition) {
    new_pos_unscaled <- new_pos_unscaled +
      unname(chr_pos_info$start_pos[as.character(chr)])
  }
  return(new_pos_unscaled)
}

#' Calculate new x-position of each point
#'
#' Calculate the actual x-positions of each point used for
#' the manhattan plot. \code{MPdata} object contains the unscaled positions
#' that has not been positioned according to the relative position and width
#' of each chromosome.
#'
#' @param mpdata an \code{MPdata} object. 
#'
#' @details
#' This is used calculate the actual positions used for the
#' inside \code{manhattan_plot} function. It was designed this way should the 
#' scaling and relative positioning of each chromosome be changed (e.g. gap
#' between the )
#'
#' @return a \code{numeric} vector containing the scaled x-positions.
#'
#' @examples
#'
#' gwasdat <- data.frame(
#'   "chromosome" = rep(1:5, each = 30),
#'   "position" = c(replicate(5, sample(1:300, 30))),
#'   "pvalue" = rbeta(150, 1, 1)^5
#' )
#'
#' mpdata <- manhattan_data_preprocess(
#'   gwasdat, pval.colname = "pvalue", chr.colname = "chromosome", pos.colname = "position",
#'   chr.order = as.character(1:5)
#' )
#'
#' calc_new_pos(mpdata)
#'
#' @rdname calc_new_pos
#' @export
calc_new_pos <- function(mpdata) {
  calc_new_pos_(mpdata$data$new_pos_unscaled, mpdata$data[[mpdata$chr.colname]], mpdata$chr.pos.info)
}

get_background_panel_df <- function(chr_pos_info, bg_colors = c("grey90", "white")) {
  chr_gap_center <- chr_pos_info$end_pos + chr_pos_info$chr_gap/2
  chr_gap_center <- chr_gap_center[-length(chr_gap_center)]
  panel_pos_df <- data.frame(
    xmin = c(-Inf, chr_gap_center),
    xmax = c(chr_gap_center, Inf),
    ymin = -Inf,
    ymax = Inf,
    panel_col = rep(bg_colors, length.out = length(chr_pos_info$chr_width))
  )
  
  return(panel_pos_df)
}


# utility function for converting a list of formulas to a list of expressions
# this is used for allowing the user to calculate custom summary function
# for each bin
# convert list of formula  to list of expressions
convert_formula_list <- function(lst) {
  expr_list <- alist()
  for (i in 1:length(lst)) {
    expr_list[[rlang::f_lhs(lst[[i]])]] <- rlang::f_rhs(lst[[i]])
  }
  return(expr_list)
}

# use summarise function with a list of formulas
summarise_with_list <- function(dat, formula_list) {
  do.call(
    function(...) dplyr::summarise(dat, ..., .groups = "drop"),
    convert_formula_list(formula_list)
  )
}

check_formula_list <- function(dat_cnames, lst) {
  for (i in 1:length(lst)) {
    if (!rlang::is_formula(lst[[i]])) {
      stop("All elements in the expression list should be a formula.")
    } else if (length(lst[[i]]) != 3) {
      stop("All formulas in the expression list should be two sided.")
    } else if (!all(all.vars(lst[[i]][-2]) %in% dat_cnames)) {
      stop("some variables in the RHS do not exist.\n", '  formula: ', deparse1(lst[[i]]))
    }
  }
}

