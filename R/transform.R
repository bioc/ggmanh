get_transform_jump <- function(x) {

  # max of significance threshold, then capped by even number
  jump <- ceiling(max(x))
  jump <- if ((jump %% 2) == 0) {jump} else {jump + 1}
  return(jump)
}

get_transform <- function(
  assoc, jump, pval.colname = "log10pval", jump.rel.pos = 0.4
) {

  # create a transformation object for transforming the scales
  # "signif" method: rescale according to where you want the significant line to be
  # will be used for manipulation of scale.

  # this method only works when the minimum value is 0
  if (min(assoc[[pval.colname]]) < 0) stop("-log10(pvalue) should be positive.")

  # This is the total number of "ticks" that's used to measure the scale of the plots
  ticks <- 100

  # highest value capped by a factor of 10
  max_tick <- ceiling(max(assoc[[pval.colname]])/10)*10
  
  # levels breaks function to get the appropriate breaks
  # lower half still increments by 2
  breaks_upper <- scales::breaks_extended(10, only.loose = TRUE)(c(jump, max_tick))
  breaks_upper <- breaks_upper[breaks_upper > jump]
  
  if (max(breaks_upper) < max_tick) {
    breaks_upper <- c(breaks_upper, max(breaks_upper) + breaks_upper[2] - breaks_upper[1])
  }
  
  max_tick <- max(breaks_upper)
  breaks_manh <- c(seq(0, jump, 2), breaks_upper)
  
  # calculate the scales for bottom & top portion
  scale_1_ticks <- ticks * jump.rel.pos
  scale_unit_1 <- jump / scale_1_ticks
  scale_unit_2 <- (max_tick - jump) / (ticks - scale_1_ticks)

  scale_transform <- function(x) {
    ifelse(x <= jump, x / scale_unit_1, (x - jump) / scale_unit_2 + jump/scale_unit_1)
  }
  scale_inverse <- function(x) {
    ifelse(x <= jump/scale_unit_1, x * scale_unit_1, (x - jump/scale_unit_1) * scale_unit_2 + jump)
  }
  
  # remove y-axis tick @ the jump
  y_axis_ticks_linetype <- rep("solid", length(breaks_manh))
  y_axis_ticks_linetype[breaks_manh == jump] <- "blank"

  # create transformation object
  trans_manh <- scales::trans_new(
    name = "manhattan_scale",
    transform = scale_transform,
    inverse = scale_inverse
  )

  return(list(
    "trans" = trans_manh,
    "breaks" = breaks_manh,
    "y_axis_linetype" = y_axis_ticks_linetype,
    "jump" = jump
  ))
}
