library(ggmanh)

test_that("Check that good arguments are provided", {
  df <- data.frame("A" = 1:3, "B" = c('a', 'b','c'), "C" = c(0.05, 0.005, 0.0005), "D" = c("A", "A", "B"))
  expect_error(manhattan_data_preprocess(df))
  expect_error(manhattan_data_preprocess(df, pval.colname = "C", chr.colname = "B", pos.colname = "B"))
  expect_error(manhattan_data_preprocess(
    df, pval.colname = "C", chr.colname = "B", pos.colname = "A",
    highlight.colname = "D", highlight.col = c("A" = "grey")
    )
  )

  tmpdf <- df
  tmpdf[2,"B"] <- NA
  expect_warning(manhattan_data_preprocess(tmpdf, pval.colname = "C", pos.colname = "A", chr.colname = "B"))

  tmpdf <- df
  tmpdf[3,"C"] <- NA
  expect_warning(manhattan_data_preprocess(tmpdf, pval.colname = "C", pos.colname = "A", chr.colname = "B"))

  tmpdf <- df
  tmpdf[1,"A"] <- NA
  expect_warning(manhattan_data_preprocess(tmpdf, pval.colname = "C", pos.colname = "A", chr.colname = "B"))

  tmpdf <- df
  tmpdf[1,"A"] <- NA; tmpdf[2,"B"] <- NA; tmpdf[3,"C"] <- NA
  expect_warning(expect_error(manhattan_data_preprocess(tmpdf, pval.colname = "C", pos.colname = "A", chr.colname = "B")))
})

test_that("Check that preprocess works as intended", {
  df <- data.frame("pos" = c(1,3,4,5,2,4), "chr" = c('a','a','a','b','b','c'), "pval" = c(0.05,0.05,0.0005,0.005,0.000005,0.0005))
  mpdat1 <- manhattan_data_preprocess(
    x = df, pval.colname = "pval", chr.colname = "chr", pos.colname = "pos",
    chr.order = c('a','b','c'), preserve.position = FALSE, thin = FALSE
  )
  mpdat2 <- manhattan_data_preprocess(
    x = df, pval.colname = "pval", chr.colname = "chr", pos.colname = "pos",
    chr.order = c('a','b','c'), preserve.position = TRUE, thin = FALSE
  )
  lg <- 0.15 / 26 * 3 # gap between chromosomes - hard coded in manhattan_preprocess function
  expect_equal(mpdat1$data$new_pos, c(0, 1/2, 1, 1 + lg, 2 + lg, 2 + lg*2 + 1/2))
  expect_equal(mpdat2$data$new_pos, c(0.375, 1.125, 1.5, 0.4 + lg + 1.5, 1 + lg + 1.5, 0.5 + lg*2 + 2.5))
  expect_equal(mpdat1$data$pval, mpdat2$data$pval)
  expect_equal(mpdat1$data$pval, c(0.05,0.05,0.0005,0.000005,0.005,0.0005))
})

test_that("Test that the thinPoint function subsets correctly.", {
  dat <- data.frame(
    A1 = c(1:20, 20, 20),
    A2 = c(rep(1, 12), rep(1,5), rep(20, 3), 20, 20),
    B = rep(c("a", "b", "c", "d"), times = c(5, 7, 8, 2))
  )
  expect_equal(nrow(thinPoints(dat, value = "A1", n = 6, nbins = 2)), 12)
  expect_equal(nrow(thinPoints(dat, value = "A2", n = 6, nbins = 2)), 11)
  expect_equal(nrow(thinPoints(dat, value = "A1", n = 3, nbins = 2, groupBy = "B")), 13)

  tmp <- thinPoints(dat, value = "A2", n = 3, nbins = 2, groupBy = "B")
  expect_equal(nrow(tmp), 14)
  expect_true(with(tmp, (sum(B == "a") == 3) && (sum(B == "b") == 3) && (sum(B == "c") == 6) && (sum(B == "d") == 2)))
  
  dat <- data.frame(
    A = c(rep(1, 3), rep(10, 3), rep(20, 4), rep(20, 6)),
    B = c(rep("A", 10), rep("B", 6))
  )
  expect_equal(nrow(thinPoints(dat, value = "A", n = 3, nbins = 1, groupBy = "B")), 6)
  expect_equal(nrow(thinPoints(dat, value = "A", n = 3, nbins = 1)), 3)
  expect_equal(nrow(thinPoints(dat, value = "A", n = 3, nbins = 2, groupBy = "B")), 9)
})

test_that("Test that sort checking works correctly.", {
  dat <- data.frame(
    A1 = factor(c(rep("A", 3), rep("B", 4)), levels = c("A","B")),
    A2 = factor(rep(c("A","B"), length.out = 7), levels = c("A","B")),
    A3 = factor(rep("A",7), levels = c("A","B")),
    B1 = 1:7,
    B2 = c(1,2,3,6,4,5,7)
  )
  
  expect_false(data_is_unsorted(dat, "A1", "B1"))
  expect_false(data_is_unsorted(dat, "A3", "B1"))
  
  expect_true(data_is_unsorted(dat, "A2", "B1"))
  expect_true(data_is_unsorted(dat, "A2", "B2"))
  expect_true(data_is_unsorted(dat, "A1", "B2"))
  expect_true(data_is_unsorted(dat, "A3", "B2"))
})

test_that("Test that sort checking works correctly.", {
  dat <- data.frame(
    A1 = factor(c(rep("A", 3), rep("B", 4)), levels = c("A","B")),
    A2 = factor(rep(c("A","B"), length.out = 7), levels = c("A","B")),
    A3 = factor(rep("A",7), levels = c("A","B")),
    B1 = 1:7,
    B2 = c(1,2,3,6,4,5,7)
  )
  
  expect_false(data_is_unsorted(dat, "A1", "B1"))
  expect_false(data_is_unsorted(dat, "A3", "B1"))
  
  expect_true(data_is_unsorted(dat, "A2", "B1"))
  expect_true(data_is_unsorted(dat, "A2", "B2"))
  expect_true(data_is_unsorted(dat, "A1", "B2"))
  expect_true(data_is_unsorted(dat, "A3", "B2"))
})
