library(ggplot2)
library(tidyverse)
library(devtools)

load_all()

# create simulation data
set.seed(1000)

nsim <- 50000

simdata <- data.frame(
  "chromosome" = sample(c(1:22,"X"), size = nsim, replace = TRUE),
  "position" = sample(1:100000000, size = nsim),
  "P.value" = rbeta(nsim, shape1 = 5, shape2 = 1)^7
)

dat <- read.table("../GCST90018926_buildGRCh37.tsv", sep = "\t", header = TRUE)

# mpdat <- manhattan_data_preprocess(
#   simdata, pval.colname = "P.value", chr.colname = "chromosome", pos.colname = "position",
#   chr.order = c(1:22,"X"), thin = FALSE, preserve.position = TRUE
# )

x <- manhattan_data_preprocess(
  dat, pval.colname = "p_value", chr.colname = "chromosome", pos.colname = "base_pair_location",
  chr.order = c(1:22, "X"),
  signif.col = NULL, preserve.position = TRUE
)

x <- binned_manhattan_preprocess.data.frame(
  dat, bins_x = 15, bins_y = 70, pval.colname = "p_value",
  chr.colname = "chromosome", pos.colname = "base_pair_location", chr.order = c(1:22,"X"),
  signif.col = NULL, preserve.position = TRUE, pval.log.transform = TRUE,
  chr.gap.scaling = 0.3, summarise_expression_list = list(
    abs_max_beta ~ max(abs(beta)),
    mean_beta ~ mean(beta),
    main_direction ~ ifelse(sum(beta > 0) > sum(beta < 0), "positive", "negative")
  )
)

x2 <- binned_manhattan_preprocess.data.frame(
  dat, bins_x = 50, bins_y = 100, pval.colname = "p_value",
  chr.colname = "chromosome", pos.colname = "base_pair_location", chr.order = c(1:5,"X"),
  signif.col = NULL, preserve.position = TRUE, pval.log.transform = TRUE,
  chr.gap.scaling = 3
)

# mpdat <- manhattan_data_preprocess(
#   dat, pval.colname = "p_value", chr.colname = "chromosome", pos.colname = "base_pair_location",
#   chr.order = c(1:5,"X"), thin = FALSE, preserve.position = TRUE
# )
# 
# x2 <- binned_manhattan_preprocess.MPdata(
#   mpdat, bins_x = 40, bins_y = 100
# )

binned_manhattan_plot.MPdataBinned(
  x, highlight.counts = FALSE, highlight.colname = "abs_max_beta", 
  bin.palette = "viridis::inferno", bin.alpha = 1, bin.outline = TRUE,
  bin.outline.alpha = 0.2, nonsignif_default = 0
)
binned_manhattan_plot.MPdataBinned(x2)

manhattan_plot(
  dat, pval.colname = "p_value",
  chr.colname = "chromosome", pos.colname = "base_pair_location", chr.order = c(1:5,"X"),
  signif.col = NULL, preserve.position = TRUE, pval.log.transform = TRUE
)

ggplot(new_data) +
  geom_rect(aes(xmin = new_pos, xmax = new_pos2, ymin = log10pval, ymax = log10pval2, fill = log10(n))) +
  scale_y_continuous(
    expand = c(0.02, 0.01)
  ) +
  scale_x_continuous(
    name = "Chromosome",
    breaks = new_center_pos_tmp,
    labels = chr.order_tmp,
    expand = c(0.01, 0.01),
    limits = c(min(new_start_pos_tmp), max(new_end_pos_tmp))
  ) +
  geom_hline(
    yintercept = -log10(mpdat$signif),
    linetype = "dashed",
    color = mpdat$signif.col
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  ylab(expression(-log[10](p))) +
  ggtitle(label = waiver())

# compare to manhattan plot

mpdat2 <- manhattan_data_preprocess(
  dat, pval.colname = "p_value", chr.colname = "chromosome", pos.colname = "base_pair_location",
  chr.order = c(1:22,"X"), thin = TRUE, preserve.position = TRUE
)

manhattan_plot(mpdat2, outfn = "../test.png", rescale = FALSE)

x <- dat
head(x)
bins_x <- 20
bins_y <- 200
signif <- c(5e-8, 1e-5)
pval.colname <- "p_value"
chr.colname <- "chromosome"
pos.colname <- "base_pair_location"
highlight.colname <- FALSE
chr.order <- c(1:22,"X")
signif.col <- NULL
chr.col <- NULL
preserve.position <- FALSE
pval.log.transform <- TRUE
thin <- FALSE

### 

head(starwars)

# pass a formula to create a new column in summarise function

