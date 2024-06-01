mydata <- data.frame(
  x = c(2,3,5),
  y = c(1,2,5)
)

ggplot(mydata) +
  geom_raster(aes(x = x,y = y))


head(faithfuld)

library(reshape2)
library(plotly)

df <- melt(volcano)

p <- ggplot(df, aes(Var1, Var2)) +
  geom_raster(aes(fill=value)) +
  labs(x="West to East",
       y="North to South",
       title = "Elevation map of Maunga Whau")
