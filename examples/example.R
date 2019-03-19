library(bds)
library(ggplot2)

options(mc.cores=parallel::detectCores())

raw_data <- read.table("test.csv", sep=" ", header = TRUE)

mlds_data <- data.frame(Resp=raw_data$Response, S1=raw_data$i1, S2=raw_data$i2, S3=raw_data$i3)
scale <- bds(mlds_data, fit.lapses = FALSE)

plot.df <- data.frame(stimulus=scale$stimulus,
                      scale=get_scale_values(scale),
                      ci.low=get_scale_credible_interval(scale)$ci.low,
                      ci.high=get_scale_credible_interval(scale)$ci.high)

ggplot(plot.df, aes(x=stimulus,y=scale)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(x=stimulus,ymin=ci.low,ymax=ci.high))

diagnostic_plots(scale)
