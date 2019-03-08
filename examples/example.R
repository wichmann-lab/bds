library(bds)
library(ggplot2)

options(mc.cores=parallel::detectCores())

raw_data <- read.table("test.csv", sep=" ", header = TRUE)

mlds_data <- data.frame(S1=raw_data$i1, S2=raw_data$i2, S3=raw_data$i3, Resp=raw_data$Response)
scale <- bds(mlds_data, fit.lapses = FALSE)

plot.df <- data.frame(stimulus=scale$stimulus,
                      scale=scale$scale,
                      ci.low=scale$scale_summary[,"2.5%"],
                      ci.high=scale$scale_summary[,"97.5%"])

ggplot(plot.df, aes(x=stimulus,y=scale)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(x=stimulus,ymin=ci.low,ymax=ci.high))

diagnostic_plots(scale)
