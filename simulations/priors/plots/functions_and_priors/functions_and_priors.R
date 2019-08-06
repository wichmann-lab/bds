setwd("../../")
source("plot_preprocess.R")
setwd("plots/functions_and_priors")

library(latex2exp)

fn.df <- data.frame(x=numeric(), y=numeric(), fn=factor())
point.df <- data.frame(x=numeric(), y=numeric(), fn=factor())

plotlist <- list()

for (fname in names(function.zoo)) {
  pt.x <- c(0, seq(0.025,0.975, len=9), 1)
  x <- seq(0, 1, len=1000)
  y <- function.zoo[[fname]](x)
  pt.y <- function.zoo[[fname]](pt.x)
  fn.df <-data.frame(x=x, y=y, fn=fname)
  point.df <- data.frame(x=pt.x, y=pt.y, fn=fname)
  
  fn.plt <- ggplot(fn.df, aes(x=x,y=y)) + geom_line() +
    geom_point(data = point.df) +
    ggtitle(fname) +
    scale_y_continuous(breaks= c(0, 1)) +
    scale_x_continuous(breaks= c(0,1)) +
    theme(axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank()) +
    NULL
  
  plotlist[[fname]] <- fn.plt
}

save_plot("functions.pdf", plot_grid(plotlist = plotlist, nrow = 2), base_height = 4, base_width = 10)

### prior forms

#halfnormal.df <- data.frame(precision=c(-5, 0, seq(0, 60, len=1000)), density = c(0, 0, dnorm(seq(0,60, len=1000), sd=20)), lb = rep(0, times=1002))
#hn.plt <- ggplot(halfnormal.df, aes(x=precision, y=density, ymin=lb, ymax=density)) +
#  geom_line() +
#  geom_ribbon(alpha = 0.3, linetype= "blank") +
#  ggtitle(TeX("$\\sigma^{-2} \\sim half-normal(0, 20)$")) +
#  theme(axis.title = element_blank(),
#        axis.text.y = element_blank(),
#        axis.ticks.y = element_blank(),
#       axis.line.y = element_blank())

#uni.plt <- ggplot(data.frame(x=c(-0.5, 0, 0, 1, 1, 1.5), y=c(0, 0, 1, 1, 0, 0), lb = c(0, 0, 0, 0, 0, 0)), aes(x=x,y=y,ymin=lb, ymax=y)) +
#  geom_line() + geom_ribbon(alpha = 0.3, linetype= "blank") +
#  ggtitle(TeX("$\\Psi_x \\sim Uniform(0,1)$")) +
#  ylim(0, 1.5) +
#  theme(axis.title = element_blank(),
#        axis.text.y = element_blank(),
#        axis.ticks.y = element_blank(),
#        axis.line.y = element_blank())

raised.df <- data.frame(x=c(-2, seq(0, 2, len=100), seq(15, 30, len=1000), 32), y=c(0, 0.5*c(-cos(pi*seq(0,2, len=100)/2), -cos(pi*seq(15, 30, len=1000)/15)) + 0.5, 0), lb = rep(0, times=1102))
raised.plt <- ggplot(raised.df, aes(x=x,y=y,ymin=lb, ymax=y)) +
  geom_line() + geom_ribbon(alpha = 0.3, linetype= "blank") +
  ggtitle(TeX("$\\sigma^{-1} \\sim raised-cosine(0,2;15,30)$")) +
  ylim(0, 1.5) +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

beta.plt <- ggplot(data.frame(x=seq(0,1,len=1000), y = dbeta(seq(0,1, len=1000), 1, 10), lb = rep(0, times=1000)), aes(x=x,y=y,ymin=lb, ymax=y)) +
  geom_line() + geom_ribbon(alpha = 0.3, linetype= "blank") +
  ggtitle(TeX("$\\lambda \\sim Beta(1, 10)$")) +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

#save_plot(paste0("prior_halfnormal.pdf"), hn.plt, base_width = 3, base_height = 3)
#save_plot(paste0("prior_uniform.pdf"), uni.plt, base_width = 3, base_height = 3)
save_plot(paste0("prior_raised_cosine.pdf"), raised.plt, base_width = 3, base_height = 2)
save_plot(paste0("prior_beta.pdf"), beta.plt, base_width = 3, base_height = 2)

asym.plt <- ggplot(data.frame(x=seq(-3, 3)))