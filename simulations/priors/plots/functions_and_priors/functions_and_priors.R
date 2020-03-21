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
    geom_point(data = point.df, size=0.5) +
    ggtitle(fname) +
    scale_y_continuous(breaks= c(0, 1)) +
    scale_x_continuous(breaks= c(0,1)) +
    theme(axis.title=element_blank(),
          plot.margin = margin(0,0,0,0),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          plot.title = element_text(size=8, face='plain')) +
    NULL
  
  plotlist[[fname]] <- fn.plt
}

save_plot("functions.pdf", plot_grid(plotlist = plotlist, nrow = 2), base_height = 1.5, base_width = 3.2)

### prior forms

halfnormal.df <- data.frame(x=c(-2, 0, seq(0, 60, len=1000)), y = c(0, 0, dnorm(seq(0,60, len=1000), sd=20)), lb = rep(0, times=1002), method = "half-normal")
#hn.plt <- ggplot(halfnormal.df, aes(x=precision, y=density, ymin=lb, ymax=density)) +
#  geom_line() +
#  geom_ribbon(alpha = 0.3, linetype= "blank") +
#  ggtitle(TeX("$\\sigma^{-2} \\sim half-normal(0, 20)$")) +
#  theme(axis.title = element_blank(),
#        axis.text.y = element_blank(),
#        axis.ticks.y = element_blank(),
#       axis.line.y = element_blank())

sens.uni.df <- data.frame(x=c(-2, 0, 0, 30, 30, 60), y=1/30*c(0, 0, 1, 1, 0, 0), lb = c(0, 0, 0, 0, 0, 0), method = "uniform(0,30)")
lps.uni.df <- data.frame(x=c(0, 0.2, 1), y=c(1, 1, 1), lb = c(0, 0, 0), method = "uniform(0,1)")
#uni.plt <- ggplot(sens.uni.df, aes(x=x,y=y,ymin=lb, ymax=y)) +
#  geom_line() + geom_ribbon(alpha = 0.3, linetype= "blank") +
#  ggtitle(TeX("$\\Psi_x \\sim Uniform(0,1)$")) +
#  ylim(0, 1.5) +
#  theme(axis.title = element_blank(),
#        axis.text.y = element_blank(),
#        axis.ticks.y = element_blank(),
#        axis.line.y = element_blank())

raised.df <- data.frame(x=c(-2, seq(0, 2.5, len=100), seq(25, 50, len=1000), 60), y=c(0, 0.5*c(-cos(pi*seq(0,2.5, len=100)/2.5), -cos(pi*seq(25, 50, len=1000)/25)) + 0.5, 0), lb = rep(0, times=1102), method = "raised cosine")
raised.plt <- ggplot(raised.df, aes(x=x,y=y,ymin=lb, ymax=y)) +
  geom_line() + geom_ribbon(alpha = 0.3, linetype= "blank") +
#  ggtitle(TeX("$\\sigma^{-1} \\sim raised-cosine(0,2;15,30)$")) +
#  ggtitle("raised-cosine(0, 2.5, 25, 50)") +
  ylim(0, 2.0) +
  xlab("sensitivity") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size=10),
        axis.text.x = element_text(size=7),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        title = element_text(size=12, face="plain"))

beta110.df <- data.frame(x=seq(0,1,len=1000), y = dbeta(seq(0,1, len=1000), 1, 10), lb = rep(0, times=1000), method="Beta(1,10)")
beta15.df <- data.frame(x=seq(0,1,len=1000), y = dbeta(seq(0,1, len=1000), 1, 5), lb = rep(0, times=1000), method="Beta(1,5)")
beta.plt <- ggplot(beta110.df, aes(x=x,y=y,ymin=lb, ymax=y)) +
  geom_line() + geom_ribbon(alpha = 0.3, linetype= "blank") +
#  ggtitle(TeX("$\\lambda \\sim Beta(1, 10)$")) +
#  ggtitle("Beta(1, 10)") +
  xlab('lapse rate') +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size=10),
        axis.text.x = element_text(size=7),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        title = element_text(size=12, face="plain"))

sensitivity.priors.plt <- 
  ggplot(raised.df, aes(x=x, y=y, ymin=lb, ymax=y, colour=method, fill=method)) +
  geom_line() +
  geom_line(data = halfnormal.df) +
  geom_line(data = sens.uni.df) +
  geom_ribbon(alpha = 0.1, linetype= "blank") +
  geom_ribbon(data = halfnormal.df, alpha = 0.1, linetype= "blank") +
  geom_ribbon(data = sens.uni.df, alpha = 0.1, linetype= "blank") +
  scale_colour_manual(values = solpal5) +
  scale_fill_manual(values = solpal5) +
  xlab("sensitivity") +
  ylab("density") +
#  xlim(0, 40) +
  NULL

lps.priors.plt <- 
  ggplot(data.frame(x=c(0, 0.025, 0.025, 0.025, 0.2, 1), y=c(0, 0, 10.0, 0, 0, 0), lb=c(0,0,0,0,0,0), method="fixed(0.025)"), aes(x=x, y=y, ymin=lb, ymax=y, colour=method, fill=method)) +
  geom_line() +
  geom_line(data =lps.uni.df) +
  geom_line(data = beta15.df) +
  geom_line(data = beta110.df) +
  geom_ribbon(alpha = 0.1, linetype= "blank") +
  geom_ribbon(data = lps.uni.df, alpha = 0.1, linetype= "blank") +
  geom_ribbon(data = beta15.df, alpha = 0.1, linetype= "blank") +
  geom_ribbon(data = beta110.df, alpha = 0.1, linetype= "blank") +
  scale_colour_manual(values = solpal5) +
  scale_fill_manual(values = solpal5) +
  xlab("lapse rate") +
  ylab("density") +
  #  xlim(0, 40) +
  NULL

#save_plot(paste0("prior_halfnormal.pdf"), hn.plt, base_width = 3, base_height = 3)
#save_plot(paste0("prior_uniform.pdf"), uni.plt, base_width = 3, base_height = 3)
save_plot(paste0("prior_raised_cosine.pdf"), raised.plt, base_width = 3, base_height = 2)
save_plot(paste0("prior_beta.pdf"), beta.plt, base_width = 3, base_height = 2)

priors.plt <- plot_grid(raised.plt, beta.plt, labels=c("A", "B"))
save_plot(paste0("priors.pdf"), priors.plt, base_width = 3.2, base_height = 2)

save_plot("sensitivity_priors.pdf", sensitivity.priors.plt, base_width = 6, base_height = 2.5)
save_plot("lapserate_priors.pdf", lps.priors.plt, base_width = 6, base_height = 2.5)
save_plot("lapserate_priors_zoomed.pdf", lps.priors.plt + xlim(0, 0.2), base_width = 6, base_height = 2.5)
asym.plt <- ggplot(data.frame(x=seq(-3, 3)))

#

xval <- seq(0, 1, len=5)
yval <- sqrt(xval)
y.min <- 0.8 * sqrt(xval)
y.max <- 1.2 * sqrt(xval)
fn.sqrt <- function(x) 2 * sqrt(x)

dat <- data.frame(xval=xval, yval=yval, ymin=y.min, ymax=y.max, mt='estimate')

sens.err.plt <- ggplot(dat, aes(x=xval,y=yval,ymin=ymin,ymax=ymax,colour=mt)) +
  stat_function(fun=fn.sqrt, aes(colour='ground truth')) +
  geom_line() +
  geom_ribbon(alpha=0.3, linetype='blank') +
  scale_y_continuous(breaks= c(0, 2), labels = c('0', 'S')) +
  scale_x_continuous(breaks= c(0, 1)) +
  scale_colour_manual(values = c(solarized[['violet']], 'black')) +
  scale_fill_manual(values = solpal5) +
  geom_point(data=data.frame(xval=c(1.0, 1.0),ymin=c(0.8, 1.6), ymax=c(1.2, 2.4), yval=c(1.0, 2.0), mt=c('estimate', 'estimate')), size=2) +
  geom_errorbar(data=data.frame(xval=c(1.0, 1.0),ymin=c(0.8, 1.6), ymax=c(1.2, 2.4), yval=c(1.0, 2.0), mt=c('estimate', 'estimate')), width=.05) +
  theme(legend.position = 'none',
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "pt"),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6)) +
#        axis.ticks.y = element_blank(),
#        axis.line.y = element_blank())
  NULL

save_plot("sensitivity_error_metrics.svg", sens.err.plt, base_width = 1.5, base_height = 1.5)

#

xval <- seq(0, 1, len=5)
yval <- xval
y.min <- xval^1.2
y.max <- xval+(xval - y.min)
fn.sqrt <- function(x) sqrt(x)

dat <- data.frame(xval=xval, yval=yval, ymin=y.min, ymax=y.max, mt='estimate')

sc.err.plt <- ggplot(dat, aes(x=xval,y=yval,ymin=ymin,ymax=ymax,colour=mt)) +
  stat_function(fun=fn.sqrt, aes(colour='ground truth')) +
  geom_line() +
  geom_point(size=2) +
  geom_errorbar(width=0.05) +
  scale_y_continuous(breaks= c(0, 1)) +
  scale_x_continuous(breaks= c(0, 1)) +
  scale_colour_manual(values = c(solarized[['violet']], 'black')) +
  scale_fill_manual(values = solpal5) +
#  geom_point(data=data.frame(xval=c(1.0, 1.0),ymin=c(0.8, 1.6), ymax=c(1.2, 2.4), yval=c(1.0, 2.0), mt=c('estimate', 'estimate')), size=2) +
#  geom_errorbar(data=data.frame(xval=c(1.0, 1.0),ymin=c(0.8, 1.6), ymax=c(1.2, 2.4), yval=c(1.0, 2.0), mt=c('estimate', 'estimate')), width=.05) +
  theme(legend.position = 'none',
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "pt"),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6)) +
  #        axis.ticks.y = element_blank(),
  #        axis.line.y = element_blank())
  NULL
sc.err.plt

save_plot("scale_value_error_metrics.svg", sc.err.plt, base_width = 1.5, base_height = 1.5)
