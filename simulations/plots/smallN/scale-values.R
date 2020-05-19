setwd("../")
source("plot_preprocess.R")
setwd("smallN")

coverage <- function(gt,l,h) {
  ifelse (gt >= l & gt <= h, 1, 0)
}

bds.byfn <- smallN.df %>%
  filter(is.element(pos, factor(1:9))) %>% #11 scale values, the first and last are fixed
  mutate(cvc = coverage(gt, ci.low, ci.high)) %>%
  group_by(method, fn, pos, lps, gt, trials) %>%
  summarize(sd=sd(sc, na.rm = TRUE), m=mean(sc, na.rm = TRUE), cv=mean(cvc)) %>%
  mutate(bias = abs(m - gt))

bds.bias.df <- bds.byfn %>%
  group_by(lps, method, trials) %>%
  summarize(bm=mean(bias), b.low=quantile(bias, probs = 0.025), b.high=quantile(bias, probs = 0.975),
            sdm=mean(sd), sd.low=quantile(sd, probs = 0.025), sd.high=quantile(sd, probs = 0.975),
            cvm = mean(cv))

###

bds.bias.plt <- ggplot(bds.bias.df, aes(x=lps, y=bm, colour=method, shape=method)) +
  facet_wrap(~trials) +
  stat_function(fun=function(x) 0, linetype='dashed', colour=solpal5[[5]]) +
  geom_line(aes(x=lps, y=b.low), linetype = "dotted") +
  geom_line(aes(x=lps, y=b.high), linetype = "dotted") +
  geom_ribbon(aes(x=lps, ymin=b.low, ymax=b.high), alpha =0.05, linetype = "blank") +
  geom_line() + geom_point() +
  scale_colour_manual(values = solpal5, labels=c("MLDS", "BDS")) +
  scale_shape_manual(values = c(15, 17), labels=c("MLDS", "BDS")) +
  scale_x_continuous(breaks = c(0, 0.1, 0.2), labels = c("0", "0.1", "0.2")) +
#  ylim(c(0, 0.4)) +
  ylab("mean absolute deviation") +
  xlab("lapse rate") +
  theme(legend.position = "none",
        strip.text = element_text(size=8),
        axis.text = element_text(size=6),
        legend.title = element_text(size=8),
        legend.text = element_text(size=7),
        axis.title = element_text(size=8),
        plot.margin = unit(c(0,0,0,0), "pt"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-10,0,0,0)) +
  NULL

save_plot("smallN_scale_bias.pdf", bds.bias.plt, base_width = 3, base_height = 3)

###

prec.df <- bds.byfn %>%
  group_by(lps, trials, method) %>%
  summarize(sdm=mean(sd), prm=mean(1/sd))

prec.plt <- ggplot(prec.df, aes(x=lps, y=prm, colour=method, shape=method)) +
  geom_point() + geom_line() +
  facet_wrap(~trials) +
  scale_x_continuous(breaks = c(0, 0.1, 0.2), labels = c("0", "0.1", "0.2")) +
  scale_colour_manual(values = solpal5, labels=c("MLDS", "BDS")) +
  scale_shape_manual(values = c(15, 17), labels=c("MLDS", "BDS")) +
  xlab("lapse rate") +
  ylab("average precision") +
  ylim(0, NA) +
  theme(legend.position = "right",#c(0.55, 0.0),
#        legend.justification = c(0.0, 0.0),
        axis.text = element_text(size=6),
        legend.title = element_text(size=8),
        strip.text = element_text(size=8),
        legend.text = element_text(size=7),
        axis.title = element_text(size=8),
        plot.margin = unit(c(0,0,0,0), "pt"),
        legend.margin = margin(0,0,0,0)
        #        legend.box.margin = margin(0,0,0,0)
  ) + 
  NULL

save_plot("smallN_scale_precision.pdf", prec.plt, base_width = 3.2, base_height = 2)

combined.plt <- plot_grid(plotlist=list(bds.bias.plt, prec.plt), ncol=1)
save_plot("smallN_scale_values_combined.pdf", combined.plt, base_width = 3, base_height = 4)

###

bds.coverage.plt <- ggplot(bds.bias.df, aes(x=lps, y=cvm, colour=method, shape=method)) +
  facet_wrap(~trials) +
  stat_function(fun=function(x) 0.95, linetype='dashed', colour=solpal5[[5]]) +
  geom_line() + geom_point() +
  scale_colour_manual(values = solpal5, labels=c("MLDS", "BDS")) +
  scale_shape_manual(values = c(15, 17), labels=c("MLDS", "BDS")) +
  scale_x_continuous(breaks = c(0, 0.1, 0.2), labels = c("0", "0.1", "0.2")) +
  ylab("coverage - scale values") +
  xlab("lapse rate") +
  ylim(0, 1) +
  theme(legend.position = c(1, -0.15),
        legend.justification = c(1, 0.0),
        strip.text = element_text(size=8),
        axis.text = element_text(size=6),
        legend.title = element_blank(),#element_text(size=8),
        legend.text = element_text(size=7),
        axis.title = element_text(size=8),
        axis.title.x = element_blank(),
        plot.margin = unit(c(0,0,3,0), "pt"),
        legend.margin = margin(0,0,10,0),
        legend.box.margin = margin(-10,0,0,0)) +
  NULL

save_plot("smallN_scale_coverage.pdf", bds.coverage.plt, base_width = 3.2, base_height = 2)

###

prec.byfn.df <- bds.byfn %>%
  group_by(lps, trials, fn, method) %>%
  summarize(sdm=mean(sd), prm=mean(1/sd))

prec.byfn.plt <- ggplot(prec.byfn.df, aes(x=lps, y=prm, colour=method)) +
  geom_point() + geom_line() +
  facet_grid(trials~fn) +
  scale_x_continuous(breaks = c(0, 0.1, 0.2), labels = c("0", "0.1", "0.2")) +
  scale_colour_manual(values=solpal5) +
  xlab("lapse rate") +
  ylab("average precision") +
  ylim(0, NA) +
  NULL