setwd("../")
source("plot_preprocess.R")
setwd("sensitivity")

coverage <- function(gt,l,h) {
  ifelse (gt >= l & gt <= h, 1, 0)
}

sensitivity.df <- data.df %>% filter(pos == "sigma") %>% filter(is.element(lps, c(0, 0.005, 0.01, 0.03, 0.05, 0.1))) %>%
  mutate(cvc = coverage(gt, ci.low, ci.high)) %>%
  group_by(method, trials, lps) %>% summarize(pb = mean(sc/gt),
           low = quantile(sc/gt, probs=0.025)[[1]],
           high = quantile(sc/gt, probs=0.975)[[1]],
           prec = 1/sd(sc/gt),
           nprec = 1/(sd(sc)/mean(sc)),
           cv = mean(cvc))

sensitivity.byfn.df <- data.df %>% filter(pos == "sigma") %>% filter(is.element(lps, c(0, 0.005, 0.01, 0.03, 0.05, 0.1))) %>%
  mutate(cvc = coverage(gt, ci.low, ci.high)) %>%
  group_by(method, trials, lps, fn) %>%
  summarize(pb = mean(sc/gt),
  prec = 1/sd(sc/gt),
  cv = mean(cvc),
  low = quantile(sc, probs=0.025)[[1]]/mean(gt) - 1,
  high = quantile(sc, probs=0.975)[[1]]/mean(gt) - 1)

sensitivity.byfn.df$lps <- as.factor(sensitivity.byfn.df$lps)

data_summary <- function(x) {
  m <- median(x)
  ymin <- quantile(x, probs = c(0.25))[[1]]
  ymax <- quantile(x, probs=c(0.75))[[1]]
  return(c(y=m,ymin=ymin,ymax=ymax))
}

bds.sens.plt <- ggplot(sensitivity.byfn.df, aes(x=lps, y=pb, colour=method, shape=method)) +
  stat_function(fun=function(x) 1.0, linetype='dashed', colour=solpal5[[5]]) +
  geom_violin(aes(fill=method), scale='width', position = position_dodge(0.75), lwd=0.1, colour='black') +
  geom_boxplot(position = position_dodge(0.75), width=0.2, outlier.shape = 32, lwd=0.2, coef = 0, colour='black', fill='white') +
  facet_wrap(~trials, ncol=1) +
  scale_colour_manual(values = solpal5, labels=c("MLDS", "BDS")) +
  scale_fill_manual(values = solpal5, labels=c("MLDS", "BDS")) +
  scale_shape_manual(values = c(15, 17), labels=c("MLDS", "BDS")) +
  xlab("lapse rate") +
  ylab("normalised sensitivity estimate") +
  theme(legend.position = "none",#c(1.0, 0.95),
#        legend.justification = c(1, 1),
        axis.text = element_text(size=6),
#        legend.title = element_text(size=8),
        strip.text = element_blank(),
        strip.background = element_blank(),
#        legend.text = element_text(size=7),
        axis.title.y = element_text(size=8),
        axis.title.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "pt"),
        legend.margin = margin(0,0,0,0)) +
  NULL

save_plot("bds_sim_sensitivity_bias.pdf", bds.sens.plt, base_width = 3, base_height=3)

###

bds.sens.prec.plt <- ggplot(sensitivity.byfn.df, aes(x=lps, y=prec, colour=method, shape=method)) +
  facet_wrap(~trials, ncol=1) +
  geom_violin(aes(fill=method), scale='width', position = position_dodge(0.75), lwd=0.1, colour='black') +
  geom_boxplot(position = position_dodge(0.75), width=0.2, outlier.shape = 32, lwd=0.2, coef = 0, colour='black', fill='white') +
  xlab("lapse rate") +
  ylab("normalised precision") +
  scale_colour_manual(values = solpal5, labels=c("MLDS", "BDS")) +
  scale_shape_manual(values = c(15, 17), labels=c("MLDS", "BDS")) +
  scale_fill_manual(values = solpal5, labels=c("MLDS", "BDS")) +
#  scale_x_continuous(breaks = c(0, 0.1, 0.2), labels = c("0", "0.1", "0.2")) +
  ylim(0, NA) +
  theme(legend.position = c(1,0), legend.justification = c(1,0),
        axis.text = element_text(size=6),
        legend.title = element_blank(),#element_text(size=8),
        strip.text = element_blank(), #element_text(size=8),
        strip.background = element_blank(),
        legend.text = element_text(size=7),
        axis.title = element_text(size=8),
        plot.margin = unit(c(0,0,0,0), "pt"),
        legend.margin = margin(0,0,0,0)) +
  NULL

save_plot("bds_sim_sensitivity_precision.pdf", bds.sens.prec.plt, base_width = 3, base_height=3)

combined.plt <- plot_grid(plotlist=list(bds.sens.plt, bds.sens.prec.plt), ncol=1)
save_plot("bds_sensitivity_combined.pdf", combined.plt, base_width = 3, base_height = 5)
###
###

bds.sens.coverage.plt <- ggplot(sensitivity.byfn.df, aes(x=lps, y=cv, colour=method, shape=method)) +
  facet_wrap(~trials, ncol=1) +
  stat_function(fun=function(x) 0.95, linetype='dashed', colour=solpal5[[5]]) +
  geom_violin(aes(fill=method), scale='width', position = position_dodge(0.75), lwd=0.1, colour='black') +
  geom_boxplot(position = position_dodge(0.75), width=0.2, outlier.shape = 32, lwd=0.2, coef = 0, colour='black', fill='white') +
  scale_colour_manual(values = solpal5, labels=c("MLDS", "BDS")) +
  scale_fill_manual(values = solpal5, labels=c("MLDS", "BDS")) +
  scale_shape_manual(values = c(15, 17), labels=c("MLDS", "BDS")) +
  #  scale_x_continuous(breaks = c(0, 0.1, 0.2), labels = c("0", "0.1", "0.2")) +
  ylab("coverage - sensitivity") +
  xlab("lapse rate") +
  ylim(0, 1) +
  theme(legend.position = "none",#c(1,.95), legend.justification = c(1,1),
        axis.text = element_text(size=6),
        legend.title = element_blank(),#element_text(size=8),
        strip.text = element_blank(), #element_text(size=8),
        strip.background = element_blank(),
        legend.text = element_text(size=7),
        axis.title = element_text(size=8),
        plot.margin = unit(c(0,0,0,0), "pt"),
        legend.margin = margin(0,0,0,0)) +
  NULL

save_plot("bds_sim_sensitivity_coverage.pdf", bds.sens.coverage.plt, base_width = 3, base_height=2)
coverage.combined.plt <- plot_grid(plotlist=list(bds.coverage.plt, bds.sens.coverage.plt), ncol=1)
save_plot("bds_coverage_combined.pdf", coverage.combined.plt, base_width = 3, base_height = 5)
###

ggplot(sensitivity.byfn.df, aes(x=sens, y=pb, colour=method)) +
  facet_grid(lps~fn) +
  stat_function(fun=function(x) 0, linetype='dashed', colour=solpal5[[5]]) +
#  geom_errorbar(aes(x=sens, ymin=low, ymax=high)) +
  scale_colour_manual(values = solpal5) +
  geom_point() +
  geom_line() +
  NULL