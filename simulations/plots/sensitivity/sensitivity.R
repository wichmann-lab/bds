setwd("../")
source("plot_preprocess.R")
setwd("sensitivity")

coverage <- function(gt,l,h) {
  ifelse (gt >= l & gt <= h, 1, 0)
}

sensitivity.df <- data.df %>% filter(pos == "sigma") %>%
  mutate(cvc = coverage(gt, ci.low, ci.high)) %>%
  group_by(method, trials, lps) %>% summarize(pb = mean(sc/gt),
           low = quantile(sc/gt, probs=0.025)[[1]],
           high = quantile(sc/gt, probs=0.975)[[1]],
           prec = 1/sd(sc/gt),
           nprec = 1/(sd(sc)/mean(sc)),
           cv = mean(cvc))

sensitivity.byfn.df <- data.df %>% filter(pos == "sigma" & trials == 330) %>% group_by(method, lps, fn, sens) %>% summarize(pb = mean(sc)/mean(gt) - 1,
                                                                                              low = quantile(sc, probs=0.025)[[1]]/mean(gt) - 1,
                                                                                              high = quantile(sc, probs=0.975)[[1]]/mean(gt) - 1)

bds.sens.plt <- ggplot(sensitivity.df, aes(x=lps, y=pb, colour=method, shape=method)) +
  stat_function(fun=function(x) 1.0, linetype='dashed', colour=solpal5[[5]]) +
  geom_line(aes(x=lps, y=low), linetype = "dotted") +
  geom_line(aes(x=lps, y=high), linetype = "dotted") +
  geom_ribbon(aes(x=lps, ymin=low, ymax=high), alpha =0.05, linetype = "blank") +
  facet_wrap(~trials) +
  scale_colour_manual(values = solpal5, labels=c("MLDS", "BDS")) +
  scale_shape_manual(values = c(15, 17), labels=c("MLDS", "BDS")) +
  geom_line() +
  geom_point() +
  xlab("lapse rate") +
  ylab("normalised sensitivity estimate") +
  scale_x_continuous(breaks = c(0, 0.1, 0.2), labels = c("0", "0.1", "0.2")) +
  theme(legend.position = c(1.0, 0.95),
        legend.justification = c(1, 1),
        axis.text = element_text(size=6),
        legend.title = element_text(size=8),
        strip.text = element_text(size=8),
        legend.text = element_text(size=7),
        axis.title = element_text(size=8),
        plot.margin = unit(c(0,0,0,0), "pt"),
        legend.margin = margin(0,0,0,0)) +
  NULL

save_plot("bds_sim_sensitivity_bias.pdf", bds.sens.plt, base_width = 3, base_height=2)

###

bds.sens.prec.plt <- ggplot(sensitivity.df, aes(x=lps, y=prec, colour=method, shape=method)) +
  facet_wrap(~trials) +
  geom_line() +
  geom_point() +
  xlab("lapse rate") +
  ylab("normalised precision") +
  scale_colour_manual(values = solpal5, labels=c("MLDS", "BDS")) +
  scale_shape_manual(values = c(15, 17), labels=c("MLDS", "BDS")) +
  scale_x_continuous(breaks = c(0, 0.1, 0.2), labels = c("0", "0.1", "0.2")) +
  ylim(0, NA) +
  theme(legend.position = c(1,0), legend.justification = c(1,0),
        axis.text = element_text(size=6),
        legend.title = element_text(size=8),
        strip.text = element_text(size=8),
        legend.text = element_text(size=7),
        axis.title = element_text(size=8),
        plot.margin = unit(c(0,0,0,0), "pt"),
        legend.margin = margin(0,0,0,0)) +
  NULL

save_plot("bds_sim_sensitivity_precision.pdf", bds.sens.prec.plt, base_width = 3, base_height=2)

###

bds.sens.coverage.plt <- ggplot(sensitivity.df, aes(x=lps, y=cv, colour=method, shape=method)) +
  stat_function(fun=function(x) 0.95, linetype='dashed', colour=solpal5[[5]]) +
  facet_wrap(~trials) +
  geom_line() +
  geom_point() +
  xlab("lapse rate") +
  ylab("coverage") +
  scale_colour_manual(values = solpal5, labels=c("MLDS", "BDS")) +
  scale_shape_manual(values = c(15, 17), labels=c("MLDS", "BDS")) +
  scale_x_continuous(breaks = c(0, 0.1, 0.2), labels = c("0", "0.1", "0.2")) +
  ylim(0, 1) +
  theme(legend.position = c(1,.95), legend.justification = c(1,1),
        axis.text = element_text(size=6),
        legend.title = element_text(size=8),
        strip.text = element_text(size=8),
        legend.text = element_text(size=7),
        axis.title = element_text(size=8),
        plot.margin = unit(c(0,0,0,0), "pt"),
        legend.margin = margin(0,0,0,0)) +
  NULL

save_plot("bds_sim_sensitivity_coverage.pdf", bds.sens.coverage.plt, base_width = 3, base_height=2)

###

ggplot(sensitivity.byfn.df, aes(x=sens, y=pb, colour=method)) +
  facet_grid(lps~fn) +
  stat_function(fun=function(x) 0, linetype='dashed', colour=solpal5[[5]]) +
#  geom_errorbar(aes(x=sens, ymin=low, ymax=high)) +
  scale_colour_manual(values = solpal5) +
  geom_point() +
  geom_line() +
  NULL