setwd("../../")
source("plot_preprocess.R")
setwd("plots/sensitivity")

sensitivity.df <- sensitivity.priors.df %>% filter(pos == "sigma") %>% group_by(method, trials, prec) %>% summarize(pb = mean(sc)/mean(gt) - 1,
                                                                                      low = quantile(sc, probs=0.025)[[1]]/mean(gt) - 1,
                                                                                      high = quantile(sc, probs=0.975)[[1]]/mean(gt) - 1)

sensitivity.bias.plt <- ggplot(sensitivity.df, aes(x=prec, y=pb, colour=method)) +
  stat_function(fun=function(x) 0, linetype='dashed', colour=solpal5[[5]]) +
  facet_wrap(~trials) +
  scale_colour_manual(values = solpal5) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(x=prec, ymin=low, ymax=high)) +
  xlab("sensitivity") +
  ylab("normalized bias of estimate") +
  NULL

save_plot("sens_prior_sim_bias.pdf", sensitivity.bias.plt, base_aspect_ratio = 1.2, base_height = 6)
###

gtfn <- function(fn, gt) function.zoo[[as.character(fn)]](gt)
gtfnv <- function(fnv, gtv) mapply(gtfn, fnv, gtv, SIMPLIFY = TRUE)

lapses.byfn <- sensitivity.priors.df %>%
  filter(is.element(pos, factor(1:9))) %>%
  group_by(method, fn, pos, gt, trials) %>%
  summarize(sd=sd(sc, na.rm = TRUE), m=mean(sc, na.rm = TRUE)) %>%
  mutate(bias = abs(m - gtfnv(fn, gt)))

all.bias.df <- lapses.byfn %>%
  group_by(trials, method) %>%
  summarize(bm=mean(bias), b.low=quantile(bias, probs = 0.025), b.high=quantile(bias, probs = 0.975),
            sdm=mean(sd), sd.low=quantile(sd, probs = 0.025), sd.high=quantile(sd, probs = 0.975))

all.bias.plt <- ggplot(all.bias.df, aes(x=trials, y=bm, colour=method)) +
  geom_point() + geom_line() +
  geom_errorbar(aes(x=trials, ymin=b_low, ymax=b_high)) +
  scale_colour_manual(values= solpal5) +
  ylab("absolute bias") +
  xlab("#trials") +
  NULL

save_plot("sens_prior_sim_scale_bias.pdf", all.bias.plt, base_aspect_ratio = 1.2, base_height = 6)

grouped.bias.df <- lapses.byfn %>%
  group_by(trials, method, fn) %>%
  summarize(bm=mean(bias), b.low=quantile(bias, probs=0.025), b.high=quantile(bias, probs=0.975))

grouped.bias.plt <- ggplot(grouped.bias.df, aes(x=trials, y=bm, colour=method)) +
  geom_point(size=2) + geom_line() +
  facet_wrap(~fn, nrow=2) +
  #  coord_cartesian(ylim=c(0,1)) +
  geom_errorbar(aes(x=trials, ymin=b.low, ymax=b.high)) +
  scale_colour_manual(values= solpal5) +
  ylab("average bias") +
  xlab("#trials") +
  theme(#legend.justification=c(1,0), legend.position=c(1,.2),
    legend.background = element_rect(colour = 'black', linetype = 'solid'),
    text = element_text(size=16),
    strip.text = element_text(size=16),
    axis.text = element_text(size=14),
    axis.title = element_text(size=18)) +
  NULL

save_plot("sens_prior_sim_scale_bias_byfn.pdf", grouped.bias.plt, base_aspect_ratio = 1.5, base_height = 6)

prec.plt <- ggplot(all.bias.df, aes(x=trials, y=1/sdm, ymin=1/sd.low, ymax=1/sd.high, colour=method)) +
  stat_smooth(aes(y=1/sdm), method = "lm", formula = y ~ I(sqrt(x)), se = F, fullrange = T) +
#  stat_smooth(aes(y=1/sd.low), method = "lm", formula = y ~ I(sqrt(x)) + x, se = F, fullrange = T) +
#  stat_smooth(aes(y=1/sd.high), method = "lm", formula = y ~ I(sqrt(x)) + x, se = F, fullrange = T) +
  geom_point() + geom_errorbar() +
#  geom_line() +
  scale_colour_manual(values = solpal5) +
  ylim(0, NA) +
  xlim(0, 1500) +
  ylab("average precision") +
  xlab("#trials") +
  NULL

save_plot("sens_prior_sim_prec.pdf", prec.plt, base_aspect_ratio = 1.2, base_height = 6)

prec.byfn.df <- lapses.byfn %>%
  group_by(trials, fn, method) %>%
  summarize(sdm=mean(sd), sd.low = quantile(sd, probs=0.025), sd.high = quantile(sd, probs = 0.975))

prec.byfn.plt <- ggplot(prec.byfn.df, aes(x=trials, y=1/sdm, colour=method)) +
  geom_point() + geom_line() +
  geom_errorbar(aes(ymin=1/sd.low, ymax=1/sd.high)) +
  facet_wrap(~fn, nrow = 2) +
  scale_x_continuous(breaks = c(330, 660, 990, 1320), labels = c("330", "", "990", "")) +
  scale_colour_manual(values=solpal5) +
  xlab("lapse rate") +
  ylab("average precision") +
  ylim(0, NA) +
  NULL 

save_plot("sens_prior_sim_prec_byfn.pdf", grouped.bias.plt, base_aspect_ratio = 1.5, base_height = 6)