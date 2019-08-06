setwd("../../")
source("plot_preprocess.R")
setwd("plots/lapserate")

sens.df <- lps.priors.df %>% filter(pos == "sigma") %>% group_by(method, trials, lps, prec) %>% summarize(pb = mean(sc/gt) - 1,
                                                                                      low = quantile(sc/gt, probs=0.025)[[1]] - 1,
                                                                                      high = quantile(sc/gt, probs=0.975)[[1]] - 1)

byfn.df <- lps.priors.df %>% filter(pos == "sigma", lps == 0.0) %>% group_by(method, trials, fn, prec) %>% summarize(pb = mean(sc)/mean(gt) - 1,
                                                                                              low = quantile(sc, probs=0.025)[[1]]/mean(gt) - 1,
                                                                                              high = quantile(sc, probs=0.975)[[1]]/mean(gt) - 1)

ggplot(sens.df, aes(x=lps, y=pb, colour=method)) +
  facet_grid(trials~prec) +
  scale_colour_manual(values = solpal5) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(x=lps, ymin=low, ymax=high)) +
  NULL

ggplot(byfn.df, aes(x=prec, y=pb, colour=method)) +
  facet_grid(trials~fn) +
  scale_colour_manual(values = solpal5) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(x=prec, ymin=low, ymax=high)) +
  NULL


###

gtfn <- function(fn, gt) function.zoo[[as.character(fn)]](gt)
gtfnv <- function(fnv, gtv) mapply(gtfn, fnv, gtv, SIMPLIFY = TRUE)

lapses.byfn <- all.df %>%
  filter(is.element(pos, factor(1:9))) %>%
  group_by(method, fn, pos, lps, gt, trials) %>%
  summarize(sd=mad(sc, na.rm = TRUE), m=median(sc, na.rm = TRUE)) %>%
  mutate(bias = abs(m - function.zoo[[as.character(fn)]](gt)))

all.bias.df <- lapses.byfn %>%
  group_by(lps, method) %>%
  summarize(bm=mean(bias))

all.bias.plt <- ggplot(all.bias.df, aes(x=lps, y=bm, colour=method)) +
  geom_point() + geom_line() +
  #  geom_errorbar(aes(x=lps, ymin=bias-sd, ymax=bias+sd)) +
  scale_colour_manual(values= solpal5) +
  ylab("absolute bias") +
  xlab("lapse rate") +
  NULL

save_plot(paste0("plots/", "avg_scale_bias.pdf"), all.bias.plt, base_aspect_ratio = 1.2, base_height = 8)

grouped.bias.df <- lapses.byfn %>%
  group_by(lps, method, fn) %>%
  summarize(bm=mean(bias))

grouped.bias.plt <- ggplot(grouped.bias.df, aes(x=lps, y=bm, colour=method)) +
  geom_point(size=2) + geom_line() +
  facet_wrap(~fn, nrow=2) +
  #  coord_cartesian(ylim=c(0,1)) +
  #  geom_errorbar(aes(x=lps, ymin=bias-sd, ymax=bias+sd)) +
  scale_colour_manual(values= solpal5) +
  ylab("average bias") +
  xlab("lapse rate") +
  theme(#legend.justification=c(1,0), legend.position=c(1,.2),
    legend.background = element_rect(colour = 'black', linetype = 'solid'),
    text = element_text(size=16),
    strip.text = element_text(size=16),
    axis.text = element_text(size=14),
    axis.title = element_text(size=18)) +
  NULL

prec.byfn.df <- lapses.byfn %>%
  group_by(lps, trials, fn, method) %>%
  summarize(sdm=mean(sd))

prec.byfn.plt <- ggplot(prec.byfn.df, aes(x=lps, y=1/sdm, colour=factor(trials))) +
  geom_point() + geom_line() +
  facet_grid(method~fn) +
  scale_x_continuous(breaks = c(0, 0.1, 0.2), labels = c("0", "0.1", "0.2")) +
  scale_colour_manual(values=solpal5) +
  xlab("lapse rate") +
  ylab("average precision") +
  ylim(0, NA) +
  NULL 
