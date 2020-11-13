setwd("../")
source("plot_preprocess.R")
setwd("scale-values")

coverage <- function(gt,l,h) {
  ifelse (gt >= l & gt <= h, 1, 0)
}

bds.byfn <- data.df %>%
  filter(is.element(pos, factor(1:9)), is.element(lps, c(0, 0.005, 0.01, 0.03, 0.05, 0.1))) %>% #11 scale values, the first and last are fixed
  mutate(cvc = coverage(gt, ci.low, ci.high)) %>%
  group_by(method, fn, pos, lps, gt, trials) %>%
  summarize(sd=sd(sc, na.rm = TRUE), m=mean(sc, na.rm = TRUE), cv=mean(cvc)) %>%
  mutate(bias = abs(m - gt))

bds.byfn$lps <- as.factor(bds.byfn$lps)

bds.bias.df <- bds.byfn %>%
  group_by(lps, method, trials) %>%
  summarize(bm=mean(bias), b.low=quantile(bias, probs = 0.025), b.high=quantile(bias, probs = 0.975),
            sdm=mean(sd), sd.low=quantile(sd, probs = 0.025), sd.high=quantile(sd, probs = 0.975),
            cvm = mean(cv))

bds.bias.df$lps <- as.factor(bds.bias.df$lps)
###

max_bias <- bds.byfn %>% group_by(lps, method, trials) %>% summarize(mb=max(bias))

bds.bias.plt <- ggplot(bds.byfn, aes(x=lps, y=bias, colour=method, shape=method)) +
  facet_wrap(~trials, ncol=1) +
#  stat_function(fun=function(x) 0, linetype='dashed', colour=solpal5[[5]]) +
  geom_violin(aes(fill=method), scale='width', position = position_dodge(0.75), lwd=0.1, colour='black') +
  geom_boxplot(position = position_dodge(0.75), width=0.2, outlier.shape = 32, lwd=0.2, coef = 0, colour='black', fill='white') +
#  geom_point(aes(y=mb), max_bias, position = position_dodge(1), shape=4) + 
  scale_colour_manual(values = solpal5, labels=c("MLDS", "BDS")) +
  scale_fill_manual(values = solpal5, labels=c("MLDS", "BDS")) +
  scale_shape_manual(values = c(15, 17), labels=c("MLDS", "BDS")) +
#  scale_x_discrete(breaks = c(0, 0.1, 0.2), labels = c("0", "0.1", "0.2")) +
  ylab("absolute deviation") +
  xlab("lapse rate") +
  theme(legend.position = "none",#c(0.5, .90),
#        legend.justification = c(0.0, 1.0),
        strip.text = element_blank(),#element_text(size=8),
        strip.background = element_blank(),
        axis.text = element_text(size=6),
#        legend.title = element_blank(),#element_text(size=8),
#        legend.text = element_text(size=7),
        axis.title.y = element_text(size=8),
        axis.title.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "pt"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-10,0,0,0)) +
  NULL

save_plot("bds_scale_bias.pdf", bds.bias.plt, base_width = 3, base_height = 3)

###

prec.df <- bds.byfn %>%
  group_by(lps, trials, method) %>%
  summarize(sdm=mean(sd), prm=mean(1/sd))

prec.plt <- ggplot(bds.byfn, aes(x=lps, y=1/sd, colour=method, shape=method)) +
  geom_violin(aes(fill=method), scale='width', position = position_dodge(0.75), lwd=0.1, colour='black') +
  geom_boxplot(position = position_dodge(0.75), width=0.2, outlier.shape = 32, lwd=0.2, coef = 0, colour='black', fill='white') +
  facet_wrap(~trials, ncol=1) +
#  scale_x_continuous(breaks = c(0, 0.1, 0.2), labels = c("0", "0.1", "0.2")) +
  scale_colour_manual(values = solpal5, labels=c("MLDS", "BDS")) +
  scale_shape_manual(values = c(15, 17), labels=c("MLDS", "BDS")) +
  scale_fill_manual(values = solpal5, labels=c("MLDS", "BDS")) +
  xlab("lapse rate") +
  ylab("precision") +
  ylim(0, NA) +
  theme(legend.position = c(0.5, 1.0),
        legend.justification = c(1.0, 1.0),
        axis.text = element_text(size=6),
        legend.title = element_blank(),#element_text(size=8),
        strip.text = element_blank(),#element_text(size=8),
        strip.background = element_blank(),
        legend.text = element_text(size=7),
        axis.title = element_text(size=8),
        plot.margin = unit(c(0,0,0,0), "pt"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(0,0,0,0)
  ) + 
  NULL

save_plot("bds_scale_precision.pdf", prec.plt, base_width = 3, base_height = 3)

combined.plt <- plot_grid(plotlist=list(bds.bias.plt, prec.plt), ncol=1)
save_plot("bds_scale_values_combined.pdf", combined.plt, base_width = 3, base_height = 5)
###

bds.coverage.plt <- ggplot(bds.byfn, aes(x=lps, y=cv, colour=method, shape=method)) +
  facet_wrap(~trials, ncol=1) +
  stat_function(fun=function(x) 0.95, linetype='dashed', colour=solpal5[[5]]) +
  geom_violin(aes(fill=method), scale='width', position = position_dodge(0.75), lwd=0.1, colour='black') +
  geom_boxplot(position = position_dodge(0.75), width=0.2, outlier.shape = 32, lwd=0.2, coef = 0, colour='black', fill='white') +
  scale_colour_manual(values = solpal5, labels=c("MLDS", "BDS")) +
  scale_fill_manual(values = solpal5, labels=c("MLDS", "BDS")) +
  scale_shape_manual(values = c(15, 17), labels=c("MLDS", "BDS")) +
#  scale_x_continuous(breaks = c(0, 0.1, 0.2), labels = c("0", "0.1", "0.2")) +
  ylab("coverage - scale values") +
  xlab("lapse rate") +
#  ylim(0, 1) +
  theme(legend.position = c(1, 0.0),
        legend.justification = c(1.0, 0.0),
        strip.text = element_blank(),#element_text(size=8),
        strip.background = element_blank(),
        axis.text = element_text(size=6),
        legend.title = element_blank(),#element_text(size=8),
        legend.text = element_text(size=7),
        axis.title.y = element_text(size=8),
        axis.title.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "pt"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-10,0,0,0)) +
  NULL

save_plot("bds_scale_coverage.pdf", bds.coverage.plt, base_width = 3, base_height = 2)

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

stimulus <- c(0, seq(0.025, 0.975, len=9), 1)
names(stimulus) <- sapply(0:10, toString)

fn.df <- data.frame(stim=numeric(), m=numeric(), sd=numeric(), fn=factor(), method=factor())

for (fname in names(function.zoo)) {
  x <- seq(0, 1, len=1000)
  y <- function.zoo[[fname]](x)
  fn.df <- rbind(fn.df, data.frame(stim=x, m=y, sd=0, fn=fname, method="z"))
}



fn.data.df <- data.df %>% filter(trials == 330, lps == 0, sens==10) %>% filter(is.element(pos, factor(0:10))) %>% group_by(pos, method, fn) %>% summarize(m=mean(sc), sd=sd(sc))
fn.data.df$stim <- sapply(fn.data.df$pos, function(x) stimulus[[toString(x)]])

fn.plt <- ggplot(fn.data.df, aes(x=stim, y=m, ymin=m-sd, ymax=m+sd, colour=method, shape=method)) +
  facet_wrap(~fn, ncol=5) +
  geom_line(aes(x=stim, y=m), data=fn.df, linetype="dashed") +
  geom_pointrange() +
  scale_colour_manual(values = c(solpal5[[1]], solpal5[[2]], "grey"), labels=c("MLDS", "BDS", "ground truth")) +
  scale_shape_manual(values = c(15, 17, 32), labels=c("MLDS", "BDS", "ground truth")) +
  scale_y_continuous(breaks= c(0, 1)) +
  scale_x_continuous(breaks= c(0,1)) +
  xlab("stimulus") +
  ylab("scale") +
  theme(legend.position = "none",
        plot.margin = margin(0,0,0,0),
        strip.text = element_text(size=8))

save_plot("example_scales.pdf", fn.plt, base_width = 7, base_height = 3)

fn2.data.df <- data.df %>% filter(trials == 330, lps == 0.05, sens==10) %>% filter(is.element(pos, factor(0:10))) %>% group_by(pos, method, fn) %>% summarize(m=mean(sc), sd=sd(sc))
fn2.data.df$stim <- sapply(fn.data.df$pos, function(x) stimulus[[toString(x)]])

fn2.plt <- ggplot(fn2.data.df, aes(x=stim, y=m, ymin=m-sd, ymax=m+sd, colour=method, shape=method)) +
  facet_wrap(~fn, ncol=5) +
  geom_line(aes(x=stim, y=m), data=fn.df, linetype="dashed") +
  geom_pointrange() +
  scale_colour_manual(values = c(solpal5[[1]], solpal5[[2]], "grey"), labels=c("MLDS", "BDS", "ground truth")) +
  scale_shape_manual(values = c(15, 17, 32), labels=c("MLDS", "BDS", "ground truth")) +
  scale_y_continuous(breaks= c(0, 1)) +
  scale_x_continuous(breaks= c(0,1)) +
  xlab("stimulus") +
  ylab("scale") +
  theme(legend.position = "none",
        plot.margin = margin(0,0,0,0),
        strip.text = element_text(size=8))

save_plot("example_scales_lapses.pdf", fn2.plt, base_width = 7, base_height = 3)

##
#
# Trying to test the influences of divergences in HMC
data.wide <- dcast(data.df, lvl + sens + trials + lps + fn + simid ~ method + pos, value.var = "sc")

data.wide.nodiv <- data.wide %>% filter(mixture_divergent < 0.03, lps < 0.1)

data.intersect <- intersect(data.df[c("lvl", "sens", "trials", "lps", "fn", "simid")], data.wide.nodiv[c("lvl", "sens", "trials", "lps", "fn", "simid")])
