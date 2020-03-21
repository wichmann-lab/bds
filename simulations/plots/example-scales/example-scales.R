setwd("../")
source("plot_preprocess.R")
setwd("sensitivity")

dat <- read.csv('../../data/bds-sqrt-11-330-10-0.05-sim.csv', sep = '\t')
dat[['id']] <- c(as.vector(sapply(1:144, rep, times=17)), as.vector(sapply(1:144, rep, times=19)))

stimulus <- c(0, seq(0.025, 0.975, len=9), 1)
names(stimulus) <- sapply(0:10, toString)

datcast <- dcast(dat, id + method ~ pos, value.var = 'sc')
datcast[sapply(1:10, toString)] <- datcast[sapply(1:10, toString)] * datcast[['sigma']]

datmelt <- melt(datcast, id.vars = c('id', 'method'), measure.vars = sapply(0:10, toString))

single_scale <- datmelt %>% filter(id == 97)
single_scale$stim <- sapply(single_scale$variable, function(x) stimulus[[toString(x)]])

single_scales.plt <- ggplot(single_scale, aes(x=stim, y=value, colour=method, shape=method)) + geom_line() + geom_point() +
  scale_colour_manual(values=solpal5[c(1,2)]) +
  scale_shape_manual(values = c(15, 17)) +
  xlab("stimulus") + ylab("estimated scale") +
  scale_y_continuous(limits = c(0, max(single_scale$value)), breaks = c(10)) +
  theme(legend.justification=c(1,0), legend.position=c(1,0),
        axis.text.y = element_blank(),#element_text(size=6),
        axis.title.x = element_blank(),#element_text(size=18),
        axis.title.y = element_blank(),#element_text(size=10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size=7),
        #            legend.background = element_rect(colour = 'black', linetype = 'solid'),
        legend.title = element_text(size=8),
        plot.margin = unit(c(0,0,10,0), "pt"))

save_plot("single_scales.pdf", single_scales.plt, base_width = 1.4, base_height = 1.4)

datsummary <- datmelt %>% group_by(method, variable) %>% summarize(low = quantile(value, probs=0.025)[[1]],
                                                                   high = quantile(value, probs=0.975)[[1]],
                                                                   m=median(value, na.rm = TRUE))
datsummary$stim <- sapply(datsummary$variable, function(x) stimulus[[toString(x)]])

summary.plt <- ggplot(datsummary, aes(x=stim, y=m, colour=method, shape=method)) +
  stat_function(aes(colour='ground truth', shape='ground truth'), fun=function(x) function.zoo[['sqrt']](x)*10) +
  geom_line(aes(x=stim, y=low), linetype = "dotted") +
  geom_line(aes(x=stim, y=high), linetype = "dotted") +
  geom_ribbon(aes(x=stim, ymin=low, ymax=high), alpha =0.05, linetype = "blank") +
  geom_line() + geom_point() +
  scale_colour_manual(values=solpal5[c(1,4,2)], labels=c('MLDS', 'ground truth', 'BDS')) +
  scale_shape_manual(values = c(15, 32, 17), labels=c('MLDS', 'ground truth', 'BDS')) +
  xlab("stimulus") + ylab("estimated scale") +
  scale_y_continuous(limits = c(0, NA), breaks = c(10)) +
  theme(legend.justification=c(0,1), legend.position=c(0,1),
        axis.text.y = element_blank(),#element_text(size=6),
        axis.title.x = element_blank(),#element_text(size=18),
        axis.title.y = element_blank(),#element_text(size=10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size=7),
        #            legend.background = element_rect(colour = 'black', linetype = 'solid'),
        legend.title = element_text(size=8),
        plot.margin = unit(c(0,0,10,0), "pt"))

save_plot("summary.pdf", summary.plt, base_width = 1.6, base_height = 1.6)
