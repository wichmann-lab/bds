library(bds)
library(reshape2)
library(dplyr)
library(parallel)
library(cowplot)

options(mc.cores=parallel::detectCores())

solarized <- c(base03 = '#002b36',
               base02 = '#073642',
               base01 = '#486e75',
               base00 = '#657b83',
               base0  = '#839496',
               base1  = '#93a1a1',
               base2  = '#eee8d5',
               base3  = '#fdf6e3',
               yellow = '#b58900',
               orange = '#cb4b16',
               red    = '#dc322f',
               magenta= '#d33682',
               violet = '#6c71c4',
               blue   = '#268bd2',
               cyan   = '#2aa198',
               green  = '#859900')

solpal5 <- solarized[c('red','blue','yellow','base02','base1')]
solpal5 <- unname(solpal5)


# 

setwd("~/Projects/bds/simulations/priors/plots/divergent_transitions/")
load("~/Projects/bds/simulations/priors/plots/divergent_transitions/demo-square-11-330-10-0.01-virt_exp.RData")

sim.div <- sim.lst$simulations[[20]]

adapt_delta = c(0.8, 0.9, 0.99, 0.999, 0.9999, 0.99999, 0.999999)

ad.div.data <- data.frame()

for (ad in adapt_delta) {
  diff_scale <- bds(sim.div, adapt_delta=ad)
  cc <- convergence.check(diff_scale$stanfit)
  
  df = data.frame(adapt_delta=ad,
                  divergence=mean(cc$chain.divergence),
                  max_rhat=max(cc$diagnostics$rhat),
                  min_ess=min(cc$diagnostics$ess.bulk),
                  min_tail_ess=min(cc$diagnostics$ess.tail))
  ad.div.data <- rbind(ad.div.data, df)
}

ad.div.plot <- ggplot(ad.div.data, aes(x=adapt_delta, y=divergence)) +
  geom_point() + geom_line() + theme_classic()
  
ggsave('adapt_delta_divergence.pdf', ad.div.plot)

ad.rhat.plot <- ggplot(ad.div.data, aes(x=adapt_delta, y=max_rhat)) +
  geom_point() + geom_line() + theme_classic()
  
ggsave('adapt_delta_rhat.pdf', ad.rhat.plot)

ad.ess.plot <- ggplot(ad.div.data, aes(x=adapt_delta, y=min_ess)) +
  geom_point() + geom_line() + theme_classic()
  
ggsave('adapt_delta_ess.pdf', ad.ess.plot)

ad.ess.tail.plot <- ggplot(ad.div.data, aes(x=adapt_delta, y=min_tail_ess)) +
  geom_point() + geom_line() + theme_classic()
  
ggsave('adapt_delta_tail_ess.pdf', ad.ess.tail.plot)
