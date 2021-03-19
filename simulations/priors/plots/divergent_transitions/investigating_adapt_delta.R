source("analysis_tools.R")

# 

setwd("~/Projects/bds/simulations/priors/plots/divergent_transitions/")
load("~/Projects/bds/simulations/priors/plots/divergent_transitions/demo-square-11-330-10-0.01-virt_exp.RData")

adapt_delta = c(0.8, 0.9, 0.99, 0.999, 0.9999, 0.99999, 0.999999, 0.9999999)
sims = 1:length(sim.lst$simulations)
ad.div.data <- data.frame()

for (s in sims) {
  sim.div <- sim.lst$simulations[[s]]

  diff_scale <- bds(sim.div)
  
#  sc.lp.div <- load.grid(paste0('square-trials=330-sens=10-lps=0.01-sim=', s), diff_scale)
  
  for (ad in adapt_delta) {


    diff_scale <- bds(sim.div, adapt_delta=ad)
    cc <- convergence.check(diff_scale$stanfit)
  
    df = data.frame(adapt_delta=ad,
                    divergence=mean(cc$chain.divergence),
                    max_rhat=max(cc$diagnostics$rhat),
                    min_ess=min(cc$diagnostics$ess.bulk),
                    min_tail_ess=min(cc$diagnostics$ess.tail),
                    sim=s)
    ad.div.data <- rbind(ad.div.data, df)
    
#    summarize.posterior(paste0('square-trials=330-sens=10-lps=0.01-sim=', s), diff_scale, sc.lp.div, postfix=paste0('-adapt_delta=', ad))
  }
}

ad.div.plot <- ggplot(ad.div.data, aes(x=factor(adapt_delta), y=divergence/10000)) +
  geom_violin() +
  ylim(0, NA) +
#  scale_x_continuous(trans='log') +
  theme_classic()
  
ggsave('adapt_delta_divergence.pdf', ad.div.plot)

ad.rhat.plot <- ggplot(ad.div.data, aes(x=adapt_delta, y=max_rhat)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(trans='log') +
  theme_classic()
  
ggsave('adapt_delta_rhat.pdf', ad.rhat.plot)

ad.ess.plot <- ggplot(ad.div.data, aes(x=adapt_delta, y=min_ess)) +
  geom_point() +
  geom_line() +
  ylim(0, NA) +
  scale_x_log10() +
  theme_classic()
  
ggsave('adapt_delta_ess.pdf', ad.ess.plot)

ad.ess.tail.plot <- ggplot(ad.div.data, aes(x=adapt_delta, y=min_tail_ess)) +
  geom_point() +
  geom_line() +
  ylim(0, NA) +
  scale_x_log10() +
  theme_classic()
  
ggsave('adapt_delta_tail_ess.pdf', ad.ess.tail.plot)
