source("analysis_tools.R")

setwd("~/Projects/bds/simulations/priors/plots/divergent_transitions/")
load("~/Projects/bds/simulations/priors/plots/divergent_transitions/demo-square-11-330-10-0.01-virt_exp.RData")

#no divergent transitions
sim.nodiv <- sim.lst$simulations[[6]]

diff_scale.nodiv <- bds(sim.nodiv)

plot.diagnostics("nodiv", diff_scale.nodiv)

sc.lp.nodiv <- grid.eval("nodiv", diff_scale.nodiv)
summarize.posterior("nodiv", diff_scale.nodiv, sc.lp.nodiv)
#animate.posterior_landscape("nodiv", diff_scale.nodiv, sc.lp.nodiv)

# divergent transitions

sim.div <- sim.lst$simulations[[20]]
#sim <- sim.lst$simulations[[10]] 

#load("~/Projects/bds/simulations/priors/plots/divergent_transitions/demo-square-11-330-10-0.05-virt_exp.RData")

#sim <- sim.lst$simulations[[32]]
#sim <- sim.lst$simulations[[56]]
#sim <- sim.lst$simulations[[73]]
#sim <- sim.lst$simulations[[84]]
#sim <- sim.lst$simulations[[95]]
#sim <- sim.lst$simulations[[28]]

diff_scale.div <- bds(sim.div)

plot.diagnostics("div", diff_scale.div)

sc.lp.div <- grid.eval("div", diff_scale.div)
summarize.posterior("div", diff_scale.div, sc.lp.div)
#animate.posterior_landscape("div", diff_scale.div, sc.lp.div)
