source("simulate.R")

# Influence of lapse rate
# State lapse rates to simulate data for
lapses <- c(0, 0.03, 0.05)
num.sims <- 2
num.levels <- 10
num.trials <- 8 * choose(num.levels, 3)
precision <- 10

run.simulations('demo', lapses=lapses, levels=num.levels, num.trials=num.trials, precisions=precision, num.sims=num.sims, sdt=FALSE)
