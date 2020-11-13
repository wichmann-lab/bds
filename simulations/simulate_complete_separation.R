source("simulate.R")
rstan_options(auto_write = TRUE)

# Influence of lapse rate
# State lapse rates to simulate data for
lapses <- 0
num.sims <- 144
st.list <- list('6'= seq(0, 1, len=6), '11'=c(0, seq(0.025, 0.975, len=9), 1), '15'=seq(0, 1, len=15))
num.trials <- c(1,2,6) * choose(11, 3)
precision <- 0

simulate('sep', lapses=lapses, stimulus=st.list, num.trials=num.trials, sensitivities=precision, num.sims=num.sims, sdt=FALSE, complete.separation=TRUE)