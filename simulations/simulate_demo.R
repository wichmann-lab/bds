source("simulate.R")
rstan_options(auto_write = TRUE)

# Influence of lapse rate
# State lapse rates to simulate data for
lapses <- c(0, 0.03, 0.05)
num.sims <- 2
st.list <- list('11'=c(0, seq(0.025, 0.975, len=9), 1))
num.trials <- 8 * choose(length(st.list[[1]]), 3)
precision <- 10

# simulate('demo', lapses=lapses, stimulus=st.list, num.trials=num.trials, sensitivities=precision, num.sims=num.sims, sdt=FALSE)
run.simulations(function.zoo[['id']], st.list[['11']], 'id', '11', num.trials, 10.0, 0.05, num.sims=2)