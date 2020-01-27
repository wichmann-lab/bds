source("simulate.R")
rstan_options(auto_write = TRUE)

# Influence of lapse rate
# State lapse rates to simulate data for
lapses <- c(0, 0.05, 0.1, 0.15, 0.2)
num.sims <- 144
st.list <- list('11'=c(0, seq(0.025, 0.975, len=9), 1))
num.trials <- c(40, 80, 120, 165) #* choose(length(st.list[[1]]), 3)
precision <- 10

simulate('smallN', lapses=lapses, stimulus=st.list, num.trials=num.trials, sensitivities=precision, num.sims=num.sims, sdt=FALSE)