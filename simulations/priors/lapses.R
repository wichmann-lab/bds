library(bds)
library(doParallel)

source('../simulate.R')

options(mc.cores=parallel::detectCores())
rstan_options(auto_write=TRUE)

registerDoParallel()

# Influence of lapse rate
# State lapse rates to simulate data for
num.sims <- 144
stimulus <- c(0, seq(0.025, 0.975, len=9), 1)
levels <- length(stimulus)
num.trials <- c(2, 6) * choose(length(stimulus), 3)
precision <- c(2.5, 5, 10, 15, 20)

#lapserates <- c(0.0, 0.05, 0.1, 0.15, 0.2)

lapserates <- c(0.05, 0.1, 0.15, 0.2)

beta <- build_model(priors=list(psi.uniform, prec.raised_cosine, lapses.beta),
                    model=bds.model,
                    extractor_function= default_extractor)
beta.params <- beta$default_params
beta.params$lpsBeta <- 10

beta.model <- stan_model(model_code = beta$model_code)

fixed <- build_model(priors=list(psi.uniform, prec.raised_cosine, lapses.const),
                          model=bds.model,
                          extractor_function = extractor_fixed_lapserate)

fixed.model <- stan_model(model_code = fixed$model_code)
fixed$default_params$lapses <- 0.02

uniform <- build_model(priors=list(psi.uniform, prec.raised_cosine, lapses.uniform),
                          model=bds.model,
                          extractor_function = default_extractor)
uniform.model <- stan_model(model_code = uniform$model_code)

init_fun <- function() {
  list(psi = stimulus[2:(length(stimulus)-1)],
       precision = 4)
}
init_list <- rep(list(init_fun()), times=4)

run.stan.lps <- function(model_obj, prior_params, extractor, sim.lst, tr, pr, lps, function.name, method) {

  df <- foreach(sim = sim.lst, .combine = rbind, .options.multicore=list(preschedule=FALSE, silent=TRUE)) %dopar% {

    time.hmc <- system.time({
      fitobj <- sample_bds_model(model_obj,
                      sim,
                      prior_params,
                      init_list=init_list,
                      .cores=1)

      scale <- extractor(fitobj$stanfit, stimulus, fitobj$data)
      disjoint <- ppc_ordered_residuals(scale)$disjoint
      pval <- ppc_residual_run(scale)$pval
    })

    utime <- time.hmc[1] + time.hmc[4]
    stime <- time.hmc[2] + time.hmc[5]
    rtime <- time.hmc[3]

    sampler_params <- get_sampler_params(fitobj$stanfit, inc_warmup = FALSE)
    divergent <- mean(sapply(sampler_params, function(x) mean(x[, "divergent__"])))

    sc <- c(get_scale_values(scale), get_precision(scale), get_lapserate(scale), pval, disjoint, utime, stime, rtime, divergent)
    gt <- c(stimulus, pr, lps, rep(NA, times=6))
    ci.low <- c(get_scale_credible_interval(scale)$ci.low, get_precision_credible_interval(scale)$ci.low, get_lapserate_credible_interval(scale)$ci.low, rep(NA, times=6))
    ci.high <- c(get_scale_credible_interval(scale)$ci.high, get_precision_credible_interval(scale)$ci.high, get_lapserate_credible_interval(scale)$ci.high, rep(NA, times=6))
    pos <- c(0:(levels - 1), 'sigma', 'lambda', 'p-value', 'disjoint', 'utime', 'stime', 'rtime', 'divergent')

    data.frame(sc=sc,
                               gt=gt,
                               ci.low=ci.low,
                               ci.high=ci.high,
                               pos=pos,
                               method=rep(method, times=levels+8),
                               fn=rep(function.name, times=levels+8),
                               lps=rep(lps, times=levels+8),
                               lvl=rep(11, times=levels+8),
                               trials=rep(tr, times=levels+8),
                               prec=rep(pr, times=levels+8))
  }

  df
}

run.lps <- function(fun, fn, tr, pr, lps) {
  if (! file.exists(paste('data/lps', fn, tr, pr, lps, 'sim.csv', sep = '-'))) {
    sl <- simulate.responses(intensities=stimulus,
                             trials=tr,
                             simulations=num.sims,
                             precision=pr,
                             scalefun=fun,
                             lapserate=lps,
                             sdt=FALSE)
    sim.lst <- sl[['simulations']]
    
    lapses.df <- run.stan.lps(beta.model, beta$default_params, beta$extractor, sim.lst, tr, pr, lps, fn, 'beta15')
    print('done beta15')
    lapses.df <- rbind(lapses.df, run.stan.lps(beta.model, beta.params, beta$extractor, sim.lst, tr, pr, lps, fn, 'beta110'))
    print('done beta110')
    lapses.df <- rbind(lapses.df, run.stan.lps(fixed.model, fixed$default_params, fixed$extractor, sim.lst, tr, pr, lps, fn, 'fixed'))
    print('done fixed')
    lapses.df <- rbind(lapses.df, run.stan.lps(uniform.model, uniform$default_params, uniform$extractor, sim.lst, tr, pr, lps, fn, 'uniform'))
    print('done uniform')
    write.table(lapses.df, paste('data/lps', fn, tr, pr, lps, 'sim.csv', sep = '-'), row.names=FALSE, sep='\t')

    # run garbage collection to remove memory-intensive stan fits
    gc()
  }

  TRUE
}

sim_params <- expand.grid(fn=names(function.zoo),
                          tr=num.trials,
                          pr=precision,
                          lps=lapserates)

nsims <- nrow(sim_params)

for (i in 1:nsims) {
  run.lps(
         function.zoo[[sim_params$fn[[i]]]],
         sim_params$fn[[i]],
         sim_params$tr[[i]],
         sim_params$pr[[i]],
         sim_params$lps[[i]])
}