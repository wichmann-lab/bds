library(bds)
source('../simulate.R')

rstan_options(auto_write=TRUE)

# Influence of lapse rate
# State lapse rates to simulate data for
lapses <- c(0, 0.05, 0.1)
num.sims <- 2
stimulus <- list('11'=c(0, seq(0.025, 0.975, len=9), 1))
num.trials <- c(2) * choose(length(stimulus[[1]]), 3)
precision <- c(10)

raised.cos <- build_model(priors=list(psi.uniform, prec.raised_cosine, lapses.const),
                          model=bds.model,
                          extractor_function = extractor_fixed_lapserate)
raised.cos.model <- stan_model(model_code = raised.cos$model_code)

half.gauss <- build_model(priors=list(psi.uniform, prec.halfgauss, lapses.const),
                          model=bds.model,
                          extractor_function = extractor_fixed_lapserate)
half.gauss.model <- stan_model(model_code = half.gauss$model_code)

uniform <- build_model(priors=list(psi.uniform, prec.uniform, lapses.const),
                          model=bds.model,
                          extractor_function = extractor_fixed_lapserate)
uniform.model <- stan_model(model_code = uniform$model_code)

init_fun <- function() {
  list(psi = stimulus[[1]][2:(length(stimulus[[1]])-1)],
       precision = 4)
}
init_list <- rep(list(init_fun()), times=4)

run.stan.prec <- function(model_obj, prior_params, simlist, lps, levels, function.name, method) {
  df <- data.frame(sc=numeric(), ci.low=numeric(), ci.high=numeric(), pos=factor(), method=factor(), fn=factor(), lps=numeric(), lvl=numeric(), trials=numeric(), prec=numeric())
  for (sim in simlist[['simulations']]) {

    prior_params$lapses <- lps

    time.hmc <- system.time({
      fitobj <- sample_bds_model(model_obj,
                      sim,
                      prior_params,
                      init_list=init_list,
                      .cores=1)

      scale <- extractor_fixed_lapserate(fitobj$stanfit, simlist[['scale']], fitobj$data)
      disjoint <- ppc_ordered_residuals(scale)$disjoint
      pval <- ppc_residual_run(scale)$pval
    })

    utime <- time.hmc[1] + time.hmc[4]
    stime <- time.hmc[2] + time.hmc[5]
    rtime <- time.hmc[3]

    sampler_params <- get_sampler_params(fitobj$stanfit, inc_warmup = FALSE)
    divergent <- mean(sapply(sampler_params, function(x) mean(x[, "divergent__"])))

    sc <- c(get_scale_values(scale), get_precision(scale), pval, disjoint, utime, stime, rtime, divergent)
    gt <- c(simlist[['scale']], simlist[['prec']], rep(NA, times=6))
    ci.low <- c(get_scale_credible_interval(scale)$ci.low, get_precision_credible_interval(scale)$ci.low, rep(NA, times=6))
    ci.high <- c(get_scale_credible_interval(scale)$ci.high, get_precision_credible_interval(scale)$ci.high, rep(NA, times=6))
    pos <- c(0:(levels-1), 'sigma', 'p-value', 'disjoint', 'utime', 'stime', 'rtime', 'divergent')

    df <- rbind(df, data.frame(sc=sc,
                               gt=gt,
                               ci.low=ci.low,
                               ci.high=ci.high,
                               pos=pos,
                               method=rep(method, times=levels+7),
                               fn=rep(function.name, times=levels+7),
                               lps=rep(lps, times=levels+7),
                               lvl=rep(simlist[['lvl']], times=levels+7),
                               trials=rep(simlist[['trials']], times=levels+7),
                               prec=rep(simlist[['prec']], times=levels+7)))
  }

  df
}

run.prec <- function(fun, stim, fn, lvl, tr, pr, lapse, sdt=FALSE, num.sims=144) {
  if (! file.exists(paste('data/prec', fn, lvl, tr, pr, lapse, 'sim.csv', sep = '-'))) {
    lapses.df <- data.frame(sc=numeric(), gt=numeric(), ci.low=numeric(), ci.high=numeric(),
                            pos=factor(), method=factor(), fn=factor(),
                            lps=numeric(), lvl=numeric(), trials=numeric(), prec=numeric())

    # generate simulations for current lapserate
    sim.lst <- simulate.responses(intensities=stim, trials=tr, simulations=num.sims, precision=pr,
                                  scalefun=fun, lapserate=lapse,
                                  sdt=FALSE)

    num.lvl <- length(stim)
    lapses.df <- rbind(lapses.df, run.stan.prec(raised.cos.model, raised.cos$default_params, sim.lst, lapse, num.lvl, fn, 'raised cosine'))
    lapses.df <- rbind(lapses.df, run.stan.prec(half.gauss.model, half.gauss$default_params, sim.lst, lapse, num.lvl, fn, 'half-normal'))
    lapses.df <- rbind(lapses.df, run.stan.prec(uniform.model, uniform$default_params, sim.lst, lapse, num.lvl, fn, 'uniform'))
    write.table(lapses.df, paste('data/prec', fn, lvl, tr, pr, lapse, 'sim.csv', sep = '-'), row.names=FALSE, sep='\t')

    # run garbage collection to remove memory-intensive stan fits
    gc()
  }

  TRUE
}

sim_params <- expand.grid(fn=names(function.zoo),
                          lvl=names(stimulus),
                          tr=num.trials,
                          pr=precision,
                          lps=lapses)

nsims <- nrow(sim_params)

mcmapply(run.prec,
         function.zoo[sim_params$fn],
         stimulus[sim_params$lvl],
         sim_params$fn,
         sim_params$lvl,
         sim_params$tr,
         sim_params$pr,
         sim_params$lps,
         MoreArgs = list(num.sims=num.sims),
         mc.silent = TRUE,
         mc.preschedule = FALSE)
