library(bds)
library(parallel)
source('../simulate.R')

options(mc.cores=parallel::detectCores())
rstan_options(auto_write=TRUE)

# Influence of lapse rate
# State lapse rates to simulate data for
num.sims <- 100
stimulus <- c(0, seq(0.025, 0.975, len=9), 1)
levels <- length(stimulus)
num.trials <- c(2,4,6,8) * choose(length(stimulus), 3)
precision <- c(2.5, 5, 10, 15, 20)

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
  list(psi = stimulus[2:(length(stimulus)-1)],
       precision = 4)
}
init_list <- rep(list(init_fun()), times=4)

run.stan.prec <- function(model_obj, prior_params, tr, pr, function.name, method) {
  df <- data.frame(sc=numeric(), ci.low=numeric(), ci.high=numeric(), pos=factor(), method=factor(), fn=factor(), lps=numeric(), lvl=numeric(), trials=numeric(), prec=numeric())
  for (i in 1:num.sims) {

    lps <- runif(1, min=0.0, max=0.2)
    prior_params$lapses <- lps

    sl <- simulate.responses(intensities=stimulus, trials=tr, simulations=1, precision=pr,
                                  scalefun=function.zoo[[function.name]], lapserate=lps,
                                  sdt=FALSE)
    sim <- sl[['simulations']][[1]]
    time.hmc <- system.time({
      fitobj <- sample_bds_model(model_obj,
                      sim,
                      prior_params,
                      init_list=init_list,
                      .cores=1)

      scale <- extractor_fixed_lapserate(fitobj$stanfit, stimulus, fitobj$data)
      disjoint <- ppc_ordered_residuals(scale)$disjoint
      pval <- ppc_residual_run(scale)$pval
    })

    utime <- time.hmc[1] + time.hmc[4]
    stime <- time.hmc[2] + time.hmc[5]
    rtime <- time.hmc[3]

    sampler_params <- get_sampler_params(fitobj$stanfit, inc_warmup = FALSE)
    divergent <- mean(sapply(sampler_params, function(x) mean(x[, "divergent__"])))

    sc <- c(get_scale_values(scale), get_precision(scale), pval, disjoint, utime, stime, rtime, divergent)
    gt <- c(stimulus, pr, rep(NA, times=6))
    ci.low <- c(get_scale_credible_interval(scale)$ci.low, get_precision_credible_interval(scale)$ci.low, rep(NA, times=6))
    ci.high <- c(get_scale_credible_interval(scale)$ci.high, get_precision_credible_interval(scale)$ci.high, rep(NA, times=6))
    pos <- c(0:(levels - 1), 'sigma', 'p-value', 'disjoint', 'utime', 'stime', 'rtime', 'divergent')

    df <- rbind(df, data.frame(sc=sc,
                               gt=gt,
                               ci.low=ci.low,
                               ci.high=ci.high,
                               pos=pos,
                               method=rep(method, times=levels+7),
                               fn=rep(function.name, times=levels+7),
                               lps=rep(lps, times=levels+7),
                               lvl=rep(11, times=levels+7),
                               trials=rep(tr, times=levels+7),
                               prec=rep(pr, times=levels+7)))
  }

  df
}

run.prec <- function(fun, fn, tr, pr) {
  if (! file.exists(paste('data/prec', fn, tr, pr, 'sim.csv', sep = '-'))) {
    lapses.df <- data.frame(sc=numeric(), gt=numeric(), ci.low=numeric(), ci.high=numeric(),
                            pos=factor(), method=factor(), fn=factor(),
                            lps=numeric(), lvl=numeric(), trials=numeric(), prec=numeric())

    num.lvl <- length(stimulus)
    lapses.df <- rbind(lapses.df, run.stan.prec(raised.cos.model, raised.cos$default_params, tr, pr, fn, 'raised cosine'))
    lapses.df <- rbind(lapses.df, run.stan.prec(half.gauss.model, half.gauss$default_params, tr, pr, fn, 'half-normal'))
    lapses.df <- rbind(lapses.df, run.stan.prec(uniform.model, uniform$default_params, tr, pr, fn, 'uniform'))
    write.table(lapses.df, paste('data/prec', fn, tr, pr, 'sim.csv', sep = '-'), row.names=FALSE, sep='\t')

    # run garbage collection to remove memory-intensive stan fits
    gc()
  }

  TRUE
}

sim_params <- expand.grid(fn=names(function.zoo),
                          tr=num.trials,
                          pr=precision)

nsims <- nrow(sim_params)

mcmapply(run.prec,
         function.zoo[sim_params$fn],
         sim_params$fn,
         sim_params$tr,
         sim_params$pr,
         mc.silent = TRUE,
         mc.preschedule = FALSE)
