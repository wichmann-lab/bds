library(MLDS)
library(bds)
library(psyphy)

#' Simulate MLDS responses
#'
#'\code{sum} Simulates responses of a MLDS experiment with a given lapserate.
#'
#'@param intensities  number of input intensities
#'@param              trials number of simulated trials in one experiment
#'@param simulations  number of experiments to simulate
#'@param precision    height of response function / inverse of the standard
#'                    deviation of the Gaussian noise
#'@param lapserate    proportion of trials to be drawn from a Bernoulli
#'                    distribution with p=0.5, indepent from the response function
#'@param scalefun     function to represent observer response
#'@param sdt          Should the noise be applied before or after
#'
#'@return Returns either a dataframe or a list of dataframes.
#'        Dataframe has fields resp, S1, S2, S3.
simulate.responses <- function(intensities,
                               trials,
                               simulations=100,
                               precision=10,
                               lapserate=0.0,
                               scalefun=function(x) x^(1/2),
                               sdt=TRUE) {

  # Generate triads corresponding to the number of trials
  num_triads <- choose(intensities, 3)
  repl <- ceiling(trials / num_triads)
  Tr <- do.call(rbind, replicate(repl, t(combn(intensities,3)), simplify=FALSE))

  # If the number of unique triads doesn't divide trials, downsample triads
  if (nrow(Tr) > trials) {
    Tr <- Tr[sample(1:nrow(Tr), trials, replace=FALSE),]
  }

  Sc <- scalefun(seq(0, 1, len=intensities))
  sim.lst <- list()
  sim.lst[['simulations']] <- list()

  # Replace responses with randomly sampled coin flips according to lapse rate
  for (n in 1:simulations) {
    if (sdt) {
      det_noise <- 1/(2*precision)
      decision <- mapply(function(a, b, c) abs(rnorm(1, Sc[c], det_noise) - rnorm(1, Sc[b], det_noise)) -
                                              abs(rnorm(1, Sc[b], det_noise) - rnorm(1, Sc[a], det_noise)),
                         Tr[,1], Tr[,2], Tr[,3])
    } else {
      decision <- mapply(function(a, b, c) rnorm(1, Sc[a] - 2*Sc[b] + Sc[c], 1/precision),
                         Tr[,1], Tr[,2], Tr[,3])
    }

    responses <- decision > 0

    if (lapserate > 0) {
      # draw positions to be replaced
      lapse.pos <- rbinom(trials, 1, lapserate)

      # draw results of possible lapses
      possible.lapses <- rbinom(trials, 1, 0.5)

      # replace drawn positions with drawn lapses
      responses <- (1-lapse.pos) * responses + lapse.pos * possible.lapses
    }
    sim.lst[['simulations']][[n]] <- data.frame(resp=responses,S1=Tr[,1], S2=Tr[,2], S3=Tr[,3])
  }
  sim.lst[["scale"]] = Sc
  sim.lst[["prec"]] = precision
  # return list of simulation results
  sim.lst
}

asym.bootstrap <- function(asym.fit) {
  lvl <- ncol(asym.fit$data)
  tr <- nrow(asym.fit$data)
  sim.coef <- coef(asym.fit)
  n <- length(sim.coef)
  prec <- sim.coef[lvl-1]
  scale <- c(0, sim.coef[1:(lvl-2)]/prec, 1)

  fn <- function(x) ifelse(x<1.0, scale[floor(x*lvl)+1], scale[lvl])
  sim.lst <- simulate.responses(lvl, tr, simulations=1000, precision=prec, lapserate=0.0, scalefun=fn)

  bt.fits <- lapply(sim.lst, function(sim) tryCatch(psyfun.2asym(asym.fit$formula,
                                                        data = as.data.frame(as.dm(sim)), link = probit.2asym), error=function(e) NA))
  bt.filtered <- bt.fits[!is.na(bt.fits)]
  bt.coef <- lapply(bt.filtered, function(fit) c(fit$coefficients, fit$lambda, fit$gam, use.names=FALSE))
  print(bt.coef)
  if (!is.list(bt.coef) | length(bt.coef) == 0) {
    obs.low <- rep(NA, times=n+2)
    obs.high <- rep(NA, times=n+2)
  } else {
    mt.coef <- matrix(unlist(bt.coef, use.names = FALSE), ncol=n+2, byrow=TRUE)

    samples <- apply(mt.coef, 1, function(x) c(x[1:(n-1)]/x[n], x[n:(n+2)]))

    obs.low <- apply(samples, 1, quantile, probs = 0.025, na.rm=TRUE)
    obs.high <- apply(samples, 1, quantile, probs = 0.975, na.rm=TRUE)
  }
  list(low=obs.low, high=obs.high)
}

run.mlds <- function(simlist, lps, levels, function.name, fac=-1) {
  df <- data.frame(sc=numeric(), ci.low=numeric(), ci.high=numeric(), pos=factor(), method=factor(), fn=factor(), lps=numeric())
  for (sim in simlist[['simulations']]) {

    time.mlds <- system.time({
      fit <- mlds(sim, seq(0,1,len=levels))
      obs.bt <- boot.mlds(fit, 10000)
      obs.diag <- binom.diagnostics(fit, nsim=10000)

      pval <- obs.diag$p
      resid <- apply(obs.diag$resid, 2, quantile, probs=c(0.025, 0.975))

      overlap <- mean(resid[1,] < sort(obs.diag$Obs.resid) & resid[2,] > sort(obs.diag$Obs.resid))

      n <- nrow(obs.bt$boot.samp)
      samples <- apply(obs.bt$boot.samp, 2, function(x) c(x[1:(n-2)], x[n]))
      obs.low <- apply(samples, 1, quantile, probs = 0.025)
      obs.high <- apply(samples, 1, quantile, probs = 0.975)
    })

    time.asym <- system.time({
      asym <- tryCatch(psyfun.2asym(cbind(resp, 1 - resp) ~ . - 1,
                       data = fit$obj$data, link = probit.2asym),
                       error=function(e) NA)

      if (!is.na(asym)) {
        as.co <- c(coef(asym), asym$lambda, asym$gam)

        asym.bt <- asym.bootstrap(asym)
      } else {
        as.co <- rep(NA, levels+1)
        asym.bt <- list(low=rep(NA, times=levels+1), high=rep(NA, times=levels+1))
      }
    })

    utime.mlds <- time.mlds[1] + time.mlds[4]
    stime.mlds <- time.mlds[2] + time.mlds[5]
    rtime.mlds <- time.mlds[3]

    utime.asym <- time.asym[1] + time.asym[4]
    stime.asym <- time.asym[2] + time.asym[5]
    rtime.asym <- time.asym[3]

    sc <- c(0, 0, 1, 1, fit$pscale[2:(levels-1)]/fit$pscale[levels], fit$pscale[levels], pval, overlap, utime.mlds, stime.mlds, rtime.mlds,
            as.co[1:(levels-2)]/as.co[levels-1], as.co[(levels-1):(levels+1)], utime.asym, stime.asym, rtime.asym)
    gt <- c(0, 0, 1, 1, simlist[['scale']][2:(levels-1)], simlist[['prec']], rep(NA, times=5),
            simlist[['scale']][2:(levels-1)], simlist[['prec']], lps/2, lps/2, rep(NA, times=3))
    ci.low <- c(0, 0, 1, 1, obs.low, rep(NA, times=5), asym.bt$low, rep(NA, times=3))
    ci.high <- c(0, 0, 1, 1, obs.high, rep(NA, times=5), asym.bt$high, rep(NA, times=3))
    pos <- c(0, 0, levels-1, levels-1, 1:(levels-2), 'sigma', 'p-value', 'overlap', 'utime', 'stime', 'rtime', 1:(levels-2), 'sigma', 'lambda', 'gamma', 'utime', 'stime', 'rtime')
    method <- c('glm', 'asym', 'glm', 'asym', rep('glm', times=levels+4), rep('asym', times=levels+4))

    df <- rbind(df, data.frame(sc=sc,
                               gt=gt,
                               ci.low=ci.low,
                               ci.high=ci.high,
                               pos=pos,
                               method=method,
                               fn=rep(function.name, times=2*levels+12),
                               lps=rep(lps, times=2*levels+12), fac=rep(fac, times=2*levels+12)))
  }

  df
}

#' Run simple stan model on simulated MLDS experiment
run.stan <- function(simlist, lps, levels, function.name, fac=-1) {
  df <- data.frame(sc=numeric(), ci.low=numeric(), ci.high=numeric(), pos=factor(), method=factor(), fn=factor(), lps=numeric())
  for (sim in simlist[['simulations']]) {
    time.hmc <- system.time({
      fit <- bds(sim, fit.lapses = FALSE)
      disjoint <- ppc_ordered_residuals(fit)$disjoint
      pval <- ppc_residual_run(fit)$pval
    })

    utime <- time.hmc[1] + time.hmc[4]
    stime <- time.hmc[2] + time.hmc[5]
    rtime <- time.hmc[3]

    sampler_params <- get_sampler_params(fit$stanfit, inc_warmup = FALSE)
    divergent <- mean(sapply(sampler_params, function(x) mean(x[, "divergent__"])))

    sc <- c(get_scale_values(fit), get_precision(fit), pval, disjoint, utime, stime, rtime, divergent)
    gt <- c(simlist[['scale']], simlist[['prec']], NA, NA, NA, NA, NA, NA)
    ci.low <- c(get_scale_credible_interval(fit)$ci.low, get_precision_credible_interval(fit)$ci.low, NA, NA, NA, NA, NA, NA)
    ci.high <- c(get_scale_credible_interval(fit)$ci.high, get_precision_credible_interval(fit)$ci.high, NA, NA, NA, NA, NA, NA)
    pos <- c(0:(levels-1), 'sigma', 'p-value', 'disjoint', 'utime', 'stime', 'rtime', 'divergent')

    df <- rbind(df, data.frame(sc=sc,
                               gt=gt,
                               ci.low=ci.low,
                               ci.high=ci.high,
                               pos=pos,
                               method=rep('stan', times=levels+7),
                               fn=rep(function.name, times=levels+7),
                               lps=rep(lps, times=levels+7),
                               fac=rep(fac, times=levels+7)))
  }

  df
}

#' Run mixture model on simulated MLDS experiment
run.stan.lapse <- function(simlist, lps, levels, function.name, fac=-1) {
  df <- data.frame(sc=numeric(), ci.low=numeric(), ci.high=numeric(), pos=factor(), method=factor(), fn=factor(), lps=numeric())

  for (sim in simlist[['simulations']]) {
    time.hmc <- system.time({
      fit <- bds(sim, fit.lapses = TRUE)
      overlap <- ppc_ordered_residuals(fit)$pval
      pval <- ppc_residual_run(fit)$pval
    })

    utime <- time.hmc[1] + time.hmc[4]
    stime <- time.hmc[2] + time.hmc[5]
    rtime <- time.hmc[3]

    sampler_params <- get_sampler_params(fit$stanfit, inc_warmup = FALSE)
    divergent <- mean(sapply(sampler_params, function(x) mean(x[, "divergent__"])))

    sc <- c(get_scale_values(fit), get_precision(fit), get_lapserate(fit), pval, overlap, utime, stime, rtime, divergent)
    gt <- c(0, simlist[['scale']][2:levels-1], 1, simlist[['prec']], lps, NA, NA, NA, NA, NA, NA)
    ci.low <- c(get_scale_credible_interval(fit)$ci.low, get_precision_credible_interval(fit)$ci.low, get_lapserate_credible_interval(fit)$ci.low, NA, NA, NA, NA, NA, NA)
    ci.high <- c(get_scale_credible_interval(fit)$ci.high, get_precision_credible_interval(fit)$ci.high, get_lapserate_credible_interval(fit)$ci.high, NA, NA, NA, NA, NA, NA)
    pos <- c(0:(levels-1), 'sigma', 'lambda', 'p-value', 'overlap', 'utime', 'stime', 'rtime', 'divergent')

    df <- rbind(df, data.frame(sc=sc,
                               gt=gt,
                               ci.low=ci.low,
                               ci.high=ci.high,
                               pos=pos,
                               method=rep('mixture', times=levels+8),
                               fn=rep(function.name, times=levels+8),
                               lps=rep(lps, times=levels+8),
                               fac=rep(fac, times=levels+8)))
  }

  df
}

munsell <- function(x, ref) ifelse( (x/ref) <= (6. / 29) ^ 3,
                                     11.6 * 841. / 108 * x/ref + 4. / 29 - 1.6,
                                     11.6* (x/ref)^(1./3)                - 1.6)

function.zoo <- list('square'    = function(x) x^2,
                     'id'        = function(x) x,
                     'sqrt'      = function(x) x^0.5,
                     'saddle'    = function(x) 4*(x-0.5)^3 + 0.5,
                     'invsaddle' = function(x) sign(x/4 - 0.125) * (abs(x/4 - 0.125))^(1/3) + 0.5,
                     'log'       = function(x) (log(x+0.05) - log(0.05))/(log(1.05) - log(0.05)),
                     'exp'       = function(x) 0.05 * ( (1.05/0.05)^x-1),
                     'logit'     = function(x) (tanh(5*(x-0.5)) - tanh(-0.5*5))/(tanh(5*0.5) - tanh(-0.5*5)),
                     'invlogit'  = function(x) atanh((tanh(5*0.5) - tanh(-0.5*5))*x + tanh(-0.5*5))/5 + 0.5,
                     'munsell'   = function(x) munsell(x, 1.0)/10
                    )

run.simulations <- function(factor.name,
                            functions=function.zoo,
                            lapses=seq(0,0.05,len=6),
                            levels=list(10),
                            num.trials=list(1000),
                            precisions=list(10),
                            num.sims=100,
                            sdt=FALSE) {


  if (! file.exists(paste0('data/', factor.name, '.Rdata'))) {
    sim_params <- expand.grid(fn=names(functions),
                              lvl=levels,
                              tr=num.trials,
                              pr=precisions,
                              lapse=lapses)
    next.row <- 1
    fcon <- file(paste0('data/sim-lapses-x-', factor.name, '.csv'), open='w')
    write.table(data.frame(sc=numeric(), ci.low=numeric(), ci.high=numeric(),
                           pos=factor(), method=factor(), fn=factor(),
                           lps=numeric(), fac=factor()),
                fcon, sep='\t')
  } else {
    fcon <- file(paste0('data/sim-lapses-x-', factor.name, '.csv'), open='a')
    load(paste0('data/', factor.name, '.Rdata'))
  }

  nsims <- nrow(sim_params)
  con <- txtProgressBar(min=1,max=nsims,style=3)

  setTxtProgressBar(con, next.row)

  # Run simulations for every lapse rate.
  # We only look at the last scale value, which should recover the precision
  # used in the simulation.
  for (row in next.row:nsims) {
    fn <- sim_params[row, 'fn']
    lvl <- sim_params[row, 'lvl']
    tr <- sim_params[row, 'tr']
    pr <- sim_params[row, 'pr']
    lapse <- sim_params[row, 'lapse']

    lapses.df <- data.frame(sc=numeric(), gt=numeric(), ci.low=numeric(), ci.high=numeric(),
                            pos=factor(), method=factor(), fn=factor(),
                            lps=numeric(), fac=factor())

    # generate simulations for current lapserate
    sim.lst <- simulate.responses(lvl, tr, simulations=num.sims, precision=pr,
                                  scalefun=functions[[fn]], lapserate=lapse,
                                  sdt=sdt)

    fac <- switch(factor.name,
                  level=lvl,
                  abs=lvl,
                  trials=tr,
                  precision=pr,
                  -1)

    invisible(capture.output(lapses.df <- rbind(lapses.df, run.mlds(sim.lst, lapse, lvl, fn, fac=fac))))
    invisible(capture.output(lapses.df <- rbind(lapses.df, run.stan(sim.lst, lapse, lvl, fn, fac=fac))))
    invisible(capture.output(lapses.df <- rbind(lapses.df, run.stan.lapse(sim.lst, lapse, lvl, fn, fac=fac))))

    write.table(lapses.df, fcon, row.names=FALSE, col.names = FALSE, sep='\t')
    flush(fcon)

    next.row <- next.row + 1
    save(next.row, sim_params, file=paste0('data/', factor.name, '.Rdata'))

    setTxtProgressBar(con, row)

    # run garbage collection to remove memory-intensive stan fits
    gc()
  }

  close(fcon)
  close(con)
}
