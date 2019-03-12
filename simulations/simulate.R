library(MLDS)
library(bds)
library(reticulate)
library(psyphy)

source('bds.R')

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
    sim.lst[[n]] <- data.frame(resp=responses,S1=Tr[,1], S2=Tr[,2], S3=Tr[,3])
  }
  sim.list[["scale"]] = Sc
  sim.list[["prec"]] = precision
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
  for (sim in simlist) {
    fit <- mlds(sim, seq(0,1,len=levels))

    obs.bt <- boot.mlds(fit, 1000)

    n <- nrow(obs.bt$boot.samp)
    samples <- apply(obs.bt$boot.samp, 2, function(x) c(x[1:(n-2)], x[n]))
    obs.low <- apply(samples, 1, quantile, probs = 0.025)
    obs.high <- apply(samples, 1, quantile, probs = 0.975)

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

    df <- rbind(df, data.frame(sc=c(0, 0, 1, 1, fit$pscale[2:(levels-1)]/fit$pscale[levels], fit$pscale[levels],
                                    0, 0, 1, 1, as.co[1:(levels-2)]/as.co[levels-1], as.co[(levels-1):(levels+1)]),
                               gt=c(0, 0, 1, 1, simlist[['scale']][2:(levels-1)], prec,
                                    simlist[['scale']][2:(levels-1)], simlist['prec'], lps/2, lps/2),
                               ci.low=c(0, 0, 1, 1, obs.low, asym.bt$low),
                               ci.high=c(0, 0, 1, 1, obs.high, asym.bt$high),
                               pos=c(0, 0, levels-1, levels-1, 1:(levels-2), 'sigma', 1:(levels-2), 'sigma', 'lambda', 'gamma'),
                               method=c('glm', 'asym', 'glm', 'asym', rep('glm', times=levels-1), rep('asym', times=levels+1)),
                               fn=rep(function.name, times=2*levels+4),
                               lps=rep(lps, times=2*levels+4), fac=rep(fac, times=2*levels+4)))
  }

  df
}

#' Run simple stan model on simulated MLDS experiment
run.stan <- function(simlist, lps, levels, function.name, fac=-1) {
  df <- data.frame(sc=numeric(), ci.low=numeric(), ci.high=numeric(), pos=factor(), method=factor(), fn=factor(), lps=numeric())
  for (sim in simlist) {
    fit <- bds.stan(sim, "bds")

    df <- rbind(df, data.frame(sc=c(0, fit$point_est[1:(levels-1)], 1),
                               gt=c(simlist[['scale']], simlist[['prec']]),
                               ci.low=c(0, fit$hdi.low[1:(levels-1)], 1),
                               ci.high=c(0, fit$hdi.high[1:(levels-1)], 1),
                               pos=c(0:levels-1, 'sigma'),
                               method=rep('stan', times=levels+1),
                               fn=rep(function.name, times=levels+1),
                               lps=rep(lps, times=levels+1), fac=rep(fac, times=levels+1)))
  }

  df
}

#' Run mixture model on simulated MLDS experiment
run.stan.lapse <- function(simlist, lps, levels, function.name, fac=-1) {
  df <- data.frame(sc=numeric(), ci.low=numeric(), ci.high=numeric(), pos=factor(), method=factor(), fn=factor(), lps=numeric())

  for (sim in simlist) {
    fit <- bds.stan(sim, "lapse_bds")

    df <- rbind(df, data.frame(sc=c(0, 1, fit$point_est[1:levels]),
                               gt=c(0, 1, simlist[['scale']][2:levels-1], simlist[['prec']]),
                               ci.low=c(0, 1, fit$hdi.low[1:levels]),
                               ci.high=c(0, 1, fit$hdi.high[1:levels]),
                               pos=c(0, levels-1, 1:(levels-2), 'sigma', 'lambda'),
                               method=rep('mixture', times=levels+2),
                               fn=rep(function.name, times=levels+2),
                               lps=rep(lps, times=levels+2), fac=rep(fac, times=levels+2)))
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
