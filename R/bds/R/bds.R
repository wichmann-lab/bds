library(rstan)

order_data <- function(mlds_data) {
  if (ncol(mlds_data) == 4) {
    new.df <- data.frame(Response=numeric(), S1=numeric(), S2=numeric(), S3=numeric())
    for (n in 1:nrow(mlds_data)) {
      if (mlds_data[n,2] > mlds_data[n,4]) {
        new.df <- rbind(new.df, data.frame( Response=1-mlds_data[n,1], S1=mlds_data[n,4], S2=mlds_data[n,3], S3=mlds_data[n,2]))
      } else {
        new.df <- rbind(new.df, data.frame(Response=mlds_data[n,1], S1=mlds_data[n,2], S2=mlds_data[n,3], S3=mlds_data[n,4]))
      }
    }
  } else {
    new.df <- data.frame(Response=numeric(), S1=numeric(), S2=numeric(), S3=numeric(), S4=numeric())
    for (n in 1:nrow(mlds_data)) {
      if (mlds_data[n,2] > mlds_data[n,3]) {
        s2 <- mlds_data[n,2]
        s1 <- mlds_data[n,3]
      } else {
        s1 <- mlds_data[n,2]
        s2 <- mlds_data[n,3]
      }
      if (mlds_data[n,4] > mlds_data[n,5]) {
        s4 <- mlds_data[n,4]
        s3 <- mlds_data[n,5]
      } else {
        s3 <- mlds_data[n,4]
        s4 <- mlds_data[n,5]
      }
      if (mlds_data[n,2] > mlds_data[n,4]) {
        new.df <- rbind(new.df, data.frame(S1=s3, S2=s4, S3=s1, S4=s2, Response=1-mlds_data[n,5]))
      } else {
        new.df <- rbind(new.df, data.frame(S1=s1, S2=s2, S3=s3, S4=s4, Response=mlds_data[n,5]))
      }
    }
  }
  new.df
}

bds <- function(mlds_data,
                stimulus=NULL,
                fit.lapses=TRUE,
                .cores=getOption('mc.cores', default = 1L)) {

  if (is.null(stimulus)) {
    if (ncol(mlds_data) == 4) {
      stimulus <- seq(0,1, len=max(mlds_data[,2:4]))
    } else if (ncol(mlds_data == 5)) {
      stimulus <- seq(0,1, len=max(mlds_data[,2:5]))
    } else {
      stop("Difference scaling data should have 4 or 5 columns!")
    }
  }
  if (fit.lapses) {
    md <- build_model(priors=list(psi.dirichlet, sens.raised_cosine, lapses.beta),
                model=bds.model,
                extractor_function = default_extractor)
    init_fun <- function() {
      list(psi_diff = diff( (stimulus- min(stimulus))/max(stimulus- min(stimulus)) ),
           sensitivity = md$default_params$sensLow,
           lapses = 0.01)
    }
  } else {
    md <- build_model(priors=list(psi.dirichlet, sens.raised_cosine, lapses.const),
                      model=bds.model,
                      extractor_function = extractor_fixed_lapserate)
    init_fun <- function() {
      list(psi_diff = diff( (stimulus- min(stimulus))/max(stimulus- min(stimulus)) ),
           sensitivity = (md$default_params$sensLow + md$default_params$sensHigh)/2.0)
    }
  }

  model_obj <- stan_model(model_code=md$model_code)

  stanfit <- sample_bds_model(model_obj, mlds_data, prior_params=md$default_params, init_list = rep(list(init_fun()), times=4), .cores=.cores)

  md$extractor(stanfit$stanfit, stimulus, stanfit$data)
}

sample_bds_model <- function(model_obj,
                             mlds_data,
                             prior_params,
                             init_list,
                             .cores=getOption('mc.cores', default = 1L)) {
  mlds_data.ordered <- order_data(mlds_data)

  data = c(list(N = length(mlds_data[,1])),
           prior_params)

  if (ncol(mlds_data) == 4) {
    data$S1 <- mlds_data.ordered[,2]
    data$S2 <- mlds_data.ordered[,3]
    data$S3 <- mlds_data.ordered[,3]
    data$S4 <- mlds_data.ordered[,4]
    data$Responses <- mlds_data.ordered[,1]
    data$K <- max(mlds_data.ordered[,2:4])
  } else if(ncol(mlds_data) == 5) {
    data$S1 <- mlds_data.ordered[,2]
    data$S2 <- mlds_data.ordered[,3]
    data$S3 <- mlds_data.ordered[,4]
    data$S4 <- mlds_data.ordered[,5]
    data$Responses <- mlds_data.ordered[,1]
    data$K <- max(mlds_data.ordered[,2:5])
  } else {
    stop("Difference scaling data should have 4 or 5 columns!")
  }

  fit <- sampling(model_obj,
                  data=data,
                  warmup = 1000,
                  iter = 3500, chains=length(init_list),
#                  include = FALSE, pars = c("psi_ext", "decision"),
                  control = list(adapt_delta = 0.999),
                  init=init_list,
                  cores=.cores)

  list(stanfit=fit, data=data)
}

convergence.check <- function(stanfit) {
  arr <- extract(stanfit, permuted=FALSE, pars=c('psi', 'sensitivity', 'lapses'))
  M <- ncol(arr)
  N <- nrow(arr)

  df <- data.frame(par=factor(), rhat=numeric(), ess.bulk=numeric(), ess.tail=numeric())

  warnings <- c(rhat=FALSE, ess.bulk=FALSE, ess.tail=FALSE, div=FALSE)

  for (par in dimnames(arr)$parameters) {
    rhat <- Rhat(arr[,,par])
    ess.b <- ess_bulk(arr[,,par])
    ess.t <- ess_tail(arr[,,par])

    if (rhat > 1.05) {
      warnings['rhat'] <- TRUE
    }

    if (ess.b < M*100) {
      warnings['ess.bulk'] <- TRUE
    }

    if (ess.t < M*100) {
      warnings['ess.tail'] <- TRUE
    }

    df <- rbind(df, data.frame(par=par, rhat=rhat, ess.bulk=ess.b, ess.tail=ess.t))
  }

  div <- get_divergent_iterations(stanfit)
  chains.div <- table((which(div)-1) %/% N + 1)

  if (length(chains.div) > 0) {
    warnings['div'] = TRUE
  }

  list(diagnostics=df, chain.divergence=chains.div, warnings=warnings)
}

grid.eval <- function(diff_scale) {
  pmean <- c(lapses=diff_scale$lapserate,sensitivity=diff_scale$sensitivity,diff(diff_scale$scale))

  log_posterior <- function(x) {
    pars <- list(lapses=x[1], sensitivity=x[2], psi_diff=x[3:length(x)])
    lp <- tryCatch(-log_prob(diff_scale$stanfit, unconstrain_pars(diff_scale$stanfit, pars)),error=function(cond) return(Inf))

    return(lp)
  }

  map.estim <- optim(pmean, log_posterior)$par

#  sc.len <-
#  sc.prior <- ddirichlet(rep(0, sc.len),rep(1, sc.len), log=TRUE)
#  lps.zero.prior <- dbeta(0.0, 1, 5, log.p=TRUE)

  for (sc in scales) {

    delta.regr <- design.matrix %*% sc

    for (sens in sensitivities) {

      decision.prob <- pnorm(sens * delta.regr)
      loglik <- sum(dbinom(resp, 1, decision.prob, log = TRUE))

      density <- loglik + sens.prior + sc.prior + lps.zero.prior
      for (lps in lapses) {
        loglik.mix <- logSumExp(log(1-lps)+loglik, log(lps) + 0.5)

        density <- loglik.mix + sens.prior + lps.prior + sc.prior
      }
    }
  }
}
