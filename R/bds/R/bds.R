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
                .cores=getOption('mc.cores', default = 1L),
		adapt_delta=0.9) {

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

  stanfit <- sample_bds_model(model_obj, mlds_data, prior_params=md$default_params, init_list = rep(list(init_fun()), times=4), .cores=.cores, adapt_delta=adapt_delta)

  md$extractor(stanfit$stanfit, stimulus, stanfit$data)
}

sample_bds_model <- function(model_obj,
                             mlds_data,
                             prior_params,
                             init_list,
                             .cores=getOption('mc.cores', default = 1L),
			     adapt_delta=0.9) {
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
                  control = list(adapt_delta = adapt_delta),
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

create.design_matrix <- function(diff_scale) {
  k <- diff_scale$K

  expand_row <- function(s1, s2, s3, s4) {
    r <- rep(0, k)

    r[s1] <- 1
    r[s2] <- r[s2] - 1
    r[s3] <- r[s3] - 1
    r[s4] <- 1

    return(r[2:k])
  }

  return(t(mapply(expand_row, diff_scale$S1, diff_scale$S2, diff_scale$S3, diff_scale$S4, SIMPLIFY = TRUE)))
}

grid.eval <- function(diff_scale, sensitivities=seq(5,20,len=11), lapses=seq(0,0.2, len=11) {
  pmean <- c(lapses=diff_scale$lapserate,sensitivity=diff_scale$sensitivity,diff(diff_scale$scale))

  log_posterior <- function(x) {
    pars <- list(lapses=x[1], sensitivity=x[2], psi_diff=x[3:length(x)])
    lp <- tryCatch(-log_prob(diff_scale$stanfit, unconstrain_pars(diff_scale$stanfit, pars)),error=function(cond) return(Inf))

    return(lp)
  }

  draised.single <- function(y, start, u1, u2, end) {
    if (start < y && y < u1) {
      res = log(0.5-0.5*cos(pi/(u1-start)*(y-start)));
    } else if (y > u2 && y < end) {
      res = log(0.5+0.5*cos(pi/(end-u2)*(y-u2)));
    } else if (u1 <= y && y <= u2) {
      res = 0.0;
    } else {
      res = -Inf
    }

    res
  }

  draised <- function(x, a, b, c, d) {
    return(sapply(x, draised.single, start=a, u1=b, u2=c, end=d))
  }

  map.estim <- optim(pmean, log_posterior)$par
  
  design.matrix <- create.design_matrix(diff_scale)
  sc.len <- diff_scale$K - 1
  sc.prior <- ddirichlet(rep(0, sc.len-1),rep(1, sc.len-1), log=TRUE)
#  lps.zero.prior <- dbeta(0.0, 1, 5, log.p=TRUE)

  sens.prior <- draised(sensitivities, 2.5, 5, 25, 50)
  lps.prior <- dbeta(lapses, 1, 5, log.p=TRUE)

  densities <- expand.grid(scales, sensitivities, lapses)
  
  for (sc in scales) {

    delta.regr <- design.matrix %*% sc

    for (si in 1:len(sensitivities)) {

      decision.prob <- pnorm(sensitivities[si] * delta.regr)
      loglik <- sum(dbinom(resp, 1, decision.prob, log = TRUE))

#      density.nolapse <- loglik + sens.prior[si] + sc.prior + lps.zero.prior
      for (li in 1:len(lapses)) {
        loglik.mix <- logSumExp(log(1-lapses[li])+loglik, log(lapses[li]) + 0.5)

        density <- loglik.mix + sens.prior[si] + lps.prior[li] + sc.prior
      }
    }
  }
  
  densities
}
