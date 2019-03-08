library(rstan)

order_data <- function(mlds_data) {
  if (ncol(mlds_data) == 4) {
    new.df <- data.frame(S1=numeric(), S2=numeric(), S3=numeric(), Response=numeric())
    for (n in 1:nrow(mlds_data)) {
      if (mlds_data[n,1] > mlds_data[n,3]) {
        new.df <- rbind(new.df, data.frame(S1=mlds_data[n,3], S2=mlds_data[n,2], S3=mlds_data[n,1], Response=1-mlds_data[n,4]))
      } else {
        new.df <- rbind(new.df, data.frame(S1=mlds_data[n,1], S2=mlds_data[n,2], S3=mlds_data[n,3], Response=mlds_data[n,4]))
      }
    }
  } else {
    new.df <- data.frame(S1=numeric(), S2=numeric(), S3=numeric(), S4=numeric(), Response=numeric())
    for (n in 1:nrow(mlds_data)) {
      if (mlds_data[n,1] > mlds_data[n,2]) {
        s2 <- mlds_data[n,1]
        s1 <- mlds_data[n,2]
      } else {
        s1 <- mlds_data[n,1]
        s2 <- mlds_data[n,2]
      }
      if (mlds_data[n,3] > mlds_data[n,4]) {
        s4 <- mlds_data[n,3]
        s3 <- mlds_data[n,4]
      } else {
        s3 <- mlds_data[n,3]
        s4 <- mlds_data[n,4]
      }
      if (mlds_data[n,1] > mlds_data[n,3]) {
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
                precLowest=0,
                precLow=2,
                precHigh=20,
                precHighest=30,
                lpsAlpha=1,
                lpsBeta=5) {

  mlds_data.ordered <- order_data(mlds_data)

  data = list(
    N = length(mlds_data[,1]),
    precLowest = precLowest,
    precLow = precLow,
    precHigh = precHigh,
    precHighest = precHighest,
    lpsAlpha = lpsAlpha,
    lpsBeta = lpsBeta
  )

  if (ncol(mlds_data) == 4) {
    data$S1 <- mlds_data.ordered[,1]
    data$S2 <- mlds_data.ordered[,2]
    data$S3 <- mlds_data.ordered[,2]
    data$S4 <- mlds_data.ordered[,3]
    data$Responses <- mlds_data.ordered[,4]
    data$K <- max(mlds_data.ordered[,1:3])
  } else if(ncol(mlds_data) == 5) {
    data$S1 <- mlds_data.ordered[,1]
    data$S2 <- mlds_data.ordered[,2]
    data$S3 <- mlds_data.ordered[,3]
    data$S4 <- mlds_data.ordered[,4]
    data$Responses <- mlds_data.ordered[,5]
    data$K <- max(mlds_data.ordered[,1:4])
  } else {
    stop("Difference scaling data should have 4 or 5 columns!")
  }

  if (is.null(stimulus)) {
    stimulus <- seq(0,1, len=data$K)
  }

  pkg_folder <- system.file(package="bds")

  if (fit.lapses) {
    data$lpsAlpha <- lpsAlpha
    data$lpsBeta <- lpsBeta
    #stan.file <- system.file("stan", "bds_lps.stan", package="bds")
    modelfile <- "/stan/models/bds_lps.stan"
  } else {
    #stan.file <- system.file("stan", "bds.stan", package="bds")
    modelfile <- "/stan/models/bds.stan"
  }

  stan.file <- paste(pkg_folder, modelfile, sep="")

  bds_model <- stan_model(file=stan.file)
  fit <- sampling(bds_model,
                  data=data,
                  iter = 2000, chains=4,
#                  include = FALSE, pars = c("psi_ext", "decision"),
                  control = list(adapt_delta = 0.99))

  summ <- rstan::summary(fit, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))$summary

  if (fit.lapses) {
    lapserate <- summ['lapses', 'mean']
    lps_summ <- summ['lapses', ]

    result <- list(
      stanfit = fit,
      stimulus = stimulus,
      scale = c(0.0, summ[paste0('psi[', 1:(data$K-2),']'),'mean'], 1.0),
      precision = summ['precision','mean'],
      lapserate = lapserate,
      scale_summary = rbind(rep(0, times=ncol(summ)),
                            summ[paste0('psi[', 1:(data$K-2),']'),],
                            rep(1, times=ncol(summ))),
      prec_summary = summ['precision',],
      lps_summary = lps_summ,
      data = data
    )

    class(result) <- "LpsDifferenceScale"

  } else {
    result <- list(
      stanfit = fit,
      stimulus = stimulus,
      scale = c(0.0, summ[paste0('psi[', 1:(data$K-2),']'),'mean'], 1.0),
      precision = summ['precision','mean'],
      scale_summary = rbind(rep(0, times=ncol(summ)),
                            summ[paste0('psi[', 1:(data$K-2),']'),],
                            rep(1, times=ncol(summ))),
      prec_summary = summ['precision',],
      data = data
    )

    class(result) <- "DifferenceScale"

  }

  result
}
