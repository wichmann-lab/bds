library(rstan)

DifferenceScale <- setClass("DifferenceScale", slots = list(
  stanfit = "stanfit",
  stimulus = "numeric",
  scale = "numeric",
  precision = "numeric",
  lapserate = "numeric",
  scale_summary = "matrix",
  prec_summary = "numeric",
  lps_summary = "numeric"
))

bds <- function(mlds_data,
                stimulus=NULL,
                fit.lapses=TRUE,
                precLowest=0,
                precLow=2,
                precHigh=20,
                precHighest=30,
                lpsAlpha=1,
                lpsBeta=5) {

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
    data$S1 <- mlds_data[,1]
    data$S2 <- mlds_data[,2]
    data$S3 <- mlds_data[,2]
    data$S4 <- mlds_data[,3]
    data$Responses <- mlds_data[,4]
    data$K <- max(mlds_data[,1:3])
  } else if(ncol(mlds_data) == 5) {
    data$S1 <- mlds_data[,1]
    data$S2 <- mlds_data[,2]
    data$S3 <- mlds_data[,3]
    data$S4 <- mlds_data[,4]
    data$Responses <- mlds_data[,5]
    data$K <- max(mlds_data[,1:4])
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
                  include = FALSE, pars = c("psi_ext", "decision"),
                  control = list(adapt_delta = 0.99))

  summ <- rstan::summary(fit, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))$summary

  if (fit.lapses) {
    lapserate <- summ['lapses', 'mean']
    lps_summ <- summ['lapses', ]
  } else {
    lapserate <- 0
    lps_summ <- rep(0, times=ncol(summ))
  }

  result <- DifferenceScale(
    stanfit = fit,
    stimulus = stimulus,
    scale = c(0.0, summ[paste0('psi[', 1:(data$K-2),']'),'mean'], 1.0),
    precision = summ['precision','mean'],
    lapserate = lapserate,
    scale_summary = rbind(rep(0, times=ncol(summ)),
                          summ[paste0('psi[', 1:(data$K-2),']'),],
                          rep(1, times=ncol(summ))),
    prec_summary = summ['precision',],
    lps_summary = lps_summ
  )
}
