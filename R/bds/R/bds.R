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
                fit.lapses=TRUE) {

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
    md <- build_model(priors=list(psi.uniform, prec.raised_cosine, lapses.beta),
                model=bds.model,
                extractor_function = default_extractor)
    init_fun <- function() {
      list(psi = stimulus[2:(length(stimulus)-1)],
           precision = md$default_params$precLow,
           lapses = 0.01)
    }
  } else {
    md <- build_model(priors=list(psi.uniform, prec.raised_cosine, lapses.const),
                      model=bds.model,
                      extractor_function = extractor_fixed_lapserate)
    init_fun <- function() {
      list(psi = stimulus[2:(length(stimulus)-1)],
           precision = (md$default_params$precLow + md$default_params$precHigh)/2.0)
    }
  }

  model_obj <- stan_model(model_code=md$model_code)

  stanfit <- sample_bds_model(model_obj, mlds_data, prior_params=md$default_params, init_list = rep(list(init_fun()), times=4))

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
                  control = list(adapt_delta = 0.99),
                  init=init_list,
                  cores=.cores)

  list(stanfit=fit, data=data)
}
