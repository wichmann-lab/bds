library(bds)
library(ggplot2)

options(mc.cores=parallel::detectCores())

raw_data <- read.table("test.csv", sep=" ", header = TRUE)

mlds_data <- data.frame(Resp=raw_data$Response, S1=raw_data$i1, S2=raw_data$i2, S3=raw_data$i3)
stimulus <- sort(unique(c(raw_data$s1, raw_data$s2, raw_data$s3)))

scale <- bds(mlds_data, fit.lapses = FALSE)

plot.df <- data.frame(stimulus=scale$stimulus,
                      scale=get_scale_values(scale),
                      ci.low=get_scale_credible_interval(scale)$ci.low,
                      ci.high=get_scale_credible_interval(scale)$ci.high)

ggplot(plot.df, aes(x=stimulus,y=scale)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(x=stimulus,ymin=ci.low,ymax=ci.high))

diagnostic_plots(scale)

# Build and use model with different priors than the default models

new.model <- build_model(priors=list(psi.uniform, prec.halfgauss, lapses.const),
                          model=bds.model,
                          extractor_function = extractor_fixed_lapserate)

# Code for model with lapserate prior
# new.model <- build_model(priors=list(psi.uniform, prec.halfgauss, lapses.beta),
#                         model=bds.model,
#                         extractor_function = default_extractor)
#
# init_fun <- function() {
#   list(psi = stimulus[2:(length(stimulus)-1)],
#     precision = (md$default_params$precLow + md$default_params$precHigh)/2.0,
#     lapses = 0.01)
# }

new.model.obj <- stan_model(model_code = new.model$model_code)

# Set the prior parameters
params <- new.model$default_params

params$lapses <- 0.025
params$precSigma <- 10

init_fun <- function() {
  list(psi = stimulus[2:(length(stimulus)-1)],
       precision = params$precSigma/2.0)
}

stanfit <- sample_bds_model(new.model.obj, mlds_data, prior_params=params, init_list = rep(list(init_fun()), times=4))

new.scale <- new.model$extractor(stanfit$stanfit, stimulus, stanfit$data)
