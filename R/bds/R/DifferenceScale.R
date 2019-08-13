library(cowplot)

residuals.DifferenceScale <- function(scale, ppc=FALSE) {
  N <- scale$data$N

  if (ppc) {
    resp <- as.matrix(scale$stanfit, pars = paste0('resp_hat[', 1:N, ']'))
    log_lik <- as.matrix(scale$stanfit, pars = paste0('log_lik_hat[', 1:N, ']'))
  } else {
    log_lik <- as.matrix(scale$stanfit, pars = paste0('log_lik[', 1:N, ']'))
    resp <- matrix(rep(scale$data$Responses, each=nrow(log_lik)), ncol=N)
  }
  residuals <- (2 * resp - 1) * sqrt(- 2 * log_lik)
}

residuals.LpsDifferenceScale <- function(scale, ppc=FALSE) {
  N <- scale$data$N

  if (ppc) {
    resp <- as.matrix(scale$stanfit, pars = paste0('resp_hat[', 1:N, ']'))
    log_lik <- as.matrix(scale$stanfit, pars = paste0('log_lik_hat[', 1:N, ']'))
  } else {
    log_lik <- as.matrix(scale$stanfit, pars = paste0('log_lik[', 1:N, ']'))
    resp <- matrix(rep(scale$data$Responses, each=nrow(log_lik)), ncol=N)
  }
  residuals <- (2 * resp - 1) * sqrt(-2 * log_lik)
}

ppc_ordered_residuals <- function(scale) {
  resid.emp.sorted <- apply(residuals(scale), 1, sort, na.last=TRUE)
  resid.sim.sorted <- apply(residuals(scale, ppc=TRUE), 1, sort, na.last=TRUE)

  quantiles.emp <- apply(resid.emp.sorted, 1, quantile, probs=c(0.025, 0.5, 0.975), na.rm=TRUE)
  quantiles.sim <- apply(resid.sim.sorted, 1, quantile, probs=c(0.025, 0.5, 0.975), na.rm=TRUE)

  ppc.df <- rbind(data.frame(low = quantiles.emp[1,],
                             m = quantiles.emp[2,],
                             high = quantiles.emp[3,],
                             origin = 'empirical',
                             sortid = seq(0, 1, len=ncol(quantiles.emp))),
                  data.frame(low = quantiles.sim[1,],
                             m = quantiles.sim[2,],
                             high = quantiles.sim[3,],
                             origin = 'simulated',
                             sortid = seq(0, 1, len=ncol(quantiles.sim)))
                  )

  pval <- sum(quantiles.emp[,3] < quantiles.sim[,1] | quantiles.sim[,3] < quantiles.emp[,1]) / nrow(quantiles.emp)

  list(disjoint=pval, resid_cdf=ppc.df)
}

ppc_residual_run <- function(scale) {
  N <- scale$data$N

  decision <- as.matrix(scale$stanfit, pars = paste0('decision[', 1:N, ']'))
  resp_hat <- as.matrix(scale$stanfit, pars = paste0('resp_hat[', 1:N, ']'))
  resp <- matrix(rep(scale$data$Responses, each=nrow(resp_hat)), ncol=N)

  for (n in 1:nrow(resp_hat)) {
    ind <- order(decision[n,])
    resp[n,] <- resp[n,ind]
    resp_hat[n,] <- resp_hat[n, ind]
  }

  count_reversals <- function(x) {
    sum(x[1:(length(x)-1)] != x[2:length(x)])
  }
  rev.emp <- apply(resp, 1, count_reversals)
  rev.sim <- apply(resp_hat, 1, count_reversals)

  runs.df <- rbind(data.frame(reversals=rev.emp, origin='empirical'),
                  data.frame(reversals=rev.sim, origin='simulated'))

  pval = sum(rev.emp - rev.sim > 0) / length(rev.emp)

  list(pval = pval, runs=runs.df)
}

diagnostic_plots <- function(scale) {
  cdf <- ppc_ordered_residuals(scale)
  runs <- ppc_residual_run(scale)

  p1 <- ggplot(cdf$resid_cdf, aes(x=m,y=sortid,colour=origin)) +
    geom_line() +
    geom_line(aes(x=low), linetype='dotted') +
    geom_line(aes(x=high), linetype='dotted')

  p2 <- ggplot(runs$runs, aes(x=reversals, fill=origin)) + geom_histogram()

  plot_grid(p1, p2, ncol=2)
}

get_scale_values <- function(scale) {
  scale$scale
}

get_scale_credible_interval <- function(scale) {
  list(ci.low=scale$scale_summary[,'2.5%'], ci.high=scale$scale_summary[,'97.5%'])
}

get_precision <- function(scale) {
  scale$precision
}

get_precision_credible_interval <- function(scale) {
  list(ci.low=scale$prec_summary[["2.5%"]], ci.high=scale$prec_summary[["97.5%"]])
}

get_lapserate <- function(scale) {
  scale$lapserate
}

get_lapserate_credible_interval <- function(scale) {
  list(ci.low=scale$lps_summary[["2.5%"]], ci.high=scale$lps_summary[["97.5%"]])
}
