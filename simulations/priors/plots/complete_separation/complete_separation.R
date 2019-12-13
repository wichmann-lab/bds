library(ggplot2)

prob.complete.sep <- function(num_intens, num_repl, power, sensitivity, thr=.Machine$double.eps) {
  num_triads <- choose(num_intens, 3)
  Tr <- do.call(rbind, replicate(num_repl, t(combn(num_intens,3)), simplify=FALSE))
  
  beta <- seq(0,1,len=num_intens)^power
  
  dec <- mapply(function(a, b, c) abs(beta[a] - 2*beta[b] + beta[c]),
                Tr[,1], Tr[,2], Tr[,3])
  # since pnorm(a) = 1-pnorm(-a), we can use pnorm(abs(a))
  # we need to exclude triads, where the linear predictor is exactly zero
  # due to numerical error we use a threshold cutoff
  log_pr_vec <- log10(pnorm(sensitivity*dec[dec > thr]))
  
  10^sum(log_pr_vec)
}

find.sensitivity.threshold <- function(num_intens, num_repl, power) {
  f <- function(x) 0.5 - prob.complete.sep(num_intens, num_repl, power, x)#, thr=1e-2)
  uniroot(f, c(0, 1e6))$root
}

get.alpha <- function(num_intens, power, correct_id=FALSE) {
  beta <- seq(0,1,len=num_intens)^power
  if (power == 1.0 && correct_id) {
    alpha <- 1.0
  } else {
    alpha <- max(abs(diff(diff(beta))))
  }
  alpha
}

get.omega <- function(num_intens, power, correct_id=FALSE) {
  beta <- seq(0,1,len=num_intens)^power
  if (power == 1.0 && correct_id) {
    omega <- 1.0
  } else {
    omega <- min(abs(diff(diff(beta))))
  }
  omega
}

sens <- seq(1.0, 100.0, len=100.0)
probs <- mapply(prob.complete.sep, 10, 1, 1.0, sens)
qplot(sens, probs) + geom_line() + xlab("sensitivity") + ylab("probability of complete separation")

pow <- c(1.0/seq(1.0, 1.5, len=100), seq(1.0, 1.5, len=100))
sens_thr <- mapply(find.sensitivity.threshold, 10, 1, pow)
qplot(pow, sens_thr) + geom_line() + xlab('exponent of power function') + ylab('sensitivity threshold')

alpha <- mapply(get.alpha, 10, pow)
omega <- mapply(get.omega, 10, pow, TRUE)
alpha_c <- mapply(get.alpha, 10, pow, TRUE)
qplot(pow, alpha) + geom_line()

sens_lb <- qnorm(0.5^(1/choose(10, 3)))/alpha_c
sens_ub <- qnorm(0.5^(1/choose(10, 3)))/omega
qplot(pow, sens_lb)
qplot(pow, sens_ub)
