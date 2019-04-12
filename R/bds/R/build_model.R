bds.model <- list(
  fn = " /**
 * Convert four stimulus ranks, corresponding to one MLDS trial,
  * to a row of the design matrix. The stimulus rank corresponds
  * to the index in the row vector. The value depends on the position
  * of each value in the ordered quadtruple.
  *
  * @param a left interval start or end point
  * @param b left interval start or end point
  * @param c right interval start or end point
  * @param d right interval start or end point
  * @param K number of unique stimulus values shown in experiment
  * @return row of design matrix with length K-1.
  */
  row_vector order_trial(int a, int b, int c, int d, int K) {
  row_vector[K] rvec;
  rvec = rep_row_vector(0.0, K);

  if(a < b) {
  if (c < d) {
  if (a < c) {
  rvec[a] += 1.0;
  rvec[b] -= 1.0;
  rvec[c] -= 1.0;
  rvec[d] += 1.0;
  } else {
  rvec[c] += 1.0;
  rvec[d] -= 1.0;
  rvec[a] -= 1.0;
  rvec[b] += 1.0;
  }
  } else {
  if (a < d) {
  rvec[a] += 1.0;
  rvec[b] -= 1.0;
  rvec[d] -= 1.0;
  rvec[c] += 1.0;
  } else {
  rvec[d] += 1.0;
  rvec[c] -= 1.0;
  rvec[a] -= 1.0;
  rvec[b] += 1.0;
  }
  }
  } else {
  if (c < d) {
  if (b < c) {
  rvec[b] += 1.0;
  rvec[a] -= 1.0;
  rvec[c] -= 1.0;
  rvec[d] += 1.0;
  } else {
  rvec[c] += 1.0;
  rvec[d] -= 1.0;
  rvec[b] -= 1.0;
  rvec[a] += 1.0;
  }
  } else {
  if (b < d) {
  rvec[b] += 1.0;
  rvec[a] -= 1.0;
  rvec[d] -= 1.0;
  rvec[c] += 1.0;
  } else {
  rvec[d] += 1.0;
  rvec[c] -= 1.0;
  rvec[b] -= 1.0;
  rvec[a] += 1.0;
  }
  }
  }

  // Since the first scale value is zero, we discard the first column
  return rvec[2:K];
  }

  /**
  * If the stimulus ranks were arranged by position in the experiment,
  * i.e. left or right, up or down, some of the encoded responses have
  * to be flipped, since we now sort the intervals by stimulus rank order.
  *
  * @param L1 array of stimulus ranks for left interval
  * @param L2 array of stimulus ranks for left interval
  * @param R1 array of stimulus ranks for right interval
  * @param R2 array of stimulus ranks for right interval
  * @param Resp array of recorded observer responses
  * @return array of observer responses corrected for interval order
  */
  int[] flip_responses(int[] L1, int[] L2, int[] R1, int[] R2, int[] Resp) {
  int R[size(Resp)];

  for (n in 1:size(Resp)) {

  if (L1[n] > R1[n] || L1[n] > R2[n]) {
  R[n] = 1 - Resp[n];
  } else {
  R[n] = Resp[n];
  }
  }

  return R;
  }

  /**
  * Rearrange the data to fit the framework of a GLM.
  *
  * @param N Number of rows corresponding to data items
  * @param K Number of predictors, counting the intercept, per
  *          item.
  * @return Simulated predictor matrix.
  */
  matrix construct_design_matrix(int[] L1, int[] L2, int[] R1, int[] R2, int K) {
  matrix[size(L1), K-1] X;

  X = rep_matrix(0.0, size(L1), K-1);

  for (n in 1:size(L1)) {
  int a = L1[n];
  int b = L2[n];
  int c = R1[n];
  int d = R2[n];

  X[n] = order_trial(a,b,c,d,K);
  }

  return X;
  }\n",

  data = "  int<lower=1> N;  // total number of observations
  int<lower=1> K;  // number of stimulus levels
  int<lower=1, upper=K> S1[N];
  int<lower=1, upper=K> S2[N];
  int<lower=1, upper=K> S3[N];
  int<lower=1, upper=K> S4[N];
  int<lower=0, upper=1> Responses[N];  // response variable\n",

  trdata = "  matrix[N, K-1] X;
  int<lower=0, upper=1> Resp[N];

  X = construct_design_matrix(S1, S2, S3, S4, K);

  // Change response encoding from left or right to upper or lower interval
  Resp = flip_responses(S1, S2, S3, S4, Responses);\n",

  model = "  // likelihood including all constants
  for (n in 1:N) {
    target += log_mix(lapses,
                      bernoulli_lpmf(Resp[n] | 0.5),
                      bernoulli_lpmf(Resp[n] | decision[n]));
  }\n",

  gen = "  vector[N] log_lik;
  vector[N] lapse_log_lik;
  vector[N] model_log_lik;
  vector[N] log_lik_sat;

  vector[N] log_lik_hat;
  vector[N] lapse_log_lik_hat;
  vector[N] model_log_lik_hat;
  vector[N] log_lik_sat_hat;

  int resp_hat[N];
  int lapse_hat[N];

  for (n in 1:N) {
  lapse_log_lik[n] = bernoulli_lpmf(Resp[n] | 0.5);
  model_log_lik[n] = bernoulli_lpmf(Resp[n] | decision[n]);

  log_lik[n] = log_mix(lapses,
  lapse_log_lik[n],
  model_log_lik[n]);

  log_lik_sat[n] = log_mix(lapses, lapse_log_lik[n], 0);

  lapse_hat[n] = bernoulli_rng(lapses);
  if (lapse_hat[n] == 0) {
  resp_hat[n] = bernoulli_rng(decision[n]);
  } else {
  resp_hat[n] = bernoulli_rng(0.5);
  }
  lapse_log_lik_hat[n] = bernoulli_lpmf(resp_hat[n] | 0.5);
  model_log_lik_hat[n] = bernoulli_lpmf(resp_hat[n] | decision[n]);
  log_lik_hat[n] = log_mix(lapses,
  lapse_log_lik[n],
  model_log_lik[n]);
  log_lik_sat_hat[n] = log_mix(lapses, lapse_log_lik[n], 0);
  }\n"
  )

prec.raised_cosine <- list(
  fn = "/**
 * Return a data matrix of specified size with rows
 * corresponding to items and the first column filled
 * with the value 1 to represent the intercept and the
 * remaining columns randomly filled with unit-normal draws.
 *
 * @param N Number of rows corresponding to data items
 * @param K Number of predictors, counting the intercept, per
 *          item.
 * @return Simulated predictor matrix.
 */
  real raised_cosine_lpdf(real y, real start, real u1, real u2, real end) {
    real res;
    real norm;
    if (y < u1) {
      res = log(0.5-0.5*cos(pi()/(u1-start)*(y-start)));
    } else if (y > u2) {
      res = log(0.5+0.5*cos(pi()/(end-u2)*(y-u2)));
    } else if (u1 <= y && y <= u2) {
      res = 0.0;
    } else {
      reject(\"input must be inside the specified range!\");
    }
    norm = log(0.5*(u1-start+end-u2) + (u2-u1));

    return res-norm;
  }\n",

  data =  " // hyper-parameters for precision prior (uniform with cosine falloff)
  real precLowest;
  real<lower=precLowest> precHighest;
  real<lower=precLowest, upper=precHighest> precLow;
  real<lower=precLow, upper=precHighest> precHigh;\n",
  trdata = "",
  param = "  real<lower=precLowest, upper=precHighest> precision;\n",
  trparam = "",
  model = "  precision ~ raised_cosine(precLowest, precLow, precHigh, precHighest);\n",
  default = list(precLowest=0, precLow=2, precHigh=15, precHighest=30)
  )

prec.halfgauss <- list(
  fn = "",
  data = " // hyper-parameters for precision prior (half-normal)
  real<lower=0> precSigma;\n",
  trdata = "",
  param = "  real<lower=0> precision;\n",
  trparam = "",
  model = "  precision ~ normal(0, precSigma);\n",
  default = list(precSigma=15)
)

prec.uniform <- list(
  fn = "",
  data = " // hyper-parameters for precision prior (uniform)
  real<lower=0> precLow;
  real<lower=precLow> precHigh;\n",
  trdata = "",
  param = "  real<lower=precLow, upper=precHigh> precision;\n",
  trparam = "",
  model = "",
  default = list(precLow=0, precHigh=30)
)

prec.const <- list(
  fn = "",
  data = "  real<lower=0> precision;\n",
  trdata = "",
  param = "",
  trparam = "",
  model = "",
  default = list(precision = 10)
)

psi.uniform <- list(
  fn = "",
  data = "",
  trdata = "",
  param = "  vector<lower=0,upper=1>[K-2] psi;\n",

  trparam = "  vector[N] decision;
  vector[K-1] psi_ext;

  for (k in 1:K-2) {
    psi_ext[k] = psi[k];
  }
  psi_ext[K-1] = 1;

  decision = Phi(X * psi_ext * precision);\n",

  model = "",
  default = list()
  )

lapses.uniform <- list(
  fn = "",
  data = "",
  trdata = "",
  param = "  real<lower=0, upper=1> lapses;\n",
  trparam = "",
  model = "",
  default = list()
)

lapses.const <- list(
  fn = "",
  data = "  real<lower=0, upper=1> lapses;",
  trdata = "",
  param = "",
  trparam = "",
  model = "",
  default = list(lapses=0.0)
)

lapses.beta <- list(
  fn = "",
  data = "  // hyper-parameters for lapse rate prior (beta distribution)
  real<lower=0> lpsAlpha;
  real<lower=0> lpsBeta;\n",
  trdata = "",
  param = "  real<lower=0, upper=1> lapses;\n",
  trparam = "",
  model = "  lapses ~ beta(lpsAlpha, lpsBeta);\n",
  default = list(lpsAlpha=1, lpsBeta=5)
)

build_model <- function(priors, model, extractor_function) {
  fn.block <- paste0("functions {\n",
                     paste0(sapply(priors, function(pr) pr$fn), collapse = ""),
                     model$fn,
                     "}\n\n")

  data.block <- paste0("data {\n",
                       model$data,
                       paste0(sapply(priors, function(pr) pr$data), collapse = ""),
                       "}\n\n")

  data.transformed.block <- paste0("transformed data {\n",
                                   model$trdata,
                                   paste0(sapply(priors, function(pr) pr$trdata), collapse = ""),
                                   "}\n\n")

  param.block <- paste0("parameters {\n",
                        paste0(sapply(priors, function(pr) pr$param), collapse = ""),
                        "}\n\n")

  param.transformed.block <- paste0("transformed parameters {\n",
                                    paste0(sapply(priors, function(pr) pr$trparam), collapse = ""),
                                    "}\n\n")

  model.block <- paste0("model {\n",
                        paste0(sapply(priors, function(pr) pr$model), collapse = ""),
                        model$model,
                        "}\n\n")

  generated.block <- paste0("generated quantities {\n",
                            model$gen,
                            "}\n")

  model_code <- paste0(fn.block, data.block, data.transformed.block, param.block, param.transformed.block, model.block, generated.block)

  default_params = do.call(c, lapply(priors, function(pr) pr$default))

  list(model_code=model_code, default_params=default_params, extractor=extractor_function)
}
