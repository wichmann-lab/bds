functions {
/**
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
      reject("input must be inside the specified range!");
    }
    norm = log(0.5*(u1-start+end-u2) + (u2-u1));
      
    return res-norm;
  }

 /**
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
  }
}

/**
 * This model implements a Bayesian analysis of 
 * Maximum-likelihood difference scaling.
 * Experimental observers have to compare two intervals and
 * indicate which one is larger in the perceived dimension.
 *
 * @param N   Number of trials
 * @param K   Number of stimulus levels shown
 * @param S1, S2
 *            stimulus ranks of the left interval
 * @param S3, S4
 *            stimulus ranks of the right interval
 *            for triads, set S3=S2
 * @param Responses.
 *            observer responses; 0 for left, 1 for right
 *            interval
 * @param precLowest
 *            lowest possible value of precision (prior parameter)
 * @param precLow
 *            lowest reasonable value of precision (prior parameter)
 * @param precHigh
 *            highest reasonable value of precision (prior parameter)
 * @param precHighest
 *            highest possible value of precision (prior parameter)
 * @param lpsAlpha
 *            parameter of Beta distribution prior on lapse rate
 * @param lpsBeta
 *            parameter of Beta distribution prior on lapse rate
 */
data {
  int<lower=1> N;  // total number of observations 
  int<lower=1> K;  // number of stimulus levels
  int<lower=1, upper=K> S1[N];
  int<lower=1, upper=K> S2[N];
  int<lower=1, upper=K> S3[N];
  int<lower=1, upper=K> S4[N];
  int<lower=0, upper=1> Responses[N];  // response variable
  
  // hyper-parameters for precision prior (uniform with cosine falloff)
  real precLowest;
  real<lower=precLowest> precHighest;
  real<lower=precLowest, upper=precHighest> precLow;
  real<lower=precLow, upper=precHighest> precHigh;
  
  // hyper-parameters for lapse rate prior (beta distribution)
  real<lower=0> lpsAlpha;
  real<lower=0> lpsBeta;
}

/**
 * For the generalized linear model, the stimulus ranks have to be
 * encoded as a design matrix, where the rows correspond to the trials
 * and the columns to the stimulus rank shown in that trial.
 */
transformed data {
  matrix[N, K-1] X;
  int<lower=0, upper=1> Resp[N];
  
  X = construct_design_matrix(S1, S2, S3, S4, K);

  // Change response encoding from left or right to upper or lower interval
  Resp = flip_responses(S1, S2, S3, S4, Responses);
}

/**
 * The parameters of the model consist of the scale values corresponding to
 * the shown stimulus values, the standard deviation of the decision noise,
 * and a lapse rate to explain observer errors.
 * The first and last scale value is fixed to 0 and 1, respectively.
 *
 * @param psi 
 *        vector of scale values to be estimated
 *        constrained to lay between 0 and 1
 * @param precision
 *        inverse of the standard deviation of the
 *        decision noise
 * @param lapses
 *        rate of stimulus-independent observer errors
 */
parameters {
  vector<lower=0,upper=1>[K-2] psi;
  real<lower=precLowest, upper=precHighest> precision;
  real<lower=0, upper=1> lapses;
}

/**
 * Upper scale value is fixed to 1. Fixation of the lower scale value is
 * done by dropping the first column of the design matrix.
 * The decision probability is computed by passing the decision variable
 * through the cumulative normal function.
 * The decision variable itself is the product between the design matrix
 * and the scale values scaled by the precision of the decision noise.
 */
transformed parameters {
  vector[N] decision;
  vector[K-1] psi_ext;

  for (k in 1:K-2) {
    psi_ext[k] = psi[k];
  }
  psi_ext[K-1] = 1;

  decision = Phi(X * psi_ext * precision);
}

/**
 * The scale values have an implicit uniform prior.
 * The precision of the decision noise has a custom
 * uniform prior with cosine falloff.
 *
 * The likelihood computes as the outcome of one Bernoulli
 * trial per experimental trial.
 */
model {
  // priors
  precision ~ raised_cosine(precLowest, precLow, precHigh, precHighest);
  lapses ~ beta(lpsAlpha, lpsBeta);

  // likelihood including all constants
  for (n in 1:N) {
    target += log_mix(lapses,
                      bernoulli_lpmf(Resp[n] | 0.5),
                      bernoulli_lpmf(Resp[n] | decision[n]));
  }
}

/**
 * @param log_lik
 *        can be used for model comparisons with the loo package.
 * @param lapse_log_lik 
 *        likelihood of responses assuming fair coin
 * @param model_log_lik
 *        likelihood of responses according to the decision variable
 * @param log_lik_ratio
 *        log of the likelihood ratio to classify lapses
 * @param resp_hat
 *        posterior predictive distribution samples
 */
generated quantities {
  vector[N] log_lik;
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
  }
}
