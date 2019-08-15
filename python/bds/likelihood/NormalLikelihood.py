from .LikelihoodModel import LikelihoodModel
from dfply import *

import plotnine as gg

class NormalLikelihood(LikelihoodModel):
  def __init__(self):
    self.model_code = "  psi[Response] ~ normal(decision, 1.0/sensitivity);"

    self.function_code = (
"""  real sign(real x) {
    return x < 0 ? -1 : 1;
  }""")

    self.generator_def = (
"""  vector[N] log_lik;
  vector[N] log_lik_hat;
  
  vector[N] resid;
  vector[N] resid_hat;

  real deviance;
  real deviance_hat;

  real resp_hat[N];""")

    self.generator = (
"""  for (n in 1:N) {
    log_lik[n] = normal_lpdf(psi[Response[n]] | decision[n], 1.0/sensitivity);

    resp_hat[n] = normal_rng(decision[n], 1.0/sensitivity);
    log_lik_hat[n] = normal_lpdf(resp_hat[n] | decision[n], 1.0/sensitivity);

    resid[n] = sign(psi[Response[n]] - decision[n]) * sqrt( 2 * (normal_lpdf(psi[Response[n]] | psi[Response[n]], 1.0/sensitivity) - log_lik[n]));
    resid_hat[n] = sign(resp_hat[n] - decision[n]) * sqrt( 2 * (normal_lpdf(resp_hat[n] | resp_hat[n], 1.0/sensitivity) - log_lik_hat[n]));
  }

  deviance = sum(square(resid));
  deviance_hat = sum(square(resid_hat));""")

  def register(self, bessmodel):
    super().register(bessmodel)

    bessmodel._code_blocks['functions']['sign'] = self.function_code
    bessmodel._code_blocks['generated quantities:definitions'].append(self.generator_def)
    bessmodel._code_blocks['generated quantities'].append(self.generator)

  def results(self, result_obj):
    super().results(result_obj)

    result_obj.decision_probabilities = self.decision_probabilities(result_obj)
    result_obj.ppc_responses = self.ppc_responses(result_obj)

  def decision_probabilities(self, result_obj):
    dec_pars = ['decision[%d]' % x for x in range(1,result_obj.n+1)]
    probs = (result_obj.stan_fit.to_dataframe(pars=dec_pars, 
                                  inc_warmup=False,
                                  diagnostics=False)
                 >> select(dec_pars)
                ).values

    return probs

  def ppc_responses(self, result_obj):
    resp_pars = ['resp_hat[%d]' % x for x in range(1,result_obj.n+1)]
    responses = (result_obj.stan_fit.to_dataframe(pars=resp_pars, 
                                  inc_warmup=False,
                                  diagnostics=False)
                 >> select(resp_pars)
                ).values

    return responses
