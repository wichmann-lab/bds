from dfply import *

class LikelihoodModel:
  def __init__(self):
    self.model_code = ""

  def register(self, bdsmodel):
    bdsmodel._code_blocks['model'].append(self.model_code)

  def results(self, result_obj):
    result_obj.residuals = self.residuals(result_obj)
    result_obj.ppc_residuals = self.ppc_residuals(result_obj)
    result_obj.deviance = self.deviance(result_obj.stan_fit)
    result_obj.ppc_deviance = self.ppc_deviance(result_obj.stan_fit)

  def residuals(self, result_obj):
    res_pars = ['resid[%d]' % x for x in range(1,result_obj.n+1)]
    residuals = (result_obj.stan_fit.to_dataframe(pars=res_pars, 
                                  inc_warmup=False,
                                  diagnostics=False)
                 >> select(res_pars)
                ).values

    return residuals

  def ppc_residuals(self, result_obj):
    res_pars = ['resid_hat[%d]' % x for x in range(1,result_obj.n+1)]
    residuals = (result_obj.stan_fit.to_dataframe(pars=res_pars, 
                                  inc_warmup=False,
                                  diagnostics=False)
                 >> select(res_pars)
                ).values

    return residuals

  def deviance(self, stan_fit):
    return stan_fit.extract(pars='deviance')

  def ppc_deviance(self, stan_fit):
    return stan_fit.extract(pars='deviance_hat')
