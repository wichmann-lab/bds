from dfply import *

class ScaleModel:
  def __init__(self):
    self.param_transform_def = "  vector[K] psi;"

  def register(self, bdsmodel):
    bdsmodel._code_blocks['transformed parameters:definitions'].append(self.param_transform_def)

  def results(self, result_obj):
    result_obj.scale = self.scale_values(result_obj)

  def scale_values(self, result_obj):
    scale_pars = ['psi[%d]' % x for x in range(1,result_obj.k+1)]
    scale_vals = (result_obj.stan_fit.to_dataframe(pars=scale_pars, 
                                  inc_warmup=False,
                                  diagnostics=False)
                 >> select(scale_pars)
                ).values

    return scale_vals

