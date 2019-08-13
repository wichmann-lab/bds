from .ScaleModel import ScaleModel

from dfply import *

class GaussianProcess(ScaleModel):
  def __init__(self, length_scale_prior, magnitude_prior):
    super().__init__()
    self.length_scale_prior = length_scale_prior
    self.magnitude_prior = magnitude_prior

    self.data_code = (
"""  int<lower=1> N_predict;
  real x_predict[N_predict];
  int<lower=1, upper=N_predict> obs_idx[K-2];""")

    self.param = "  vector[N_predict-2] psi_tilde;"

    self.data_transform_def = "  vector[2] fix_points = [0, 1]';\n  int idx[2] = {1, N_predict};"

    self.param_transform_gp_def = (
"""  matrix[N_predict, N_predict] cov_all = cov_exp_quad(x_predict, """ + self.magnitude_prior.name + """, """ + self.length_scale_prior.name + """)
                     + diag_matrix(rep_vector(1e-10, N_predict));

  matrix[N_predict-2, 2] cov_ub = cov_all[2:N_predict-1, idx];
  matrix[2, 2] cov_bb = cov_all[idx, idx];

  vector[N_predict-2] y_mean = (cov_ub / cov_bb) * fix_points;

  matrix[N_predict-2, N_predict-2] L_cov = cholesky_decompose(cov_all[2:N_predict-1, 2:N_predict-1] - (cov_ub / cov_bb) * cov_ub');
  vector[N_predict-2] psi_predict = y_mean + L_cov * psi_tilde;""")

    self.param_transform = (
"""  psi[1] = 0;
  psi[K] = 1;
  for (k in 2:K-1) {
    psi[k] = psi_predict[obs_idx[k-1]-1];
  }""")

    self.model_code = "  psi_tilde ~ normal(0, 1);"

  def register(self, bdsmodel):
    self.length_scale_prior.register(bdsmodel)
    self.magnitude_prior.register(bdsmodel)

    super().register(bdsmodel)

    bdsmodel._code_blocks['data'].append(self.data_code)
    bdsmodel._code_blocks['transformed data:definitions'].append(self.data_transform_def)
    bdsmodel._code_blocks['parameters'].append(self.param)
    bdsmodel._code_blocks['transformed parameters:definitions'].append(self.param_transform_gp_def)
    bdsmodel._code_blocks['transformed parameters'].append(self.param_transform)
    bdsmodel._code_blocks['model'].append(self.model_code)

  def results(self, result_obj):
    result_obj.scale = self.scale_values(result_obj)
    result_obj.stimulus = result_obj.stan_data['x_predict']

  def scale_values(self, result_obj):
    scale_pars = ['psi[1]'] + ['psi_predict[%d]' % x for x in range(1,result_obj.stan_data['N_predict']-1)] + ['psi[%d]' % result_obj.k]
    scale_vals = (result_obj.stan_fit.to_dataframe(pars=scale_pars, 
                                  inc_warmup=False,
                                  diagnostics=False)
                 >> select(scale_pars)
                ).values

    return scale_vals
