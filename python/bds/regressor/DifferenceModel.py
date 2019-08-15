from .RegressorModel import RegressorModel

class DifferenceModel(RegressorModel):
  def __init__(self, link, sensitivity_prior):
    super().__init__(link)
    self.sns_prior = sensitivity_prior

    self.param_transform_def = "  vector[N] difference;"
    self.param_transform = (
"""
  for (n in 1:N) {
    difference[n] = fabs(psi[S4[n]] - psi[S3[n]]) - fabs(psi[S2[n]] - psi[S1[n]]);
  }""")

    self.formula = "difference * sensitivity"

    self.build_code()

  def register(self, bdsmodel):
    self.sns_prior.register(bdsmodel)

    bdsmodel._code_blocks['transformed parameters:definitions'].append(self.param_transform_def)
    bdsmodel._code_blocks['transformed parameters'].append(self.param_transform)

    super().register(bdsmodel)

  def results(self, result_obj):
    result_obj.sensitivity = self.sensitivity(result_obj.stan_fit)

  def sensitivity(self, stan_fit):
    return stan_fit.extract(pars=self.sns_prior.name)

class BESSDifferenceModel:
  def __init__(self, sensitivity_prior):
    self.sns_prior = sensitivity_prior

    self.param_transform_def = "  vector[N] decision;"

    self.param_transform = "  decision = (psi[U] - psi[L])/2.0;"

  def register(self, bdsmodel):
    self.sns_prior.register(bdsmodel)

    bdsmodel._code_blocks['transformed parameters:definitions'].append(self.param_transform_def)
    bdsmodel._code_blocks['transformed parameters'].append(self.param_transform)

  def results(self, result_obj):
    result_obj.sensitivity = self.sensitivity(result_obj.stan_fit)

  def sensitivity(self, stan_fit):
    return stan_fit.extract(pars=self.sns_prior.name)

