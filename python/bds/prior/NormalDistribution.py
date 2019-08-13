from .Prior import Prior

class NormalDistribution(Prior):
  def __init__(self, param_name, defaults=dict()):
    super().__init__(param_name, defaults)

    self.build_prior_code()

  def build_prior_code(self):
    self.param_code = "  real " + self.name + ";"
    self.model_code = "  " + self.name + " ~ normal(" + self.name + "Mu, " + self.name + "Sigma);"
    self.data_code = "  real " + self.name + "Mu;\n  real<lower=0> " + self.name + "Sigma;"

  def register(self, bdsmodel):
    super().register(bdsmodel)

    bdsmodel._code_blocks['data'].append(self.data_code)

class HalfNormalDistribution(Prior):
  def __init__(self, param_name, defaults=dict()):
    super().__init__(param_name, defaults)

    self.build_prior_code()

  def build_prior_code(self):
    self.param_code = "  real<lower=0> " + self.name + ";"
    self.model_code = "  " + self.name + " ~ normal(0, " + self.name + "Sigma);"
    self.data_code = "  real<lower=0> " + self.name + "Sigma;"

  def register(self, bdsmodel):
    super().register(bdsmodel)

    bdsmodel._code_blocks['data'].append(self.data_code)

class VectorizedNormalDistribution(Prior):
  def __init__(self, param_name, dims, defaults=dict()):
    super().__init__(param_name, defaults)

    self.dims = dims
    self.build_prior_code()

  def build_prior_code(self):
    self.param_code = "  vector[" + str(self.dims) + "] " + self.name + ";"
    self.model_code = "  " + self.name + " ~ normal(" + self.name + "Mu, " + self.name + "Sigma);"
    self.data_code = "  real " + self.name + "Mu;\n  real<lower=0> " + self.name + "Sigma;"

  def register(self, bdsmodel):
    super().register(bdsmodel)

    bdsmodel._code_blocks['data'].append(self.data_code)
