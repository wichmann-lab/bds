from .Prior import Prior

class BetaDistribution(Prior):
  def __init__(self, param_name, defaults=dict()):
    super().__init__(param_name, defaults)

    self.build_prior_code()

  def build_prior_code(self):
    self.param_code = "  real<lower=0, upper=1> " + self.name + ";"
    self.model_code = "  " + self.name + " ~ beta(" + self.name + "Alpha, " + self.name + "Beta);"
    self.data_code = "  real<lower=0> " + self.name + "Alpha;\n  real<lower=0> " + self.name + "Beta;"

  def register(self, bdsmodel):
    super().register(bdsmodel)

    bdsmodel._code_blocks['data'].append(self.data_code)
