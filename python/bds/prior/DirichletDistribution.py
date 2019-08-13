from .Prior import Prior

class DirichletDistribution(Prior):
  def __init__(self, param_name, dims, defaults=dict()):
    super().__init__(param_name, defaults)

    self.dims = dims
    self.build_prior_code()

  def build_prior_code(self):
    self.param_code = "  simplex[" + str(self.dims) + "] " + self.name + ";"
    self.model_code = "  " + self.name + " ~ dirichlet(rep_vector(" + self.name + "Alpha, " + str(self.dims) + "));"
    self.data_code = "  real<lower=0> " + self.name + "Alpha;"

  def register(self, bdsmodel):
    super().register(bdsmodel)

    bdsmodel._code_blocks['data'].append(self.data_code)
