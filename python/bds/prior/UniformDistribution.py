from .Prior import Prior

class UniformDistribution(Prior):
  def __init__(self, param_name, defaults=dict()):
    super().__init__(param_name, defaults)

    self.build_prior_code()

  def build_prior_code(self):
    self.data_code = "  real " + self.name + "Start;\n  real " + self.name + "End;"
    self.param_code = "  real<lower=" + self.name + "Start, upper=" + self.name + "End> " + self.name + ";"

  def register(self, bdsmodel):
    super().register(bdsmodel)

    bdsmodel._code_blocks['data'].append(self.data_code)

class VectorizedUniformDistribution(Prior):
  def __init__(self, param_name, dims, defaults=dict()):
    super().__init__(param_name, defaults)

    self.dims = dims
    self.build_prior_code()

  def build_prior_code(self):
    self.data_code = "  real " + self.name + "Start;\n  real " + self.name + "End;"
    self.param_code = "  vector<lower=" + self.name + "Start, upper=" + self.name + "End>[" + str(self.dims) + "] " + self.name + ";"

  def register(self, bdsmodel):
    super().register(bdsmodel)

    bdsmodel._code_blocks['data'].append(self.data_code)
