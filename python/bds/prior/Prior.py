class Prior:
  def __init__(self, param_name, defaults=dict()):
    self.name = param_name
    self.defaults = defaults

    self.param_code = ""
    self.model_code = ""

  def build_prior_code(self):
    pass

  def register(self, bdsmodel):
    bdsmodel._code_blocks['parameters'].append(self.param_code)
    bdsmodel._code_blocks['model'].append(self.model_code)

    for (k,v) in self.defaults.items():
      bdsmodel.hyper_params[k] = v

class UnconstrainedVector(Prior):
  def __init__(self, param_name, dims):
    super().__init__(param_name)

    self.dims = dims
    self.build_prior_code()

  def build_prior_code(self):
    self.param_code = "  vector[" + str(self.dims) + "] " + self.name + ";"

class UnconstrainedReal(Prior):
  def __init__(self, param_name):
    super().__init__(param_name)

    self.dims = dims
    self.build_prior_code()

  def build_prior_code(self):
    self.param_code = "  real " + self.name + ";"
