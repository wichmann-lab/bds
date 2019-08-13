from .Prior import Prior

class Constant(Prior):
  def __init__(self, param_name, type_sig, defaults=dict()):
    super().__init__(param_name, defaults)
    self.type_sig = type_sig

    self.build_prior_code()

  def build_prior_code(self):
    self.data_code = "  " + self.type_sig + " " + self.name + ";"

  def register(self, bdsmodel):
    bdsmodel._code_blocks['data'].append(self.data_code)

    for (k,v) in self.defaults.items():
      bdsmodel.hyper_params[k] = v
