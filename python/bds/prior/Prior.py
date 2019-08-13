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
