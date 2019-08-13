class LinkFunction:
  def __init__(self):
    self.param_transform_def = "  vector[N] decision;"

  def embed(self, formula):
    return formula

  def register(self, bdsmodel):
    bdsmodel._code_blocks['transformed parameters:definitions'].append(self.param_transform_def)

class ProbitLink(LinkFunction):
  def embed(self, formula):
    return '  decision = Phi(' + formula + ');'
