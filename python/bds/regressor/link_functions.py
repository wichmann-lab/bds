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

class SensoryNoiseProbitLink(LinkFunction):
  def embed(self, formula):
    return '  decision = Phi( (' + formula[0] + ')/sqrt(6) ).*Phi( (' + formula[1] + ')/sqrt(2) ) + Phi( -(' + formula[0] + ')/sqrt(6) ).*Phi( -(' + formula[1] + ')/sqrt(2) );'
