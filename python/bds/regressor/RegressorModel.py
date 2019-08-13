class RegressorModel:
  def __init__(self, link):
    self.link_function = link
    self.formula = ""

  def build_code(self):
    self.code = self.link_function.embed(self.formula)

  def register(self, bdsmodel):
    self.link_function.register(bdsmodel)
    bdsmodel._code_blocks['transformed parameters'].append(self.code)

  def results(self, result_obj):
    pass
