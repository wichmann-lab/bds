from .ScaleModel import ScaleModel

class MonotonicScale(ScaleModel):
  def __init__(self, diff_prior):
    super().__init__()
    self.diff_prior = diff_prior

    self.param_transform = (
""" psi[1] = 0;
  psi[2] = """ + self.diff_prior.name + """[1];
  for (k in 3:K) {
    psi[k] = psi[k-1] + """ + self.diff_prior.name + """[k-1];
  }""")

  def register(self, bdsmodel):
    self.diff_prior.register(bdsmodel)

    super().register(bdsmodel)

    bdsmodel._code_blocks['transformed parameters'].append(self.param_transform)
