from .ScaleModel import ScaleModel

class ParameterScale(ScaleModel):
  def __init__(self, scale_prior):
    super().__init__()
    self.scale_prior = scale_prior

    self.param_transform = (
"""  psi[1] = 0;
  psi[K] = 1;

  for (k in 2:K-1) {
    psi[k] = """ + self.scale_prior.name + """[k-1];
  }""")

  def register(self, bdsmodel):
    self.scale_prior.register(bdsmodel)

    super().register(bdsmodel)

    bdsmodel._code_blocks['transformed parameters'].append(self.param_transform)
