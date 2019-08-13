from .Prior import Prior

class RaisedCosineDistribution(Prior):
  def __init__(self, param_name, defaults):
    super().__init__(param_name, defaults)

    self.build_prior_code()

  def build_prior_code(self):
    self.param_code = "  real<lower=" + self.name + "Lowest, upper=" + self.name + "Highest> " + self.name + ";"
    self.model_code = ("  " + self.name + " ~ raised_cosine(" +
                       self.name + "Lowest, " +
                       self.name + "Low, " +
                       self.name + "High, " +
                       self.name + "Highest);")
    self.data_code = ("  real " + self.name + "Lowest;\n" +
                      "  real<lower=" + self.name + "Lowest> " + self.name + "Low;\n"
                      "  real<lower=" + self.name + "Low> " + self.name + "High;\n"
                      "  real<lower=" + self.name + "High> " + self.name + "Highest;")
    self.function_code = (
"""
  /**
   * Probability density function for raised cosine prior
   *
   * @param y     point to evaluate density on
   * @param start lowest possible value
   * @param u1    lowest expected value
   * @param u2    highest expected value
   * @param end   highest possible value
   * @return      probability density value.
   */
  real raised_cosine_lpdf(real y, real start, real u1, real u2, real end) {
    real res;
    real norm;
    if (y < u1) {
      res = log(0.5-0.5*cos(pi()/(u1-start)*(y-start)));
    } else if (y > u2) {
      res = log(0.5+0.5*cos(pi()/(end-u2)*(y-u2)));
    } else if (u1 <= y && y <= u2) {
      res = 0.0;
    } else {
      reject("input must be inside the specified range!");
    }
    norm = log(0.5*(u1-start+end-u2) + (u2-u1));
      
    return res-norm;
  }""")

  def register(self, bdsmodel):
    super().register(bdsmodel)

    bdsmodel._code_blocks['data'].append(self.data_code)
    bdsmodel._code_blocks['functions']['raised_cosine'] = self.function_code
