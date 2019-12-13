from .RegressorModel import RegressorModel

class SensoryNoiseModel(RegressorModel):
  def __init__(self, link, sensitivity_prior):
    super().__init__(link)

    self.sns_prior = sensitivity_prior

    self.formula = ("X * psi[2:K] * sensitivity", "Y * psi[2:K] * sensitivity")

    self.build_code()

    self.data_transform_def = "  matrix[N, K-1] X;\n  matrix[N, K-1] Y;"
    self.data_transform = "  X = construct_design_matrix(S1, S2, S3, S4, K);\n  Y = construct_design_matrix2(S1, S2, S3, S4, K);"

    self.functions = (
"""/**
 * Rearrange the data to fit the framework of a GLM.
 *
 * @param N Number of rows corresponding to data items
 * @param K Number of predictors, counting the intercept, per
 *          item.
 * @return Simulated predictor matrix.
 */
  matrix construct_design_matrix(int[] L1, int[] L2, int[] R1, int[] R2, int K) {
    matrix[size(L1), K-1] X;
    
    X = rep_matrix(0.0, size(L1), K-1);

    for (n in 1:size(L1)) {
      int a = L1[n];
      int b = L2[n];
      int c = R1[n];
      int d = R2[n];
    
      row_vector[K] rvec;
      rvec = rep_row_vector(0.0, K);
    
      rvec[a] += 1.0;
      rvec[b] -= 1.0;
      rvec[c] -= 1.0;
      rvec[d] += 1.0;
      X[n] = rvec[2:K];
    }
    
    return X;
  }

/**
 * Rearrange the data to fit the framework of a GLM.
 *
 * @param N Number of rows corresponding to data items
 * @param K Number of predictors, counting the intercept, per
 *          item.
 * @return Simulated predictor matrix.
 */
  matrix construct_design_matrix2(int[] L1, int[] L2, int[] R1, int[] R2, int K) {
    matrix[size(L1), K-1] X;
    
    X = rep_matrix(0.0, size(L1), K-1);

    for (n in 1:size(L1)) {
      int a = L1[n];
      int b = L2[n];
      int c = R1[n];
      int d = R2[n];
    
      row_vector[K] rvec;
      rvec = rep_row_vector(0.0, K);
    
      rvec[a] -= 1.0;
      rvec[d] += 1.0;
      X[n] = rvec[2:K];
    }
    
    return X;
  }
""")

  def register(self, bdsmodel):
    self.sns_prior.register(bdsmodel)
    super().register(bdsmodel)

    bdsmodel._code_blocks['functions']['construct_design_matrix'] = self.functions
    bdsmodel._code_blocks['transformed data:definitions'].append(self.data_transform_def)
    bdsmodel._code_blocks['transformed data'].append(self.data_transform)

  def results(self, result_obj):
    result_obj.sensitivity = self.sensitivity(result_obj.stan_fit)

  def sensitivity(self, stan_fit):
    return stan_fit.extract(pars=self.sns_prior.name)[self.sns_prior.name]
