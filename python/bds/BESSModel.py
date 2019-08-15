import pandas as pd
import numpy as np

from .DifferenceScale import DifferenceScale
from .stan_helpers import *

class BESSModel:
  def __init__(self, model_name, likelihood, regressor, scale, priors=dict()):
    self.unregister()

    self.hyper_params = dict()

    self.model_name = model_name

    self.likelihood = likelihood
    self.regressor = regressor
    self.scale = scale

    self.register_modules()
    self.build_model()

  def unregister(self):
    self._code_blocks = {'functions': dict(),
                         'data': list(),
                         'transformed data:definitions': list(),
                         'transformed data': list(),
                         'parameters': list(),
                         'transformed parameters:definitions': list(),
                         'transformed parameters': list(),
                         'model': list(),
                         'generated quantities:definitions': list(),
                         'generated quantities': list()
                        }

    self._code_blocks['data'].append(
"""  int<lower=1> N;  // total number of observations 
  int<lower=1> K;  // number of stimulus levels
  int<lower=1, upper=K> L[N];
  int<lower=1, upper=K> U[N];
  int<lower=1, upper=K> Response[N];  // response variable
  real Stimulus[K];
""")

  def register_modules(self):
    self.scale.register(self)
    self.regressor.register(self)
    self.likelihood.register(self)

  def build_model(self):
    self.model_code = (
      'functions {\n' +
      '\n'.join(self._code_blocks['functions'].values()) +
      '\n}\n\ndata {\n' +
      '\n'.join(self._code_blocks['data']) +
      '\n}\n\ntransformed data {\n' +
      '\n'.join(self._code_blocks['transformed data:definitions']) + '\n\n' +
      '\n'.join(self._code_blocks['transformed data']) +
      '\n}\n\nparameters {\n' +
      '\n'.join(self._code_blocks['parameters']) +
      '\n}\n\ntransformed parameters {\n' +
      '\n'.join(self._code_blocks['transformed parameters:definitions']) + '\n\n' +
      '\n'.join(self._code_blocks['transformed parameters']) +
      '\n}\n\nmodel {\n' +
      '\n'.join(self._code_blocks['model']) +
      '\n}\n\ngenerated quantities {\n' +
      '\n'.join(self._code_blocks['generated quantities:definitions']) + '\n\n' +
      '\n'.join(self._code_blocks['generated quantities']) +
      '\n}\n')

    self.stan_model = StanModel_cache(self.model_code, model_name=self.model_name)

  def sample(self, data, stimulus, params=dict(), **kwargs):
    if not 'init' in kwargs.keys():
      kwargs['init'] = 'random'

    if not 'control' in kwargs.keys():
      kwargs['control'] = {'adapt_delta': 0.99}

    if not 'chains' in kwargs.keys():
      kwargs['chains'] = 4

    if not 'iter' in kwargs.keys():
      kwargs['iter'] = 2000

    if not 'warmup' in kwargs.keys():
      kwargs['warmup'] = 1000

    if not 'pars' in kwargs.keys():
      pass

    stan_data = self.hyper_params.copy()

    for k, v in params.items():
      stan_data[k] = v

    stan_data['L'] = data['L']
    stan_data['U'] = data['U']
    stan_data['Response'] = data['Response']
    stan_data['N'] = data.shape[0]
    stan_data['K'] = stimulus.shape[0]

    stan_data['Stimulus'] = stimulus

    fit = self.stan_model.sampling(data=stan_data,
                       iter=kwargs['iter'],
                       warmup=kwargs['warmup'],
                       chains=kwargs['chains'],
                       control=kwargs['control'],
#                       pars=kwargs['pars'],
                       init=kwargs['init'])

    scale_result = DifferenceScale(fit, stan_data, data)

    self.likelihood.results(scale_result)
    self.regressor.results(scale_result)
    self.scale.results(scale_result)

    return scale_result
