import pandas as pd
import numpy as np

from .DifferenceScale import DifferenceScale
from .stan_helpers import *

def order_data(data):
  df = pd.DataFrame()
  df['Response'] = data['Response']
  df['S1'] = data['S1']
  df['S2'] = data['S2']
  df['S3'] = data['S3']
  
  if ('S4' not in data.columns):
    df['S3'] = data['S2']
    df['S4'] = data['S3']
    
    for i, rows in data.iterrows():
      
      if data.at[i,'S1'] > data.at[i, 'S3']:
        df.at[i, 'S1'] = data.at[i, 'S3']
        df.at[i, 'S4'] = data.at[i, 'S1']

        df.at[i, 'Response'] = 1 - data.at[i, 'Response']
  else:
    for i, rows in data.iterrows():
      if data.at[i, 'S1'] > data.at[i, 'S2']:
        df.at[i, 'S1'] = data.at[i, 'S2']
        df.at[i, 'S2'] = data.at[i, 'S1']

      if data.at[i, 'S3'] > data.at[i, 'S4']:
        df.at[i, 'S3'] = data.at[i, 'S4']
        df.at[i, 'S4'] = data.at[i, 'S3']

      if data.at[i,'S1'] > data.at[i, 'S3']:
        df.at[i, 'S1'] = data.at[i, 'S3']
        df.at[i, 'S3'] = data.at[i, 'S1']

        df.at[i, 'S2'] = data.at[i, 'S4']
        df.at[i, 'S4'] = data[i, 'S2']

        df.at[i, 'Response'] = 1 - data.at[i, 'Response']

  return df

class BDSModel:
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
  int<lower=1, upper=K> S1[N];
  int<lower=1, upper=K> S2[N];
  int<lower=1, upper=K> S3[N];
  int<lower=1, upper=K> S4[N];
  int<lower=0, upper=1> Response[N];  // response variable
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

  def sample(self, data, stimulus=None, params=dict(), **kwargs):
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

    if not 'n_jobs' in kwargs.keys():
      kwargs['n_jobs'] = -1
    
    stan_data = self.hyper_params.copy()

    for k, v in params.items():
      stan_data[k] = v

    ordered_data = order_data(data)
    stan_data['S1'] = ordered_data['S1']
    stan_data['S2'] = ordered_data['S2']
    stan_data['S3'] = ordered_data['S3']
    stan_data['S4'] = ordered_data['S4']
    stan_data['Response'] = ordered_data['Response']
    stan_data['K'] = int(max(ordered_data.max()))

    stan_data['N'] = data.shape[0]

    if stimulus is None:
      stan_data['Stimulus'] = np.linspace(0.0, 1.0, num=stan_data['K'])
    else:
      stan_data['Stimulus'] = stimulus

    fit = self.stan_model.sampling(data=stan_data,
                       iter=kwargs['iter'],
                       warmup=kwargs['warmup'],
                       chains=kwargs['chains'],
                       control=kwargs['control'],
#                       pars=kwargs['pars'],
                       init=kwargs['init'],
                       n_jobs=kwargs['n_jobs'],
                       check_hmc_diagnostics=False)

    scale_result = DifferenceScale(fit, stan_data, ordered_data)

    self.likelihood.results(scale_result)
    self.regressor.results(scale_result)
    self.scale.results(scale_result)

    return scale_result
