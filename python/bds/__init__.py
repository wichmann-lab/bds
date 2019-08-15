import numpy as np
import os

from .DifferenceScale import *
from .BDSModel import BDSModel
from .BESSModel import BESSModel
from .likelihood import *
from .regressor import *
from .prior import *
from .scale import *

default_model = BDSModel('bds',
                         likelihood = BinomialMixture(lapses_prior=BetaDistribution('lapses', {'lapsesAlpha': 1, 'lapsesBeta': 10})),
                         regressor = LinearModel(link=ProbitLink(),
                                                 sensitivity_prior=RaisedCosineDistribution('sensitivity', {'sensitivityLowest': 0, 'sensitivityLow': 2.5, 'sensitivityHigh': 25, 'sensitivityHighest': 50})),
                         scale = MonotonicScale(diff_prior=DirichletDistribution('psi_diff', 'K-1', {'psi_diffAlpha': 1})))

bess_gp_model = BESSModel('bess',
                         likelihood = NormalLikelihood(),
                         regressor = BESSDifferenceModel(sensitivity_prior=RaisedCosineDistribution('sensitivity', {'sensitivityLowest': 0, 'sensitivityLow': 2.5, 'sensitivityHigh': 25, 'sensitivityHighest': 50})),
                         scale = GaussianProcess(length_scale_prior = HalfNormalDistribution('rho', {'rhoSigma': 2}),
                                            magnitude_prior = HalfNormalDistribution('alpha', {'alphaSigma': 2})))

gp_model = BDSModel('bds_gp',
                    likelihood = BinomialMixture(lapses_prior=BetaDistribution('lapses', {'lapsesAlpha': 1, 'lapsesBeta': 10})),
                    regressor = DifferenceModel(link=ProbitLink(),
                                                 sensitivity_prior=RaisedCosineDistribution('sensitivity', {'sensitivityLowest': 0, 'sensitivityLow': 2.5, 'sensitivityHigh': 25, 'sensitivityHighest': 50})),
                    scale = GaussianProcess(length_scale_prior = HalfNormalDistribution('rho', {'rhoSigma': 2}),
                                            magnitude_prior = HalfNormalDistribution('alpha', {'alphaSigma': 2})))

def bds(data, stimulus=None, **kwargs):

  result = default_model.sample(data, stimulus, **kwargs)

  return result

def gp_bds(data, stimulus, **kwargs):
  gp_dict = {'N_predict': 100,
             'x_predict': np.linspace(stimulus[0], stimulus[-1], num=100)}

  result = gp_model.sample(data, stimulus, params=gp_dict, **kwargs)
  return result

def gp_bess(data, **kwargs):
  new_data = pd.DataFrame()

  lower = data['lower'].values
  upper = data['upper'].values
  resp = data['value'].values

  stimulus = np.unique(np.concatenate([lower, resp, upper]))
  stimulus.sort()

  l = np.zeros(lower.shape, dtype=int)
  u = np.zeros(upper.shape, dtype=int)
  r = np.zeros(resp.shape, dtype=int)

  for i in range(0, stimulus.shape[0]):
    l[lower == stimulus[i]] = i+1
    u[upper == stimulus[i]] = i+1
    r[resp == stimulus[i]] = i+1

  new_data['Response'] = r
  new_data['L'] = l
  new_data['U'] = u

  gp_dict = {'N_predict': 1,
             'x_predict': np.array([0.5])}

  result = bess_gp_model.sample(new_data, stimulus, params=gp_dict, **kwargs)
  return result
