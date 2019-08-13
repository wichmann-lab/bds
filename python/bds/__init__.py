import numpy as np
import os

from .DifferenceScale import *
from .BDSModel import BDSModel
from .likelihood import *
from .regressor import *
from .prior import *
from .scale import *

default_model = BDSModel('bds',
                         likelihood = BinomialMixture(lapses_prior=BetaDistribution('lapses', {'lapsesAlpha': 1, 'lapsesBeta': 10})),
                         regressor = LinearModel(link=ProbitLink(),
                                                 sensitivity_prior=RaisedCosineDistribution('sensitivity', {'sensitivityLowest': 0, 'sensitivityLow': 2.5, 'sensitivityHigh': 25, 'sensitivityHighest': 50})),
                         scale = MonotonicScale(diff_prior=DirichletDistribution('psi_diff', 'K-1', {'psi_diffAlpha': 1})))

diff_model = BDSModel('bds_diff',
                         likelihood = BinomialMixture(lapses_prior=BetaDistribution('lapses', {'lapsesAlpha': 1, 'lapsesBeta': 10})),
                         regressor = DifferenceModel(link=ProbitLink(),
                                                 sensitivity_prior=RaisedCosineDistribution('sensitivity', {'sensitivityLowest': 0, 'sensitivityLow': 2.5, 'sensitivityHigh': 25, 'sensitivityHighest': 50})),
                         scale = ParameterScale(scale_prior=VectorizedUniformDistribution('psi_hat', 'K-2', {'psi_hatStart': 0.0, 'psi_hatEnd': 1.0})))

gp_model = BDSModel('bds_gp',
                    likelihood = BinomialMixture(lapses_prior=BetaDistribution('lapses', {'lapsesAlpha': 1, 'lapsesBeta': 10})),
                    regressor = DifferenceModel(link=ProbitLink(),
                                                 sensitivity_prior=RaisedCosineDistribution('sensitivity', {'sensitivityLowest': 0, 'sensitivityLow': 2.5, 'sensitivityHigh': 25, 'sensitivityHighest': 50})),
                    scale = GaussianProcess(length_scale_prior = HalfNormalDistribution('rho', {'rhoSigma': 2}),
                                            magnitude_prior = HalfNormalDistribution('alpha', {'alphaSigma': 2})))

def bds(data, stimulus=None, **kwargs):

  result = default_model.sample(data, stimulus, **kwargs)

#  gp_dict = {'N_predict': 100,
#             'x_predict': np.concatenate([stimulus[:-1], np.linspace(stimulus[0], stimulus[-1], num=100-stimulus.shape[0]+2)[1:]]),
#             'obs_idx': np.array(range(2, stimulus.shape[0])),
#             'rho': 0.5*(stimulus[-1]-stimulus[0])}

#  result = gp_model.sample(data, stimulus, params=gp_dict)
  return result
