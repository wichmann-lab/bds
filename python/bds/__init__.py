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

def bds(data, stimulus=None, **kwargs):

  if not 'init' in kwargs.keys():
    k = int(max(data.loc[:, ('S1', 'S2', 'S3')].max()))
    init = lambda: {'psi_diff': np.ones(k-1) / (k-1),
                    'lapses': 0.01,
                    'sensitivity': default_model.hyper_params['sensitivityLow']}
    
    kwargs['init'] = init

  result = default_model.sample(data, stimulus, **kwargs)

  return result
