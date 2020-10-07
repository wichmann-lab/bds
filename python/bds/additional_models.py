from .DifferenceScale import *
from .BDSModel import BDSModel
from .BESSModel import BESSModel
from .likelihood import *
from .regressor import *
from .prior import *
from .scale import *

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

bds_model_constant_lapses = BDSModel('bds_constant_lapses',
                         likelihood = BinomialMixture(lapses_prior=Constant('lapses', 'real<lower=0, upper=1>')),
                         regressor = LinearModel(link=ProbitLink(),
                                                 sensitivity_prior=RaisedCosineDistribution('sensitivity', {'sensitivityLowest': 0, 'sensitivityLow': 2.5, 'sensitivityHigh': 25, 'sensitivityHighest': 50})),
                         scale = MonotonicScale(diff_prior=DirichletDistribution('psi_diff', 'K-1', {'psi_diffAlpha': 1})))

sensory_noise_model = BDSModel('bds_sensory_noise',
                         likelihood = BinomialMixture(lapses_prior=BetaDistribution('lapses', {'lapsesAlpha': 1, 'lapsesBeta': 10})),
                         regressor = SensoryNoiseModel(link=SensoryNoiseProbitLink(),
                                                 sensitivity_prior=RaisedCosineDistribution('sensitivity', {'sensitivityLowest': 0, 'sensitivityLow': 2.5, 'sensitivityHigh': 25, 'sensitivityHighest': 50})),
                         scale = MonotonicScale(diff_prior=DirichletDistribution('psi_diff', 'K-1', {'psi_diffAlpha': 1})))

abs_model = BDSModel('bds_sensory_noise',
                         likelihood = BinomialMixture(lapses_prior=BetaDistribution('lapses', {'lapsesAlpha': 1, 'lapsesBeta': 10})),
                         regressor = DifferenceModel(link=ProbitLink(),
                                                 sensitivity_prior=RaisedCosineDistribution('sensitivity', {'sensitivityLowest': 0, 'sensitivityLow': 2.5, 'sensitivityHigh': 25, 'sensitivityHighest': 50})),
                         scale = MonotonicScale(diff_prior=DirichletDistribution('psi_diff', 'K-1', {'psi_diffAlpha': 1})))

def bds_gp(data, stimulus, predictive, **kwargs):
  gp_dict = {'N_predict': len(predictive),
             'x_predict': predictive}

  init = lambda: {'psi_tilde': np.zeros(len(stimulus)+len(predictive)-2),
                  'lapses': 0.01,
                  'alpha': 0.5,
                  'rho': 0.5,
                  'sensitivity': gp_model.hyper_params['sensitivityLow']}
  
  result = gp_model.sample(data, stimulus, params=gp_dict, init=init, **kwargs)
  return result

def bess(data, predictive, **kwargs):
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

  gp_dict = {'N_predict': len(predictive),
             'x_predict': predictive}

#  init = lambda: {'psi_tilde': np.zeros((len(stimulus)+len(predictive)-2)),
#                  'alpha': 0.1,
#                  'rho': 1,
#                  'sensitivity': 5}
  
  result = bess_gp_model.sample(new_data, stimulus, params=gp_dict, **kwargs)
  return result

def bds_constant_lapserate(data, stimulus=None, lapserate=0.0, **kwargs):

  if not 'init' in kwargs.keys():
    k = int(max(data.max()))

    init = lambda: {'psi_tilde': np.ones(k - 1) / k,
                    'sensitivity': bds_model_constant_lapses.hyper_params['sensitivityLow']}
                                                      
    kwargs['init'] = init
                                                      
  result = bds_model_constant_lapses.sample(data, stimulus, params={'lapses': lapserate}, **kwargs)

  return result

def bds_sensory_noise(data, stimulus=None, **kwargs):

  if not 'init' in kwargs.keys():
    k = int(max(data.max()))

    init = lambda: {'psi_tilde': np.ones(k - 1) / k,
                    'lapses': 0.01,
                    'sensitivity': sensory_noise_model.hyper_params['sensitivityLow']}
                                                      
    kwargs['init'] = init
                                                      
  result = sensory_noise_model.sample(data, stimulus, **kwargs)

  return result

def bds_abs(data, stimulus=None, **kwargs):

  if not 'init' in kwargs.keys():
    k = int(max(data.max()))
            
    init = lambda: {'psi_tilde': np.ones(k - 1) / k,
                    'lapses': 0.01,
                    'sensitivity': abs_model.hyper_params['sensitivityLow']}
                                                      
    kwargs['init'] = init
                                                      
  result = abs_model.sample(data, stimulus, **kwargs)

  return result
