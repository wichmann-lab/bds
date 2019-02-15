import numpy as np
import os

from bds.DifferenceScale import *
from bds.stan_helpers import fit_stan_model

def bds(mldsdata,
        stimulus=None,
        lapses=True,
        precLowest=0,
        precLow=2,
        precHigh=20,
        precHighest=30,
        lpsAlpha=1,
        lpsBeta=5):

  if mldsdata.shape[1] != 4 and mlds.data.shape[1] != 5:
    raise ValueError('MLDS data should have 4 or 5 columns!')

  modelname = 'bds'
  data = dict()

  data['S1'] = mldsdata[:,0]

  # Do we have triad or quadtruple data?
  if mldsdata.shape[1] == 4:
    data['S2'] = mldsdata[:,1]
    data['S3'] = mldsdata[:,1]
    data['S4'] = mldsdata[:,2]
    data['Responses'] = mldsdata[:,3]
    data['K'] = np.amax(mldsdata[:,0:2])
  else:
    data['S2'] = mldsdata[:,1]
    data['S3'] = mldsdata[:,2]
    data['S4'] = mldsdata[:,3]
    data['Responses'] = mldsdata[:,4]
    data['K'] = np.amax(mldsdata[:,0:3])

  data['N'] = mldsdata.shape[0]

  if stimulus is None:
    stimulus = np.linspace(0.0, 1.0, num=data['K'])

  # Set hyperparameters for prior on the decision noise
  data['precLowest'] = precLowest
  data['precLow'] = precLow
  data['precHigh'] = precHigh
  data['precHighest'] = precHighest

  # find installed model files
  modeldir = os.path.dirname(__file__).split('lib/python')[0] + 'share/models/'

#  print(modeldir)

  # Should a lapse rate be fitted?
  if lapses:
    data['lpsAlpha'] = lpsAlpha
    data['lpsBeta'] = lpsBeta

    stanfit = fit_stan_model(modeldir, 'bds_lps', data)
    result = LpsDifferenceScale(stimulus, stanfit)
  else:
    stanfit = fit_stan_model(modeldir, 'bds', data)
    result = DifferenceScale(stimulus, stanfit)

  return result
