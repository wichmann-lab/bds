import numpy as np

import pystan

class DifferenceScale:
  def __init__(self, stimulus, stanfit):
    self.stanfit = stanfit
    self.stimulus = stimulus

    k = len(self.stimulus)
    par_names = ['psi[' + str(x) + ']' for x in range(1, k-1)]
#    par_names.append('precision')

    summary = stanfit.summary(pars=par_names, probs=[0.025, 0.25, 0.5, 0.75, 0.975])['summary']
    self.scale = np.concatenate([np.zeros((1,summary.shape[1])),
                                summary,
                                np.ones((1,summary.shape[1]))], axis=0)

    self.precision = stanfit.summary(pars=['precision'], probs=[0.025, 0.25, 0.5, 0.75, 0.975])['summary']

  def get_scale_values(self):
    return self.scale[:,0]

  def get_scale_credible_interval(self):
    return (self.scale[:,3], self.scale[:,7])

  def __repr__(self):
    return repr(self.get_scale_values(self))

  def __str__(self):
    return ("fitted BDS scale\nstimulus values:\n" + str(self.stimulus) +
           "\nscale values:\n" + str(self.scale[:,0]) +
           "\nprecision:\t" + str(self.precision[:,0]))

class LpsDifferenceScale(DifferenceScale):
  def __init__(self, stimulus, stanfit):
    super().__init__(stimulus, stanfit)

    self.lapserate = stanfit.summary(pars=['lapses'], probs=[0.025, 0.25, 0.5, 0.75, 0.975])['summary']

  def __str__(self):
    return super().__str__() + "\nlapse rate:\t" + str(self.lapserate[:,0])
