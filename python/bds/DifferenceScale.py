import numpy as np
import pandas as pd

import pystan

#try:
#  import .plots
#  plotting_enabled = True
#except ImportError:
#  plotting_enabled = False

class DifferenceScale:
  def __init__(self, stan_fit, stan_data, data):
    self.stan_fit = stan_fit
    self.stimulus = stan_data['Stimulus']
    self.n = stan_data['N']
    self.k = stan_data['K']

    self.stan_data = stan_data

    self.data = data
    self.pvalues = dict()
    self.summary_statistics = dict()
    self.diagnostic_plots = dict()
    self.plots = dict()

    self.lapserate = None
    self.scale = None
    self.sensitivity = None
    self.deviance = None

  def get_scale_values(self):
    return np.mean(self.scale, axis=0)

  def get_scale_credible_interval(self):
    return (np.percentile(self.scale, 2.5, axis=0),
            np.percentile(self.scale, 97.5, axis=0))

  def get_sensitivity(self):
    return np.mean(self.sensitivity)

  def get_sensitivity_credible_interval(self):
    return (np.percentile(self.sensitivity, 2.5), np.percentile(self.sensitivity, 97.5))

  def get_lapserate(self):
    if (self.lapserate):
      return np.mean(self.lapserate)
    else:
      raise ValueError('Apparently no lapse rate was fitted! lapserate is None.')

  def get_lapserate_credible_interval(self):
    if (self.lapserate):
      return (np.percentile(self.lapserate, 2.5), np.percentile(self.lapserate, 97.5))
    else:
      raise ValueError('Apparently no lapse rate was fitted! lapserate is None.')

  def get_deviance(self):
    return np.mean(self.deviance)

  def get_deviance_credible_interval(self):
    return (np.percentile(self.deviance, 2.5), np.percentile(self.deviance, 97.5))

  def __repr__(self):
    return repr(self.get_scale_values(self))

  def __str__(self):
    return ("fitted BDS scale\nstimulus values:\n" + str(self.stimulus) +
           "\nscale values:\n" + str(self.scale[:,0]) +
           "\nsensitivity:\t" + str(self.sensitivity[:,0]))
