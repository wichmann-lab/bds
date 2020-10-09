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
    self.responses = stan_data['Response']

    self.stan_data = stan_data

    self.data = data
    self.pvalues = dict()

    self.lapserate = None
    self.scale = None
    self.sensitivity = None
    self.deviance = None
    self.loo = None

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
    if self.lapserate is not None:
      return np.mean(self.lapserate)
    else:
      raise ValueError('Apparently no lapse rate was fitted! lapserate is None.')

  def get_lapserate_credible_interval(self):
    if self.lapserate is not None:
      return (np.percentile(self.lapserate, 2.5), np.percentile(self.lapserate, 97.5))
    else:
      raise ValueError('Apparently no lapse rate was fitted! lapserate is None.')

  def get_deviance(self):
    return np.mean(self.deviance)

  def get_deviance_credible_interval(self):
    return (np.percentile(self.deviance, 2.5), np.percentile(self.deviance, 97.5))

  def summary(self):
    summary = pd.DataFrame()

    summary = summary.append(pd.DataFrame({'type': 'scale',
                                           'name': [str(i) for i in range(1, self.k+1)],
                                           'stimulus': self.stimulus,
                                           'value': self.get_scale_values(),
                                           'low': self.get_scale_credible_interval()[0],
                                           'high': self.get_scale_credible_interval()[1]}),
                             ignore_index=True)

    summary = summary.append(pd.DataFrame({'type': 'unconstrained',
                                           'name': [str(i) for i in range(1, self.k+1)],
                                           'stimulus': self.stimulus,
                                           'value': np.mean(self.scale * self.sensitivity[:,None], axis=0),
                                           'low': np.percentile(self.scale * self.sensitivity[:,None], 2.5, axis=0),
                                           'high': np.percentile(self.scale * self.sensitivity[:,None], 97.5, axis=0)}),
                             ignore_index=True)

    if self.lapserate is not None:
      summary = summary.append(pd.DataFrame({'type': 'lapses',
                                             'name': 'lapses',
                                             'stimulus': -1,
                                             'value': self.get_lapserate(),
                                             'low': self.get_lapserate_credible_interval()[0],
                                             'high': self.get_lapserate_credible_interval()[1]}, index=[0]),
                               ignore_index=True)

    summary = summary.append(pd.DataFrame({'type': 'sensitivity',
                                           'name': 'sensitivity',
                                           'stimulus': -1,
                                           'value': self.get_sensitivity(),
                                           'low': self.get_sensitivity_credible_interval()[0],
                                           'high': self.get_sensitivity_credible_interval()[1]}, index=[0]),
                             ignore_index=True)

    summary = summary.append(pd.DataFrame({'type': 'gof',
                                           'name': [*self.pvalues.keys()],
                                           'stimulus': -1,
                                           'value': [*self.pvalues.values()],
                                           'low': -1,
                                           'high': -1}),
                             ignore_index=True)

    summary = summary.append(pd.DataFrame({'type': 'gof',
                                           'name': 'deviance',
                                           'stimulus': -1,
                                           'value': self.get_deviance(),
                                           'low': self.get_deviance_credible_interval()[0],
                                           'high': self.get_deviance_credible_interval()[1]}, index=[0]),
                             ignore_index=True)

    if self.loo is not None:
      summary = summary.append(pd.DataFrame({'type': 'gof',
                                             'name': 'loo',
                                             'stimulus': -1,
                                             'value': self.loo['IC'],
                                             'low': self.loo['IC'] - self.loo['SE'],
                                             'high': self.loo['IC'] + self.loo['SE']}, index=[0]),
                               ignore_index=True)

      summary = summary.append(pd.DataFrame({'type': 'gof',
                                             'name': 'eff_param',
                                             'stimulus': -1,
                                             'value': self.loo['pIC'],
                                             'low': -1,
                                             'high': -1}, index=[0]),
                               ignore_index=True)

      summary = summary.append(pd.DataFrame({'type': 'gof',
                                             'name': 'loo_warning',
                                             'stimulus': -1,
                                             'value': self.loo['warning'],
                                             'low': -1,
                                             'high': -1}, index=[0]),
                               ignore_index=True)

    return summary

  def __repr__(self):
    return repr(self.get_scale_values(self))

  def __str__(self):
    return ("fitted BDS scale\nstimulus values:\n" + str(self.stimulus) +
           "\nscale values:\n" + str(self.scale[:,0]) +
           "\nsensitivity:\t" + str(self.sensitivity[:,0]))
