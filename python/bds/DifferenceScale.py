import numpy as np
import pandas as pd

from dfply import *
import pystan

import plotnine as gg

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

  def plot(self):
    df = pd.DataFrame({'Stimulus': self.stimulus, 'Scale': np.mean(self.scale, axis=0), 'Low': np.percentile(self.scale, 2.5, axis=0), 'High': np.percentile(self.scale, 97.5, axis=0)})

    return (gg.ggplot(df, gg.aes(x='Stimulus', y='Scale', ymin='Low', ymax='High')) +
            gg.geom_point() +
            gg.geom_line() +
            gg.geom_ribbon(alpha=0.3));

  def plot_posterior_samples(self, samples=200):
    pp = self.scale.T
    pp = pp[:, np.random.choice(pp.shape[1], samples, replace=False)]

    cols = ['smp[%d]' % x for x in range(0, pp.shape[1])]

    scales_df = pd.DataFrame(data = pp,
                           columns = cols)


    df = (pd.DataFrame({'Stimulus': self.stimulus})
          >> bind_cols(scales_df)
          >> gather('smp', 'Scale', cols))

    return (gg.ggplot(df, gg.aes(x='Stimulus', y='Scale', grouping='smp')) +
            gg.geom_line(alpha=0.1))
    
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

  def __repr__(self):
    return repr(self.get_scale_values(self))

  def __str__(self):
    return ("fitted BDS scale\nstimulus values:\n" + str(self.stimulus) +
           "\nscale values:\n" + str(self.scale[:,0]) +
           "\nsensitivity:\t" + str(self.sensitivity[:,0]))
