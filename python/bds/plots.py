import pandas as pd
import numpy as np
import plotnine as gg
import arviz

from dfply import *

def plot_deviance(scale):
  dev_df = (pd.DataFrame({'deviance': scale.deviance, 'origin': 'empirical'})
         >> bind_rows(pd.DataFrame({'deviance': scale.ppc_deviance, 'origin': 'simulated'}))
        )

  return (gg.ggplot(dev_df, gg.aes(x='deviance', fill='origin')) +
          gg.geom_histogram(binwidth=1.0))

def plot_ordered_residuals(scale):
  return (gg.ggplot(scale.ppc_ordered_residuals, gg.aes(x='residuals_median', y='sortid', color='origin')) +
          gg.geom_line() +
          gg.geom_line(gg.aes(x='residuals_p250')) + 
          gg.geom_line(gg.aes(x='residuals_p975')))

def plot_reversals(scale):
  return (gg.ggplot(scale.ppc_residual_reversals, gg.aes(x='runs', fill='origin')) +
       gg.geom_histogram(binwidth=2.0))

def plot_flip_count(scale):
  return (gg.ggplot(scale.ppc_flip_count, gg.aes(x='runs')) +
       gg.geom_histogram(binwidth=2.0) +
       gg.geom_vline(gg.aes(xintercept=scale.summary_statistics['flip count'])))

def plot_scale(scale):
  df = pd.DataFrame({'Stimulus': scale.stimulus,
                     'Scale': np.mean(scale.scale, axis=0),
                     'Low': np.percentile(scale.scale, 2.5, axis=0),
                     'High': np.percentile(scale.scale, 97.5, axis=0)})

  return (gg.ggplot(df, gg.aes(x='Stimulus', y='Scale', ymin='Low', ymax='High')) +
          gg.geom_point() +
          gg.geom_line() +
          gg.geom_ribbon(alpha=0.3));

def plot_posterior_samples(scale, samples=200):
  pp = scale.scale.T
  pp = pp[:, np.random.choice(pp.shape[1], samples, replace=False)]

  cols = ['smp[%d]' % x for x in range(0, pp.shape[1])]

  scales_df = pd.DataFrame(data = pp,
                           columns = cols)


  df = (pd.DataFrame({'Stimulus': scale.stimulus})
        >> bind_cols(scales_df)
        >> gather('smp', 'Scale', cols))

  return (gg.ggplot(df, gg.aes(x='Stimulus', y='Scale', grouping='smp')) +
          gg.geom_line(alpha=0.1))
