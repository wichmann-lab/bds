import pandas as pd
import numpy as np
import plotnine as gg
import arviz as az

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
          gg.geom_pointrange() +
          gg.geom_line());

def plot_predictive_scale(scale):
  df = pd.DataFrame({'Stimulus': scale.pred_stimulus,
                     'Scale': np.mean(scale.pred_scale, axis=0),
                     'Low': np.percentile(scale.pred_scale, 2.5, axis=0),
                     'High': np.percentile(scale.pred_scale, 97.5, axis=0)})

  point_df = pd.DataFrame({'Stimulus': scale.stimulus,
                     'Scale': np.mean(scale.scale, axis=0),
                     'Low': np.percentile(scale.scale, 2.5, axis=0),
                     'High': np.percentile(scale.scale, 97.5, axis=0)})

  return (gg.ggplot(df, gg.aes(x='Stimulus', y='Scale', ymin='Low', ymax='High')) +
          gg.geom_point(data=point_df) +
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

from scipy.stats import norm

def plot_residuals(scale):
  pr = np.mean(scale.decision_probabilities, axis=0)
  lps = scale.get_lapserate()

  df = pd.DataFrame({'Delta(a,b,c)': norm.ppf(pr), 'probability': pr, 'response': scale.responses})

  return (gg.ggplot(df, gg.aes(x='Delta(a,b,c)', y='probability')) +
          gg.stat_function(fun= norm.cdf) +
          gg.geom_point(gg.aes(y='response'), size=3))

def compare_models_deviance(scales):
  dev_data = pd.DataFrame()

  for n, sc in scales.items():
    dev_mean = sc.get_deviance()
    dev_low, dev_high = sc.get_deviance_credible_interval()

    dev_data = dev_data.append({'deviance': dev_mean, 'dev_low': dev_low, 'dev_high': dev_high, 'model': n}, ignore_index=True) 

  return (gg.ggplot(dev_data, gg.aes(x='model', y='deviance', ymin='dev_low', ymax='dev_high')) +
          gg.geom_point() +
          gg.geom_errorbar())

def to_arviz(scale):
  return az.from_pystan(
    posterior=scale.stan_fit,
    posterior_predictive='resp_hat',
    observed_data=['Response'],
    log_likelihood={'Response': 'log_lik'},
    coords={'stimulus': np.arange(scale.stan_data['K']),
            'trial': np.arange(scale.stan_data['N'])},
    dims={'psi': ['stimulus'],
          'Response': ['trial'],
          'S1': ['trial'],
          'S2': ['trial'],
          'S3': ['trial'],
          'S4': ['trial'],
          'log_lik': ['trial'],
          'resp_hat': ['trial']}
    )

def compare_models_loo(scales):
  loo_data = pd.DataFrame()

  for n, sc in scales.items():
    az_data = to_arviz(sc)

    loo = az.loo(az_data, scale='deviance')
    loo['model'] = n
    loo_data = loo_data.append(loo,ignore_index=True)

  return (gg.ggplot(loo_data, gg.aes(x='model', y='loo', ymin='loo-loo_se', ymax='loo+loo_se')) +
          gg.geom_point() +
          gg.geom_text(gg.aes(label='p_loo'), nudge_x=0.2, format_string='{:.2}') +
          gg.geom_errorbar())
