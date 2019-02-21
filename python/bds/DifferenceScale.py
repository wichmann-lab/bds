import numpy as np
import pandas as pd

from dfply import *
import pystan

import plotnine as gg

class DifferenceScale:
  def __init__(self, stimulus, stanfit, data):
    self.stanfit = stanfit
    self.stimulus = stimulus
    if (data.shape[1] == 4):
      colnames = ['S1', 'S2', 'S3', 'Response']
    else:
      colnames = ['S1', 'S2', 'S3', 'S4', 'Response']

    self.data = pd.DataFrame(data, columns=colnames)
    self.n = data.shape[0]

    self.data.loc[:, 'trial'] = np.arange(1, self.n+1)

    k = len(self.stimulus)
    par_names = ['psi[' + str(x) + ']' for x in range(1, k-1)]
#    par_names.append('precision')

    summary = stanfit.summary(pars=par_names, probs=[0.025, 0.25, 0.5, 0.75, 0.975])['summary']
    self.scale = np.concatenate([np.zeros((1,summary.shape[1])),
                                summary,
                                np.ones((1,summary.shape[1]))], axis=0)

    self.precision = stanfit.summary(pars=['precision'], probs=[0.025, 0.25, 0.5, 0.75, 0.975])['summary']

  def ppc(self):
    rhat_pars = ['resp_hat[%d]' % x for x in range(1,self.n+1)]
    ll_pars =   ['log_lik[%d]' % x for x in range(1,self.n+1)]
    pp = (self.stanfit.to_dataframe(pars=rhat_pars,
                                   inc_warmup=False,
                                   diagnostics=False)
          >> select(rhat_pars)
         ).T
    ll = (self.stanfit.to_dataframe(pars=ll_pars, 
                                  inc_warmup=False,
                                  diagnostics=False)
          >> select(ll_pars)
         ).T

    ppcols = ['smp[%d]' % x for x in range(0, pp.shape[1])]
    llcols = ['smp[%d]' % x for x in range(0, pp.shape[1])]
    pp.columns = ppcols
    ll.columns = llcols

    pp.loc[:,'trial'] = np.arange(1, self.n+1)
#    ll.loc[:,'trial'] = np.arange(1, self.n+1)

    ll >>= gather('smp', 'll')

    ppc = self.data.set_index('trial') >> bind_cols(pp.set_index('trial'))

    ppc = ppc.reset_index()
    ppc >>= gather('smp', 'rhat', ppcols)

    ppc.loc[:,'ll'] = ll.loc[:,'ll']

    # make np.sqrt usable with dfply's placeholder variables
    @make_symbolic
    def np_sqrt(a):
      return np.sqrt(a)
    ppc >>= mutate(dev_res=(2.0*X.Response-1.0) * np_sqrt(-2.0 * X.ll),
                   dev_res_pred=(2.0*X.rhat-1.0) * np_sqrt(-2.0 * X.ll))
         
    def p250(x):
      return np.percentile(x,2.5)

    def p975(x):
      return np.percentile(x, 97.5)

    ppc_summary = (ppc >> group_by('trial')
                       >> summarize_each([p250, p975, np.median],
                                         X.dev_res, X.dev_res_pred)
                       >> select(X.dev_res_median, X.dev_res_pred_median,
                                 X.dev_res_p250, X.dev_res_pred_p250,
                                 X.dev_res_p975, X.dev_res_pred_p975)
                  )

    ppc_summary >>= arrange(X.dev_res_median)

    ppc_summary = ppc_summary.apply(np.sort)

    # create a column with the sorted order
    ppc_summary.reset_index(drop=True, inplace=True)
    ppc_summary['sorted'] = ppc_summary.index

    def count_revs(a):
      return np.sum(np.abs(np.diff(a)))

    reversals = (ppc >> arrange(X.S1, X.S2, X.S3)
                     >> group_by(X.smp)
                     >> summarize_each([count_revs],X.rhat)
                )
    sorted_resp = self.data >> arrange(X.S1, X.S2, X.S3)

    print(sorted_resp)

    emp_rev = count_revs(sorted_resp['Response'])
    return (ppc_summary, reversals, emp_rev)

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

  def ppc_plot(self):
    (ppc_df, rev_df, emp_rev) = self.ppc()

    print(rev_df)
    print('empirical reversal', emp_rev)
    return (gg.ggplot(ppc_df, gg.aes(x='dev_res_median', y='sorted')) +
       gg.geom_line(color='red') +
       gg.geom_line(gg.aes(x='dev_res_p250'), linetype='dashed', color='red') + 
       gg.geom_line(gg.aes(x='dev_res_p975'), linetype='dashed', color='red') +
       gg.geom_line(gg.aes(x='dev_res_pred_p250'), linetype='dashed', color='black') +
       gg.geom_line(gg.aes(x='dev_res_pred_p975'), linetype='dashed', color='black') +
       gg.geom_line(gg.aes(x='dev_res_pred_median'), linetype='solid', color='orange'),
       gg.ggplot(rev_df, gg.aes(x='rhat_count_revs')) +
       gg.geom_histogram() + gg.geom_vline(gg.aes(xintercept=emp_rev))
    )

class LpsDifferenceScale(DifferenceScale):
  def __init__(self, stimulus, stanfit, trials):
    super().__init__(stimulus, stanfit, trials)

    self.lapserate = stanfit.summary(pars=['lapses'], probs=[0.025, 0.25, 0.5, 0.75, 0.975])['summary']

  def __str__(self):
    return super().__str__() + "\nlapse rate:\t" + str(self.lapserate[:,0])
