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

  def ppc_ordered_residuals(self):
    rhat_pars = ['resp_hat[%d]' % x for x in range(1,self.n+1)]
    ll_pars =   ['log_lik[%d]' % x for x in range(1,self.n+1)]
    llh_pars =  ['log_lik_hat[%d]' % x for x in range(1,self.n+1)]
    pp = (self.stanfit.to_dataframe(pars=rhat_pars,
                                   inc_warmup=False,
                                   diagnostics=False)
          >> select(rhat_pars)
         ).T.values
    ll = (self.stanfit.to_dataframe(pars=ll_pars, 
                                  inc_warmup=False,
                                  diagnostics=False)
          >> select(ll_pars)
         ).T.values
    ll_hat = (self.stanfit.to_dataframe(pars=llh_pars, 
                                  inc_warmup=False,
                                  diagnostics=False)
          >> select(llh_pars)
         ).T.values

    resp = self.data['Response']

    respm = np.tile(resp[:, np.newaxis], (1,pp.shape[1]))

    emp_resid = (2.0*respm-1.0) * np.sqrt(-2.0 * ll)
    sim_resid = (2.0*pp-1.0) * np.sqrt(-2.0 * ll_hat)

    emp_resid_sorted = np.sort(emp_resid, axis=0)
    sim_resid_sorted = np.sort(sim_resid, axis=0)

    cols = ['smp[%d]' % x for x in range(0, pp.shape[1])]
    emp_df = pd.DataFrame(data = emp_resid_sorted,
                          columns = cols)

    emp_df['origin'] = 'empirical'
    sim_df = pd.DataFrame(data = sim_resid_sorted,
                          columns = cols)
    sim_df['origin'] = 'simulated'

    resids_df = emp_df >> bind_rows(sim_df)
    resids_df['sortid'] = resids_df.index/resids_df.shape[0]*2.0
    resids_df >>= gather('smp', 'residuals', cols) 

    def p250(x):
      return np.percentile(x, 2.5)

    def p975(x):
      return np.percentile(x, 97.5)

    ppc_summary = (resids_df >> group_by('sortid', 'origin')
                       >> summarize_each([p250, p975, np.median],
                                         X.residuals)
                  )

    @make_symbolic
    def is_disjoint(emp_low, emp_high, sim_low, sim_high):
      return np.logical_or(np.less(emp_high, sim_low),
                           np.less(sim_high, emp_low))

    pval = (ppc_summary >> mask(X.origin == 'empirical')
            >> inner_join(ppc_summary >> mask(X.origin == 'simulated')
                          >> rename(sim_p250 = X.residuals_p250,
                                    sim_p975 = X.residuals_p975)
                          >> select(X.sim_p250, X.sim_p975, X.sortid),
                         by='sortid')
            >> transmute(disjoint = is_disjoint(X.residuals_p250, X.residuals_p975,
                                                  X.sim_p250, X.sim_p975))
            >> summarize(disjoint_sum = X.disjoint.sum())
           )['disjoint_sum']/ppc_summary.shape[0]

    print('Proportion of non-overlapping intervals:', pval[0])

    return (pval[0], ppc_summary)

  def ppc_residuals_run(self):
    rhat_pars = ['resp_hat[%d]' % x for x in range(1,self.n+1)]
    dd_pars =   ['decision[%d]' % x for x in range(1,self.n+1)]
    pp = (self.stanfit.to_dataframe(pars=rhat_pars,
                                   inc_warmup=False,
                                   diagnostics=False)
          >> select(rhat_pars)
         ).T.values
    dd = (self.stanfit.to_dataframe(pars=dd_pars, 
                                  inc_warmup=False,
                                  diagnostics=False)
          >> select(dd_pars)
         ).T.values

    resp = self.data['Response']
    respm = np.tile(resp[:, np.newaxis], (1,pp.shape[1]))

    def count_revs(a):
      return np.sum(np.abs(np.diff(a, axis=0)), axis=0)

    for i in range(0, dd.shape[1]):
      sort_index = np.argsort(dd[:,i])

      respm[:,i] = respm[sort_index, i]
      pp[:,i] = pp[sort_index, i]

    emp_runs = count_revs(respm)
    sim_runs = count_revs(pp)

    run_df = (pd.DataFrame({'origin': 'empirical', 'runs': emp_runs})
              >> bind_rows(pd.DataFrame({'origin': 'simulated', 'runs': sim_runs}))
             )
    pval = np.sum(np.greater(emp_runs - sim_runs, 0.0))/sim_runs.shape[0]

    print('Proportion of runs larger in empirical data:', pval)

    return (pval, run_df)

  def ppc_flip_count(self):
    rhat_pars = ['resp_hat[%d]' % x for x in range(1,self.n+1)]
    pp = (self.stanfit.to_dataframe(pars=rhat_pars,
                                   inc_warmup=False,
                                   diagnostics=False)
          >> select(rhat_pars)
         ).T.values
    cols = ['smp[%d]' % x for x in range(0, pp.shape[1])]

    def count_revs(a):
      return np.sum(np.abs(np.diff(a)))

    rhat_df = pd.DataFrame(data = pp,
                           columns = cols)

    reversals = (self.data
                 >> bind_cols(rhat_df)
                 >> gather('smp', 'rhat', cols)
                 >> arrange(X.S1, X.S2, X.S3)
                 >> group_by(X.smp)
                 >> summarize_each([count_revs],X.rhat)
                 >> ungroup()
                 >> arrange(X.rhat_count_revs)
                )

    sorted_resp = self.data >> arrange(X.S1, X.S2, X.S3)

    emp_rev = count_revs(sorted_resp['Response'])

    pval = np.searchsorted(reversals['rhat_count_revs'], emp_rev)/reversals.shape[0]

    print('Number of reversals:', emp_rev, '-- p-value:', pval[0])

    return (pval[0], reversals, emp_rev)

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
    (_, ppc_df) = self.ppc_ordered_residuals()
    (_, rev_df, emp_rev) = self.ppc_flip_count()
    (_, run_df) = self.ppc_residuals_run()

    return (gg.ggplot(ppc_df, gg.aes(x='residuals_median', y='sortid', color='origin')) +
       gg.geom_line() +
       gg.geom_line(gg.aes(x='residuals_p250')) + 
       gg.geom_line(gg.aes(x='residuals_p975')),
       gg.ggplot(rev_df, gg.aes(x='rhat_count_revs')) +
       gg.geom_histogram() + gg.geom_vline(gg.aes(xintercept=emp_rev)),
       gg.ggplot(run_df, gg.aes(x='runs', fill='origin')) +
       gg.geom_histogram()
    )

class LpsDifferenceScale(DifferenceScale):
  def __init__(self, stimulus, stanfit, trials):
    super().__init__(stimulus, stanfit, trials)

    self.lapserate = stanfit.summary(pars=['lapses'], probs=[0.025, 0.25, 0.5, 0.75, 0.975])['summary']

  def __str__(self):
    return super().__str__() + "\nlapse rate:\t" + str(self.lapserate[:,0])
