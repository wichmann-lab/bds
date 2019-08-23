from .LikelihoodModel import LikelihoodModel
from dfply import *

import plotnine as gg

class BinomialMixture(LikelihoodModel):
  def __init__(self, lapses_prior):

    self.lps_prior = lapses_prior
    self.model_code = """
  for (n in 1:N) {
    target += log_mix(""" + self.lps_prior.name + """,
                      bernoulli_lpmf(Response[n] | 0.5),
                      bernoulli_lpmf(Response[n] | decision[n]));
  }"""

    self.generator_def = (
"""  vector[N] log_lik;
  vector[N] lapse_log_lik;
  vector[N] model_log_lik;

  vector[N] log_lik_hat;
  vector[N] lapse_log_lik_hat;
  vector[N] model_log_lik_hat;
  
  vector[N] resid;
  vector[N] resid_hat;

  real deviance;
  real deviance_hat;

  int resp_hat[N];
  int lapse_hat[N];""")

    self.generator = (
"""  for (n in 1:N) {
    lapse_log_lik[n] = bernoulli_lpmf(Response[n] | 0.5);
    model_log_lik[n] = bernoulli_lpmf(Response[n] | decision[n]);
    
    log_lik[n] = log_mix(""" + self.lps_prior.name + """,
                         lapse_log_lik[n],
                         model_log_lik[n]);

    lapse_hat[n] = bernoulli_rng(""" + self.lps_prior.name + """);
    if (lapse_hat[n] == 0) {
      resp_hat[n] = bernoulli_rng(decision[n]);
    } else {
      resp_hat[n] = bernoulli_rng(0.5);
    }
    lapse_log_lik_hat[n] = bernoulli_lpmf(resp_hat[n] | 0.5);
    model_log_lik_hat[n] = bernoulli_lpmf(resp_hat[n] | decision[n]);
    log_lik_hat[n] = log_mix(""" + self.lps_prior.name + """,
                             lapse_log_lik[n],
                             model_log_lik[n]);

    resid[n] = (2 * Response[n] - 1) * sqrt( -2 * log_lik[n]);
    resid_hat[n] = (2 * resp_hat[n] - 1) * sqrt( -2 * log_lik_hat[n]);
  }

  deviance = sum(square(resid));
  deviance_hat = sum(square(resid_hat));""")

  def register(self, bdsmodel):
    self.lps_prior.register(bdsmodel)
    super().register(bdsmodel)

    bdsmodel._code_blocks['generated quantities:definitions'].append(self.generator_def)
    bdsmodel._code_blocks['generated quantities'].append(self.generator)

  def results(self, result_obj):
    super().results(result_obj)

    result_obj.lapserate = self.lapses(result_obj.stan_fit)

    result_obj.decision_probabilities = self.decision_probabilities(result_obj)
    result_obj.ppc_responses = self.ppc_responses(result_obj)

    pval1, ordered_residuals, plt1 = self.ppc_ordered_residuals(result_obj)
    result_obj.pvalues['ordered residuals'] = pval1
    result_obj.ppc_ordered_residuals = ordered_residuals
    result_obj.diagnostic_plots['ordered residuals'] = plt1

    pval2, runs, plt2 = self.ppc_residual_reversals(result_obj)
    result_obj.pvalues['residual reversals'] = pval2
    result_obj.ppc_residual_reversals = runs
    result_obj.diagnostic_plots['residual reversals'] = plt2

#    pval3, flip_count, emp_rev, plt3 = self.ppc_flip_count(result_obj)
#    result_obj.pvalues['flip count'] = pval3
#    result_obj.summary_statistics['flip count'] = emp_rev
#    result_obj.ppc_flip_count = flip_count
#    result_obj.diagnostic_plots['flip count'] = plt3

  def decision_probabilities(self, result_obj):
    dec_pars = ['decision[%d]' % x for x in range(1,result_obj.n+1)]
    probs = (result_obj.stan_fit.to_dataframe(pars=dec_pars, 
                                  inc_warmup=False,
                                  diagnostics=False)
                 >> select(dec_pars)
                ).values

    return probs

  def lapses(self, stan_fit):
    return stan_fit.extract(pars=self.lps_prior.name)

  def classify_lapses(self, stan_fit):
    pass

  def ppc_responses(self, result_obj):
    resp_pars = ['resp_hat[%d]' % x for x in range(1,result_obj.n+1)]
    responses = (result_obj.stan_fit.to_dataframe(pars=resp_pars, 
                                  inc_warmup=False,
                                  diagnostics=False)
                 >> select(resp_pars)
                ).values

    return responses

  def ppc_ordered_residuals(self, result_obj):
    emp_resid_sorted = np.sort(result_obj.residuals.T, axis=0)
    sim_resid_sorted = np.sort(result_obj.ppc_residuals.T, axis=0)

    cols = ['smp[%d]' % x for x in range(0, emp_resid_sorted.shape[1])]
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

    ppc_plot = (gg.ggplot(ppc_summary, gg.aes(x='residuals_median', y='sortid', color='origin')) +
       gg.geom_line() +
       gg.geom_line(gg.aes(x='residuals_p250')) + 
       gg.geom_line(gg.aes(x='residuals_p975')))

    return (pval[0], ppc_summary, ppc_plot)

  def ppc_residual_reversals(self, result_obj):
    pp = result_obj.ppc_responses.T
    dd = result_obj.decision_probabilities.T

    resp = result_obj.data['Response']
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

    ppc_plot = (gg.ggplot(run_df, gg.aes(x='runs', fill='origin')) +
       gg.geom_histogram(binwidth=2.0))

    return (pval, run_df, ppc_plot)

  def ppc_flip_count(self, result_obj):
    pp = result_obj.ppc_responses.T
    cols = ['smp[%d]' % x for x in range(0, pp.shape[1])]

    def count_revs(a):
      return np.sum(np.abs(np.diff(a)))

    rhat_df = pd.DataFrame(data = pp,
                           columns = cols)

    reversals = (result_obj.data
                 >> bind_cols(rhat_df)
                 >> gather('smp', 'rhat', cols)
                 >> arrange(X.S1, X.S2, X.S3, X.S4)
                 >> group_by(X.smp)
                 >> summarize_each([count_revs],X.rhat)
                 >> ungroup()
                 >> arrange(X.rhat_count_revs)
                )


    sorted_resp = result_obj.data >> arrange(X.S1, X.S2, X.S3, X.S4)

    emp_rev = count_revs(sorted_resp['Response'])

    pval = np.searchsorted(reversals['rhat_count_revs'], emp_rev)/reversals.shape[0]

    print('Number of reversals:', emp_rev, '-- p-value:', pval)

    ppc_plot = (gg.ggplot(reversals, gg.aes(x='rhat_count_revs')) +
       gg.geom_histogram(binwidth=2.0) + gg.geom_vline(gg.aes(xintercept=emp_rev)))
    return (pval, reversals, emp_rev, ppc_plot)

