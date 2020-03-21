#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example call of BDS

@author: G. Aguilar, Feb 2019
"""

import bds
import pandas as pd
import numpy as np
import pystan

df = pd.read_csv('test.csv', sep=' ')

# need to document what is the structure of the data to pass as input
# order of columns
# all int? 

df.rename(columns={'i1': 'S1', 'i2': 'S2', 'i3': 'S3'}, inplace=True)

s =df[['s1', 's2', 's3']].values
stimulus = np.unique(s.flatten())


# how to call the functions?

fit = bds.bds(df, stimulus)
#fit = bds.gp_bds(df, stimulus, np.linspace(stimulus[0], stimulus[-1], num=100))

# gives warning
# WARNING:pystan:Maximum (flat) parameter count (1000) exceeded: skipping diagnostic tests for n_eff and Rhat.
# To run all diagnostics call pystan.check_hmc_diagnostics(fit)

#pystan.check_hmc_diagnostics(fit.stanfit)

# how to read out results? 
scale = fit.get_scale_values()

CIl, CIh = fit.get_scale_credible_interval()

from bds.plots import *

import matplotlib.pyplot as plt

# result plots

plot_scale(fit).draw();
plt.show()

plot_posterior_samples(fit).draw();
plt.show()

# posterior predictive checks

plot_deviance(fit).draw();
plt.show()

plot_ordered_residuals(fit).draw();
plt.show()

plot_reversals(fit).draw();
plt.show()

plot_flip_count(fit).draw();
plt.show()

#import bds.additional_models
