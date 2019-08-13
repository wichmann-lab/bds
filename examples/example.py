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

data = df[['Response', 'i1', 'i2', 'i3']].values

s =df[['s1', 's2', 's3']].values
stimulus = np.unique(s.flatten())


# how to call the functions?

fit = bds.bds(data, lapses=True)

# gives warning
# WARNING:pystan:Maximum (flat) parameter count (1000) exceeded: skipping diagnostic tests for n_eff and Rhat.
# To run all diagnostics call pystan.check_hmc_diagnostics(fit)

pystan.check_hmc_diagnostics(fit.stanfit)

# how to read out results? 
scale = fit.get_scale_values()

CIl, CIh = fit.get_scale_credible_interval()

import matplotlib.pyplot as plt

# posterior predictive checks

for p in fit.diagnostic_plots():
  p.draw();
  plt.show()

# how to plot results?

plt.plot(stimulus,scale,'o')
yerr = [CIh-scale, scale-CIl]
plt.errorbar(stimulus, scale, yerr=yerr, fmt='none', 
                        ecolor='b', capsize=0)

plt.show()
# others:

# list of requirements: pystan, numpy,..
# be python2 compatible too? Now it only runs in python3

# getting python 2 incompatibility error.
# TypeError: super() takes at least 1 argument (0 given)
