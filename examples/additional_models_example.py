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

fit_default = bds.bds(df, stimulus)

from bds.plots import *

import matplotlib.pyplot as plt

plot_scale(fit_default).draw();
plt.show()

# While the binomial mixture model is the default, there are other variants available
from bds.additional_models import *

fit_no_lapses = bds_constant_lapserate(df, stimulus, lapserate=0.0)

fit_sens_noise = bds_sensory_noise(df, stimulus)

fit_abs = bds_abs(df, stimulus)

fit_gp = bds_gp(df, stimulus, np.linspace(stimulus[0], stimulus[-1], num=100))

plot_scale(fit_no_lapses).draw();
plot_scale(fit_sens_noise).draw();
plot_scale(fit_abs).draw();
plot_predictive_scale(fit_gp).draw();
plt.show()

scale_fits = {'default': fit_default,
                'no lapses': fit_no_lapses,
                'sensory noise': fit_sens_noise,
                'absolute values': fit_abs,
                'Gaussian process': fit_gp
               }

compare_models_deviance(scale_fits).draw();

# further analysis with arviz

# Compares models with LOO-PSIS --- smaller values are better
compare_models_loo(scale_fits).draw();
plt.show()
