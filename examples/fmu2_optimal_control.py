"""
This script runs a multiple shooting example using
an FMI model (R2C2).

Objective: minimize energy consumption.

Author: K. Arendt
Date: Sep 2018
"""

import os
import logging

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import mshoot

# Set up logging
logging.basicConfig(filename='test.log', filemode='w', level='DEBUG')

def cfun(xdf, ydf):
    """
    :param ydf: DataFrame, model states
    :param ydf: DataFrame, model outputs
    :return: float
    """
    qout = ydf['qout'].values
    c = np.sum(qout ** 2) / qout.size
    return c

def get_model():
    fmupath = os.path.join('resources', 'fmus', 'R2C2', 'R2C2.fmu')
    parameters = {'C': 1e5, 'R': 0.001}
    model = mshoot.SimFMU(
        fmupath,
        outputs=['qout', 'Tr'],
        states=['heatCapacitor1.T', 'heatCapacitor2.T'],
        parameters=parameters,
        verbose=False)
    return model

# Inputs
t = np.arange(0, 3600 * 10, 3600)
inp = pd.DataFrame(index=pd.Index(t, name='time'), columns=['q', 'Tout'])
inp['q'] = np.full(t.size, 0.)
inp['Tout'] = np.full(t.size, 273.15)

# Bounds
ubounds = [(0., 4000.)]
xbounds = [(293.15, 296.15), (293.15, 296.15)]

# Initial state
x0 = [293.65, 293.35]

# Optimization
mpc = mshoot.MShoot(cfun=cfun)
udf, xdf = mpc.optimize(
    model=get_model(),
    inp=inp,
    free=['q'],
    ubounds=ubounds,
    xbounds=xbounds,
    x0=x0,
    uguess=None,
    ynominal=[4000, 295.],
    join=1,
    maxiter=30
)

# Show results
print(udf)
print(xdf)

ax = udf.plot(title='udf')
ax.set_ylim(0, 4000)
ax = xdf.plot(title='xdf')
ax.set_ylim(293.15, 296.15)
plt.show()
