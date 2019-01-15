import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logging
import os
import time

import mshoot

# Set up logging
logging.basicConfig(filename='mpc_case3.log', filemode='w', level='DEBUG')

# Random seed
np.random.seed(12345)

# Paths
ms_file = os.path.join('examples', 'bs2019', 'measurements.csv')
fmu = os.path.join('examples', 'bs2019', 'case3', 'models', 'r1c1co2pid_dymola_1e-11.fmu')

# Simulation period
t0 = '2018-04-05 00:00:00'
t1 = '2018-04-08 00:00:00'

# Read measurements
ms = pd.read_csv(ms_file)
ms['datetime'] = pd.to_datetime(ms['datetime'])
ms = ms.set_index('datetime')
ms = ms.loc[t0:t1]

# Resample
ms = ms.resample('1h').mean().ffill().bfill()

# Assign model inputs
inp = ms[['solrad', 'Tout', 'occ']]
inp['time'] = (inp.index - inp.index[0]).total_seconds()
inp = inp.set_index('time')
inp.index = inp.index.astype(int)

constr = pd.read_csv('examples/bs2019/case3/results/mpc/h2/constr.csv', index_col=0)
constr.index = constr.index.astype(int)

inp['tstp'] = constr['Tmin']
inp['co2stp'] = constr['CO2max']

print(inp)

# Initial state
airT0 = 20. + 273.15
CO20 = 500.
x0 = [airT0, CO20]

# Parameters and working directory
par_file = os.path.join('examples', 'bs2019', 'case3', 'results', 'est',
                        'parameters.csv')
wdir = os.path.join('examples', 'bs2019', 'case3', 'results', 'pid')

# Skip the case, if results already there
if not os.path.exists(wdir):
    os.makedirs(wdir)

# Read parameters and modify Tve to allow cooling
pm = pd.read_csv(par_file)
parameters = dict()
for p in pm:
    parameters[p] = pm.iloc[0][p]

# parameters['Tve'] = 18.  # In real building this temperature is higher.
#                             # Here, we assume that the air passes through the
#                             # rotary wheel HX (air preheated to around 18 degC),
#                             # but not through the heating coil.

# Instantiate emulation and control models
model_pid = mshoot.SimFMU(
    fmu,
    outputs=['T', 'Qr', 'vetot', 'vpos', 'dpos'],
    states=['cair.T', 'co2.balance.CO2ppmv_i'],
    parameters=parameters,
    verbose=False)

# Run
ydf, xdf = model_pid.simulate(inp, x0)

# Save results
ydf.to_csv(os.path.join(wdir, 'ydf.csv'))
xdf.to_csv(os.path.join(wdir, 'xdf.csv'))
