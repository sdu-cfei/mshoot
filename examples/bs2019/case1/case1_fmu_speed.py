import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logging
import os
import time

import mshoot

# Random seed
np.random.seed(12345)

# Paths
ms_file = os.path.join('examples', 'bs2019', 'measurements.csv')
fmu_dir = os.path.join('examples', 'bs2019', 'case1', 'models')

# FMU list
fmus = os.listdir(fmu_dir)

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
inp = ms[['solrad', 'Tout', 'occ', 'dpos', 'vpos']]
inp['time'] = (inp.index - inp.index[0]).total_seconds()
inp = inp.set_index('time')

# Modify inputs
inp['dpos'] = 0

# Initial state
airT0 = 20. + 273.15
x0 = [airT0]

# Computational time array
tols = [int(x.split('_')[-1].split('.')[0].split('-')[-1]) for x in fmus]
tols = np.array(sorted(tols))
ct = pd.Series(index=pd.Index(tols, name='tol'))  # Index: '4' means 1e-4, '11' means 1e-11

# Iterate over FMUs
for fmu in fmus:

    fmu_name = fmu.split('.')[0]
    par_file = os.path.join('examples', 'bs2019', 'case1', 'results', 'est',
                            fmu_name, 'parameters.csv')

    # Skip the case, if results already there
    fmu_file = os.path.join(fmu_dir, fmu)

    # Read parameters and modify heating power (used for heating/cooling in this example)
    pm = pd.read_csv(par_file)
    parameters = dict()
    for p in pm:
        parameters[p] = pm.iloc[0][p]
    parameters['maxHeat'] = 2000.  # [W]

    # Instantiate emulation and control models
    model = mshoot.SimFMU(
        fmu_file,
        outputs=['T', 'Qr', 'vetot'],
        states=['cair.T'],
        parameters=parameters,
        verbose=False)

    # Test speed
    t0 = time.time()
    for i in range(3000):
        ydf, xdf = model.simulate(inp, x0)
    cputime = time.time() - t0

    ct.loc[int(fmu_name.split('-')[-1])] = cputime

print(ct)

