"""
Type: MPC example
Control model: R2C2 FMU
Emulation model: R2C2 FMU
Inputs: same for control and emulation
Objective: minimize energy
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logging
import os
import time

import mshoot

# Set up logging
logging.basicConfig(filename='test.log', filemode='w', level='DEBUG')

# Random seed
np.random.seed(1)

# Setup
# =================

# MPC model definition
fmupath = os.path.join('resources', 'fmus', 'R2C2', 'R2C2.fmu')
parameters = {'C': 1e6, 'R': 0.01}
model_ctr = mshoot.SimFMU(
    fmupath,
    outputs=['qout', 'Tr'],
    states=['heatCapacitor1.T', 'heatCapacitor2.T'],
    parameters=parameters,
    verbose=False)

model_emu = mshoot.SimFMU(
    fmupath,
    outputs=['qout', 'Tr'],
    states=['heatCapacitor1.T', 'heatCapacitor2.T'],
    parameters=parameters,
    verbose=False)

# Output directory
outdir = './examples/out/mpc_fmi_2/'
if not os.path.exists(outdir):
    print('Creating directory: {}'.format(outdir))
    os.makedirs(outdir)

# Inputs
tstep = 3600.  # s
tend = 3600. * 48
t = np.arange(0., tend + 1, tstep)
h = t / 3600.
q = np.full(t.size, 0.)
Tout  = np.sin(t / 86400. * 2. * np.pi) + 273.15

data = pd.DataFrame(
    index=pd.Index(t, name='time'),
    columns=['q', 'Tout'],
    data=np.vstack((q, Tout)).T
)

# Initial state
x0 = [21.1 + 273.15, 21.1 + 273.15]

# Bounds
Tlo = np.where((h >= 8) & (h <= 17), 21. + 273.15, 21. + 273.15)
Thi = np.where((h >= 8) & (h <= 17), 25. + 273.15, 25. + 273.15)

data['Tlo'] = Tlo
data['Thi'] = Thi

# Cost function
def cfun(xdf, ydf):
    """
    :param ydf: DataFrame, model states
    :param ydf: DataFrame, model outputs
    :return: float
    """
    qout = ydf['qout']
    c = np.sum(qout ** 2) / qout.size
    return c


# Instantiate MPCEmulation
mpc = mshoot.MPCEmulation(model_emu, cfun)
step = 1
horizon = 6

t0 = time.time()
u, xctr, xemu, yemu, u_hist = mpc.optimize(
    model=model_ctr,
    inp_ctr=data[['q', 'Tout']],
    inp_emu=data[['q', 'Tout']],
    free=['q'],
    ubounds=[(0., 20000.)],
    xbounds=[(data['Tlo'].values, data['Thi'].values),
             (0. + 273.15, 50. + 273.15)],
    x0=x0,
    maxiter=30,
    ynominal=[20000., 20.],
    step=step,
    horizon=horizon
)
cputime = int(time.time() - t0)

fig, ax = plt.subplots(3, 2, figsize=(14, 8), dpi=100, sharex=True)
# Plot 1
u.plot(grid=True, title='Final optimized trajectory', ax=ax[0][0])
# Plot 2
xctr[['x0', 'x1']].plot(grid=True, title='Control states', ax=ax[0][1])
# Plot 3
xemu[['heatCapacitor1.T', 'heatCapacitor2.T']].plot(grid=True,
                                                    title='Emulation states',
                                                    ax=ax[1][1])
# Plot 4
for ui in u_hist:
    ax[1][0].plot(ui)
ax[1][0].set_title('Optimal solutions from all intervals')
ax[1][0].grid()
# Plot 5
ax[2][0].set_title('Penalty')
ax[2][0].grid()
# Plot 6
ax[2][1].plot(u.index, Tout)
ax[2][1].set_title('Outdoor temperature')
ax[2][1].grid()

fig.suptitle("{}h step, {}h horizon, CPU time {}s"
             .format(step, horizon, cputime))
fig.savefig(os.path.join(outdir, 'mpc_result.png'))

plt.show()
