"""
This script runs an MPC example using:
- an state-space model for emulation,
- the same model for control.

Author: K. Arendt
Date: Sep 2018
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logging
import os
import time
import mshoot
from mshoot.mpc import MPCEmulation
from examples.models import statespace_price

# Set up logging
logging.basicConfig(filename='test.log', filemode='w', level='DEBUG')

# Random seed
np.random.seed(1)

# Setup
# =================

# MPC model definition
class SSModel(mshoot.mshoot.SimModel):
    def __init__(self, ss):
        """
        Constructor is optional.
        """
        self.ss = ss

    def simulate(self, udf, x0, **kwargs):
        """
        Method required by SimModel.

        :param udf: DataFrame, inputs
        :param x0: vector, initial state
        :return: ydf (outputs, DataFrame), xdf (states, DataFrame)
        """
        t = udf.index.values
        t_shift = t - t[0]
        u = udf.values
        tout, yout, xout = self.ss.output(u, t_shift, X0=x0)
        ydf = pd.DataFrame(index=t, data=yout,
            columns=['Tr', 'Ti', 'Te', 'qhvac', 'price'])
        xdf = pd.DataFrame(index=t, data=xout,
            columns=['Tr', 'Ti', 'Te'])
        return ydf, xdf

# Output directory
script_name = os.path.splitext(__file__)[0].split(os.path.sep)[-1]
outdir = os.path.join('.', 'examples', 'out', script_name)
if not os.path.exists(outdir):
    print('Creating directory: {}'.format(outdir))
    os.makedirs(outdir)

# Inputs
tstep = 3600.  # s
tend = 3600. * 48
t     = np.arange(0., tend + 1, tstep)
h     = t / tstep
Tout  = np.sin(t / 86400. * 2. * np.pi)
Hglo  = np.full(t.size, 0.)
qhvac = np.full(t.size, 0.)
nocc  = np.full(t.size, 0.)
# price = np.where((h >= 24) & (h <= 28), 2, 1)
price = np.random.rand(t.size)

data = pd.DataFrame(
    index=pd.Index(t, name='time'),
    columns=['Tout', 'Hglo', 'qhvac', 'nocc', 'price'],
    data=np.vstack((Tout, Hglo, qhvac, nocc, price)).T
)

# Initial state
x0 = [21.1, 21.1, 4.2]

# Bounds
Tlo = np.where((h >= 8) & (h <= 17), 21., 21.)
Thi = np.where((h >= 8) & (h <= 17), 24., 24.)

data['Tlo'] = Tlo
data['Thi'] = Thi

# Cost function
def cfun(xdf, ydf):
    """
    :param ydf: DataFrame, model states
    :param ydf: DataFrame, model outputs
    :return: float
    """
    qhvac = ydf['qhvac']
    p = ydf['price']
    n = qhvac.size
    c = np.square(qhvac * p).sum() / n
    return c


# Instantiate model
ssm = statespace_price.get_model()
model = SSModel(ssm)


# Instantiate MPCEmulation
mpc = MPCEmulation(model, cfun)
step = 1
horizon = 6

t0 = time.time()
u, xctr, xemu, yemu, u_hist = mpc.optimize(
    model=model,
    inp=data[['Tout', 'Hglo', 'qhvac', 'nocc', 'price']],
    free=['qhvac'],
    ubounds=[(0., 5000.)],
    xbounds=[(data['Tlo'].values, data['Thi'].values), (0, 50.), (0, 50.)],
    x0=x0,
    maxiter=30,
    ynominal=[20., 20., 20., 5000., 1.],
    step=step,
    horizon=horizon
)
cputime = int(time.time() - t0)

fig, ax = plt.subplots(3, 2, figsize=(14, 8), dpi=100, sharex=True)
# Plot 1
u.plot(grid=True, title='Final optimized trajectory', ax=ax[0][0])
# Plot 2
xctr[['x0']].plot(grid=True, title='Control states', ax=ax[0][1])
# Plot 3
xemu[['Tr']].plot(grid=True, title='Emulation states', ax=ax[1][1])
# Plot 4
for ui in u_hist:
    ax[1][0].plot(ui)
ax[1][0].set_title('Optimal solutions from all intervals')
ax[1][0].grid()
# Plot 5
ax[2][0].plot(u.index, price)
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

print("Emulation outputs:")
print(yemu)

# Validate
data['qhvac'] = u['qhvac']
ydf, xdf = model.simulate(data[['Tout', 'Hglo', 'qhvac', 'nocc', 'price']], x0)

ydf.to_csv(os.path.join(outdir, 'ydf_vld.csv'))
yemu.to_csv(os.path.join(outdir, 'yemu.csv'))
print(ydf)
print(yemu)
