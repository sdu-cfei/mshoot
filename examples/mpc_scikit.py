"""
This script runs an MPC example using:
- a scikit model for control,
- a state-space model for emulation.

Author: K. Arendt
Date: Sep 2018
"""

import time
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logging
import mshoot
from mshoot import SimScikit
from mshoot import MPCEmulation
from examples.models import statespace
from resources.ou44 import ou44_data

# Set up logging
logging.basicConfig(filename='test.log', filemode='w', level='DEBUG')

# Output directory
# ================
script_name = os.path.splitext(__file__)[0].split(os.path.sep)[-1]
outdir = os.path.join('.', 'examples', 'out', script_name)
if not os.path.exists(outdir):
    os.makedirs(outdir)

# Training and validation periods
# ===============================
train_end = 86400 * 15
valid_len = 86400 * 2  # Optimization on this period

len_all = train_end + valid_len

# Get "ground truth" data using the state space model
# ===================================================

# State space model
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
            columns=['Tr', 'Ti', 'Te', 'qhvac'])
        xdf = pd.DataFrame(index=t, data=xout,
            columns=['Tr', 'Ti', 'Te'])
        return ydf, xdf

# Read inputs: period [h], dt [min]
data = ou44_data.get_inputs(period=len_all / 3600., dt=30)
data = data.reset_index()
data.index = (data['datetime'] - data['datetime'][0]).dt.total_seconds()

# Non-zero qhvac
data['qhvac'] = np.sin(data.index.values * 2 * np.pi / 86400.) * 5000.

udf_ss = data[['Tout', 'Hglo', 'qhvac', 'nocc']]

data_train = data.loc[:train_end]
data_valid = data.loc[train_end:]

# Instantiate control and emulation models
model_ss = SSModel(statespace.get_model())

ydf_ss, xdf_ss = model_ss.simulate(udf_ss, np.array([20., 20., 20.]))

# Train scikit models
# ===================
utrain = udf_ss.loc[:train_end]
xtrain = xdf_ss.loc[:train_end]
ytrain = ydf_ss.loc[:train_end]

uvalid = udf_ss.loc[train_end:]

# Choose ML model:
# ml = SimScikit(model="Linear")
ml = SimScikit(model="SVM", use_state=True, C=10.)
ml.train(utrain, xtrain, ytrain)

# Simulate scikit models
# ======================
ydf_ml, xdf_ml = ml.simulate(udf_ss, xdf_ss.iloc[0])

# Change index from [s] to [h]
# ============================
def s2h(*args):
    for df in args:
        df.index /= 3600.

s2h(xdf_ss, xdf_ml)

# Compare state space with linear regression
# ==========================================
ax = xdf_ss.plot(color=['red', 'blue', 'green'])
ylim = ax.get_ylim()
ax.set_ylim(ylim)
ax.plot(xdf_ml['Tr'], label='Tr(ML)', ls='-.', c='red')
ax.plot(xdf_ml['Ti'], label='Ti(ML)', ls='-.', c='blue')
ax.plot(xdf_ml['Te'], label='Te(ML)', ls='-.', c='green')
ax.vlines(x=train_end/3600., ymin=ylim[0], ymax=ylim[1], linestyles='dashed', color='k')
ax.set_xticks(np.arange(0, xdf_ss.index.max() + 1, 24))
ax.set_xlabel('Time [h]')
ax.set_ylabel(r'Temperature [$^\circ$C]')
ax.legend()
ax.text(x=train_end/3600.-3, y=ylim[1]-0.5, s='Training', horizontalalignment='right')
ax.text(x=train_end/3600.+3, y=ylim[1]-0.5, s='Validation', horizontalalignment='left')
ax.set_title('LTI State-Space Model vs. Machine Learning Model')

fig = ax.get_figure()
fig.set_size_inches(8, 6)
fig.savefig(os.path.join(outdir, 'model_validation.png'), dpi=200)

# udf_ss.plot(subplots=True)
plt.show()

# Define optimization problem
# ===========================

# Cost function
def cfun(xdf, ydf):
    """
    :param ydf: DataFrame, model states
    :param ydf: DataFrame, model outputs
    :return: float
    """
    qhvac = ydf['qhvac']
    n = qhvac.size
    c = np.square(qhvac).sum() / n
    return c

# Instantiate multiple shooting MPC
mpc = mshoot.MShoot(cfun=cfun)
x0 = [19., 19., 19.]
steps = uvalid.index.size

uguess = np.full((steps - 1, 1), 9999.)

# Temperature bounds
Tr_lo = np.where((data_valid['datetime'].dt.hour >= 8)
    & (data_valid['datetime'].dt.hour <= 16), 21., 17.)

xbounds = [
    (Tr_lo, 25.),  # Room temperature
    (0., 50.),  # Internal thermal mass
    (0., 50.)   # External thermal mass
]

# Run MPC
# =======
mpc = MPCEmulation(model_ss, cfun)
step = 1
horizon = 6

t0 = time.time()
u, xctr, xemu, yemu, u_hist = mpc.optimize(
    model=ml,
    inp=uvalid,
    free=['qhvac'],
    ubounds=[(0., 5000.)],
    xbounds=[(Tr_lo, 30.), (0, 50.), (0, 50.)],
    x0=x0,
    maxiter=30,
    ynominal=[20., 20., 20., 5000.],
    step=step,
    horizon=horizon
)
cputime = int(time.time() - t0)

# Change index from [s] to [h]
# ============================
s2h(u, xctr, xemu, yemu, *u_hist)

# Plots
# =====
fig, ax = plt.subplots(3, 1, figsize=(14, 8), dpi=100, sharex=True)
# Plot 1
u.plot(grid=True, title='Final optimized trajectory', ax=ax[0])
# Plot 2
for ui in u_hist:
    ax[1].plot(ui)
ax[1].set_title('Optimal solutions from all intervals')
ax[1].grid()
# Plot 3
ax[2].set_title("Indoor temperature")
ax[2].plot(xctr.index, xctr['x0'].values, label="Control model")
ax[2].plot(xemu.index, xemu['Tr'].values, label="Emulation model")
ax[2].plot(xemu.index, Tr_lo, label="Setpoint", ls=':', c='k')
ax[2].grid()
ax[2].legend(loc='upper right')
ax[2].set_xlabel('Time [h]')

fig.suptitle("{}h step, {}h horizon, CPU time {}s"
             .format(step, horizon, cputime))
fig.savefig(os.path.join(outdir, 'mpc_result.png'))

plt.show()
