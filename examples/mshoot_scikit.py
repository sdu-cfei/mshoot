"""
This script runs a multiple-shooting optimization using a scikit model.

Objective: minimize energy consumption.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logging
import mshoot
from mshoot.interfaces.scikit import SimScikit
from examples.models import statespace
from resources.ou44 import ou44_data

# Set up logging
logging.basicConfig(filename='test.log', filemode='w', level='DEBUG')

# Training and validation periods
# ===============================
train_end = 86400 * 9
valid_len = 86400 * 1  # Optimization on this period

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

    def simulate(self, udf, x0):
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
data = ou44_data.get_inputs(period=len_all / 3600., dt=60)
data = data.reset_index()
data.index = (data['datetime'] - data['datetime'][0]).dt.total_seconds()
data['qhvac'] = np.sin(data.index.values * 2 * np.pi / 86400.) * 3000.

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
ml = SimScikit(model="Linear")
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
x0 = [16., 16., 16.]
steps = uvalid.index.size

uguess = np.full((steps - 1, 1), 9999.)

# Define free input and state bounds
ubounds = [
    (0., 4000.)  # qhvac
]

# Temperature bounds
Tr_lo = np.where((data_valid['datetime'].dt.hour >= 5)
    & (data_valid['datetime'].dt.hour <= 7), 21., 16.)

xbounds = [
    (Tr_lo, 25.),  # Room temperature
    (0., 50.),  # Internal thermal mass
    (0., 50.)   # External thermal mass
]

# Optimize
# ========
udf, xdf = mpc.optimize(
    model=ml,
    inp=uvalid,
    free=['qhvac'],
    ubounds=ubounds,
    xbounds=xbounds,
    x0=x0,
    uguess=uguess,
    ynominal=[20., 20., 20., 4000.],
    join=1
)

# Plot results
# ============
xdf.index /= 3600.
udf.index /= 3600.

fig, ax = plt.subplots(2, 1, sharex=True)
ax[0].plot(xdf['x0'], label=r'$T_i$', marker='o')
ax[0].plot(udf.index, Tr_lo, 'k--', label=r'$T_{sp}$')
ax[0].legend()
ax[0].set_ylabel(r'$T$ [$^\circ$C]')
ax[1].plot(udf['qhvac'], label=r'$q_{hvac}$', marker='o')
ax[1].set_xlabel('Time [h]')
ax[1].legend()
ax[1].set_ylabel(r'$q$ [W]')

plt.show()
