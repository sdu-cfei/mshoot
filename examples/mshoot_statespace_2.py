"""
This examples presents a multiple shooting optimization
of a state-space model.

Objective: energy cost minimization.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logging
import sys
import mshoot
from examples.models import statespace_price
from resources.ou44 import ou44_data

# Set up logging
logging.basicConfig(filename='test.log', filemode='w', level='DEBUG')

# Setup
# =================

# MPC model definition
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
            columns=['Tr', 'Ti', 'Te', 'qhvac', 'price'])
        xdf = pd.DataFrame(index=t, data=xout,
            columns=['Tr', 'Ti', 'Te'])
        return ydf, xdf

# Read inputs: period [h], dt [min]
data = ou44_data.get_inputs(period=12, dt=90)
data = data.reset_index()
data.index = (data['datetime'] - data['datetime'][0]).dt.total_seconds()

inp = data[['Tout', 'Hglo', 'qhvac', 'nocc']]
inp['price'] = np.linspace(0., 100., inp.index.size)

# Bounds
Tb = data[['datetime', 'nocc']].copy()
Tb['Tr_lo'] = np.where((Tb['datetime'].dt.hour >= 5)
                     & (Tb['datetime'].dt.hour <= 7), 21., 16.)

Tr_lo = Tb['Tr_lo'].values

# Instantiate control and emulation models
ssmodel = statespace_price.get_model()
ctrmod = SSModel(ssmodel)

# Define problem
# =================

# Define cost function
def cfun(xdf, ydf):
    """
    :param ydf: DataFrame, model states
    :param ydf: DataFrame, model outputs
    :return: float
    """
    qhvac = ydf['qhvac']
    price = ydf['price']
    n = qhvac.size
    c = np.square(qhvac*price).sum() / n
    return c

# Instantiate multiple shooting MPC
mpc = mshoot.MShoot(cfun=cfun)
x0 = [16., 16., 16.]
steps = inp.index.size

# Define free input and state bounds
ubounds = [
    (0., 10000.)  # qhvac
]

xbounds = [
    (Tr_lo, 28.),  # Room temperature
    (-100., 100.),  # Internal thermal mass
    (-100., 100.)   # External thermal mass
]

# Optimize
udf, xdf = mpc.optimize(
    model=ctrmod,
    inp=inp,
    free=['qhvac'],
    ubounds=ubounds,
    xbounds=xbounds,
    x0=x0,
    uguess=None,
    ynominal=[20., 20., 20., 4000., 100.],
    join=1
)

# Plot results
# ============
xdf.index /= 3600.
udf.index /= 3600.
Tb.index /= 3600.
inp.index /= 3600.

fig, ax = plt.subplots(3, 1, sharex=True)
ax[0].plot(xdf['x0'], label=r'$T_i$', marker='o')
ax[0].plot(Tb['Tr_lo'], 'k--', label=r'$T_{sp}$')
ax[0].legend()
ax[0].set_ylabel(r'$T$ [$^\circ$C]')
ax[1].plot(udf['qhvac'], label=r'$q_{hvac}$', marker='o')
ax[1].legend()
ax[1].set_ylabel(r'$q$ [W]')
ax[2].plot(inp['price'], label='Price', marker='o')
ax[2].legend()
ax[2].set_ylabel(r'$p$ [-]')
ax[2].set_xlabel('Time [h]')

plt.show()
