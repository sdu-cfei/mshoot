import os

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import mshoot

fmupath = os.path.join('examples', 'tutorial', 'modelica', 'R3C3.fmu')

# 1) Emulation model
model_emu = mshoot.SimFMU(
    fmupath,
    outputs=['y', 'u1_y'],
    states=['heatCapacitor1.T', 'heatCapacitor2.T', 'heatCapacitor3.T'],
    parameters={'C1': 75000, 'C2': 100000, 'C3': 50000, 'R1': 0.01, 'R2': 0.01, 'R3': 0.01})

# 2) Control model
model_ctr = mshoot.SimFMU(
    fmupath,
    outputs=['y', 'u1_y'],
    states=['heatCapacitor1.T', 'heatCapacitor2.T', 'heatCapacitor3.T'],
    parameters={'C1': 75000, 'C2': 100000, 'C3': 50000, 'R1': 0.01, 'R2': 0.01, 'R3': 0.01})

# 3) Cost function
def cfun(xdf, ydf):
    cost = (ydf['u1_y'] ** 2).sum()
    return cost

# 4) Define inputs (48 hours emulation, 1h step)
t = np.arange(0, 48. * 3600., 3600.)
u1 = np.zeros(48)
u2 = np.sin(t / 86400. * 2. * np.pi) * 1000. + np.random.rand(48) * 1000. - 500.  # Noisy sinusoid

inp = pd.DataFrame(index = pd.Index(t, name = 'time'))
inp['u1'] = u1
inp['u2'] = u2

inp.plot()
plt.show()

# 5) Define bounds
Tlo = np.where((t > 86400 / 2) & (t < 86400 * 1.5), 273.15 + 23, 273.15 + 17)
Thi = 273.15 + 25

# 6) Instantiate MPCEmulation
mpc = mshoot.MPCEmulation(model_emu, cfun)

# 7) Optimize
u, xctr, xemu, yemu, u_hist = mpc.optimize(
    model = model_ctr,
    inp_ctr = inp,
    inp_emu = inp,
    free = ['u1'],
    ubounds = [(-1000, 1000)],
    xbounds = [(273.15, 333.15), (273.15, 333.15), (Tlo, Thi)],
    x0 = [293.15, 293.15, 293.15],
    maxiter = 20,
    ynominal = [293.15, 1000.],
    step = 1,
    horizon = 3
)

# 8) Plot some results
ax1 = xemu.plot()
ax1.plot(xemu.index, Tlo, color = 'black')
ax1.plot(xemu.index, np.full(48, Thi), color = 'black')

ax2 = u.plot()
ax2.plot(u.index, u2, color = 'red')

plt.show()
