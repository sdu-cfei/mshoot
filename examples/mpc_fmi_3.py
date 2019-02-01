"""
Type: MPC example
Control model: R1C1 FMU
Emulation model: R2C2 FMU
Inputs: same for control and emulation
Objective: minimize price
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

# MPC model definition
fmu_R1C1_path = os.path.join('resources', 'fmus', 'RCModels_R1C1.fmu')
fmu_R2C2_path = os.path.join('resources', 'fmus', 'RCModels_R2C2.fmu')

model_emu = mshoot.SimFMU(
    fmu_R2C2_path,
    outputs=['q', 'T', 'p'],
    states=['heatCapacitorC.T'],
    parameters={'C': 1e6, 'R': 0.01, 'Ci': 1e3},
    verbose=False)

model_ctr = mshoot.SimFMU(
    fmu_R1C1_path,
    outputs=['q', 'T', 'p'],
    states=['heatCapacitorC.T'],
    parameters={'C': 1e6, 'R': 0.01},
    verbose=False)

# Output directory
script_name = os.path.splitext(__file__)[0].split(os.path.sep)[-1]
outdir = os.path.join('.', 'examples', 'out', script_name)
if not os.path.exists(outdir):
    print('Creating directory: {}'.format(outdir))
    os.makedirs(outdir)

# Inputs
tstep = 3600.  # s
tend = 3600. * 24
t = np.arange(0., tend + 1, tstep)
h = t / 3600.
qin = np.full(t.size, 0.)
price = np.where((h > 10) & (h < 16), 10., 1.0)
Tamb  = np.sin(t / 86400. * 2. * np.pi) + 273.15

data = pd.DataFrame(
    index=pd.Index(t, name='time'),
    columns=['qin', 'Tamb', 'price'],
    data=np.vstack((qin, Tamb, price)).T
)

# Initial state
x0 = [21.1 + 273.15]

# Bounds
Tlo = np.where((h >= 8) & (h <= 17), 21. + 273.15, 21. + 273.15)
Thi = np.where((h >= 8) & (h <= 17), 24. + 273.15, 24. + 273.15)

data['Tlo'] = Tlo
data['Thi'] = Thi

# Cost function
def cfun(xdf, ydf):
    """
    :param ydf: DataFrame, model states
    :param ydf: DataFrame, model outputs
    :return: float
    """
    pfinal = ydf['p'].iloc[-1]
    return pfinal

# Instantiate MPCEmulation
mpc = mshoot.MPCEmulation(model_emu, cfun)
step = 1
horizon = 6

t0 = time.time()
u, xctr, xemu, yemu, u_hist = mpc.optimize(
    model=model_ctr,
    inp_ctr=data[['qin', 'Tamb', 'price']],
    inp_emu=data[['qin', 'Tamb', 'price']],
    free=['qin'],
    ubounds=[(0., 5000.)],
    xbounds=[(data['Tlo'].values, data['Thi'].values)],
    x0=x0,
    maxiter=30,
    ynominal=[5000., 20., 5.],
    step=step,
    horizon=horizon
)
cputime = int(time.time() - t0)

# Save data frames
u.to_csv(os.path.join(outdir, 'u.csv'))
xctr.to_csv(os.path.join(outdir, 'xctr.csv'))
xemu.to_csv(os.path.join(outdir, 'xemu.csv'))
yemu.to_csv(os.path.join(outdir, 'yemu.csv'))
for i in range(len(u_hist)):
    u_hist[i].to_csv(os.path.join(outdir, 'u{}.csv'.format(i)))

# Plot results
fig, ax = plt.subplots(3, 2, figsize=(14, 8), dpi=100, sharex=True)
# Plot 1
u.plot(grid=True, title='Final optimized trajectory', ax=ax[0][0])
# Plot 2
xctr[['x0']].plot(grid=True, title='Control states', ax=ax[0][1])
# Plot 3
xemu[['heatCapacitorC.T']].plot(grid=True, title='Emulation states', ax=ax[1][1])
# Plot 4
for ui in u_hist:
    ax[1][0].plot(ui)
ax[1][0].set_title('Optimal solutions from all intervals')
ax[1][0].grid()
# Plot 5
ax[2][0].plot(t, price)
ax[2][0].set_title('Price')
ax[2][0].grid()
# Plot 6
ax[2][1].plot(u.index, Tamb)
ax[2][1].set_title('Outdoor temperature')
ax[2][1].grid()

fig.suptitle("{}h step, {}h horizon, CPU time {}s"
             .format(step, horizon, cputime))
fig.savefig(os.path.join(outdir, 'mpc_result.png'))

plt.show()

# Validate u
inp = data[['qin', 'Tamb']]
inp['qin'] = u['qin']
yvld, xvld = model_emu.simulate(inp, x0)

# Should have only zeroes if the control model matches the emulator model
(yvld - yemu).to_csv(os.path.join(outdir, 'yvld-yemu.csv'))
(xvld - xemu).to_csv(os.path.join(outdir, 'xvld-xemu.csv'))