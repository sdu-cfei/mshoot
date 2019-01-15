import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logging
import os
import time

import mshoot

# Set up logging
logging.basicConfig(filename='mpc_case2.log', filemode='w', level='DEBUG')

# Random seed
np.random.seed(12345)

# Paths
ms_file = os.path.join('examples', 'bs2019', 'measurements.csv')
par_file = os.path.join('examples', 'bs2019', 'case2', 'results', 'est',
                        'r1c1_dymola_1e-11', 'parameters.csv')
fmu_file = os.path.join('examples', 'bs2019', 'case2', 'models', 'r1c1_dymola_1e-11_dpos0.fmu')
out_dir = os.path.join('examples', 'bs2019', 'case2', 'results', 'mpc-lin')

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

# SVM training and validation periods
trn_t0 = 0
trn_t1 = trn_t0 + 5 * 86400
vld_t0 = trn_t1
vld_t1 = trn_t0 + 9 * 86400

# MPC simulation period
mpc_t0 = '2018-04-05 00:00:00'
mpc_t1 = '2018-04-08 00:00:00'

# Read measurements
ms = pd.read_csv(ms_file)
ms['datetime'] = pd.to_datetime(ms['datetime'])
ms = ms.set_index('datetime')

# Resample
ms = ms.resample('1h').mean().ffill().bfill()

# Add time [s] column
ms['time'] = (ms.index - ms.index[0]).total_seconds()

# Modify valve position to facilitate training
ms['vpos'] = np.sin(ms['time'].values / 86400. * 2 * np.pi) * 100.

# Assign model inputs
inp = ms[['time', 'solrad', 'Tout', 'occ', 'vpos']]

# Select training data
dat_trn = ms.loc[(ms['time'] >= trn_t0) & (ms['time'] <= trn_t1)]
inp_trn = dat_trn[['time', 'solrad', 'Tout', 'occ', 'vpos']].set_index('time')
inp_trn.index = inp_trn.index - inp_trn.index[0]

# Select validation data
dat_vld = ms.loc[(ms['time'] >= vld_t0) & (ms['time'] <= vld_t1)]
inp_vld = dat_vld[['time', 'solrad', 'Tout', 'occ', 'vpos']].set_index('time')
inp_vld.index = inp_vld.index - inp_vld.index[0]

# Select MPC simulation data
dat_mpc = ms.loc[mpc_t0:mpc_t1]
inp_mpc = dat_mpc[['time', 'solrad', 'Tout', 'occ', 'vpos']].set_index('time')
inp_mpc.index = inp_mpc.index - inp_mpc.index[0]

# Initial state
airT0 = 20. + 273.15
x0 = [airT0]

# Instantiate emulation model for training
pm = pd.read_csv(par_file)
parameters = dict()
for p in pm:
    parameters[p] = pm.iloc[0][p]
parameters['maxHeat'] = 2000.  # [W]

model_emu = mshoot.SimFMU(
    fmu_file,
    outputs=['qr'],
    states=['cair.T'],
    parameters=parameters,
    verbose=False)

# Simulate to get model outputs for training
y_trn, x_trn = model_emu.simulate(inp_trn, x0)
y_vld, x_vld = model_emu.simulate(inp_vld, x0)

# Find optimal SVM model
cgrid = np.arange(0.01, 4., 0.01)
err = list()

for c in cgrid:
    print("Testing C={}".format(c))

    # Instantiate control model to be tested
    model_ctr = mshoot.SimScikit(model="SVM",
        use_state=False,
        kernel='linear',
        C=c,
        epsilon=0.25)

    # SVM training
    model_ctr.train(inp_trn, x_trn[['cair.T']], y_trn[['qr']])

    # SVM validation
    y_emu, x_emu = model_emu.simulate(inp_vld, x0)
    y_ctr, x_ctr = model_ctr.simulate(inp_vld, x0)

    # Validation error
    err_x = (x_emu - x_ctr).abs().sum() / 300.
    err_y = (y_emu - y_ctr).abs().sum() / 2000.
    err.append(err_x.loc['cair.T'] + err_y.loc['qr'])

    # Plot comparison of states and outputs
    # fig, axes = plt.subplots(2, 1, figsize=(5,4))
    # ax = axes[0]
    # ax.plot(x_emu['cair.T'], label='FMU (R1C1)')
    # ax.plot(x_ctr['cair.T'], label='SVM')
    # ax.legend()
    # ax = axes[1]
    # ax.plot(y_emu['qr'], label='FMU (R1C1)')
    # ax.plot(y_ctr['qr'], label='SVM')
    # ax.legend()
    # fig.suptitle("C={}".format(c))
    # plt.show()

# Find optimal C
imin = np.argmin(err)
copt = cgrid[imin]

# Plot C vs. error
plt.plot(cgrid, err)
plt.title("C_opt = {}".format(copt))
plt.show()

# Instantiate final control model
model_ctr = mshoot.SimScikit(model="SVM",
    use_state=False,
    kernel='linear',
    C=copt,
    epsilon=0.25)

# Plot final validation result
model_ctr.train(inp_vld, x_vld[['cair.T']], y_vld[['qr']])
y_emu, x_emu = model_emu.simulate(inp_vld, x0)
y_ctr, x_ctr = model_ctr.simulate(inp_vld, x0)

x_emu.to_csv(os.path.join(out_dir, 'vld_xemu.csv'))
x_ctr.to_csv(os.path.join(out_dir, 'vld_xctr.csv'))
y_emu.to_csv(os.path.join(out_dir, 'vld_yemu.csv'))
y_ctr.to_csv(os.path.join(out_dir, 'vld_yctr.csv'))

fig, axes = plt.subplots(2, 1, figsize=(5,4), sharex=True)
ax = axes[0]
ax.plot(x_emu['cair.T'], label='FMU (R1C1)')
ax.plot(x_ctr['cair.T'], label='SVM')
ax.legend()
ax = axes[1]
ax.plot(y_emu['qr'], label='FMU (R1C1)')
ax.plot(y_ctr['qr'], label='SVM')
ax.legend()
ax.set_xlabel('Time [s]')
fig.suptitle("C={}".format(copt))
fig.set_dpi(120)
fig.savefig(os.path.join(out_dir, 'svm_validation.pdf'))
plt.show()

# Cost function
def cfun(xdf, ydf):
    """
    :param ydf: DataFrame, model states
    :param ydf: DataFrame, model outputs
    :return: float
    """
    Qr = (ydf['qr'] ** 2).sum()
    return Qr

# Loop over horizons
horizons = [2, 4, 6, 8, 10]

# Slice measurements (for easy access to index in constraint definition)
ms = ms.loc[mpc_t0:mpc_t1]

for hrz in horizons:

    # Working directory
    wdir = os.path.join(out_dir, 'h{}'.format(hrz))

    if os.path.exists(wdir):
        pass
    else:
        os.makedirs(wdir)

        # Instantiate MPCEmulation
        mpc = mshoot.MPCEmulation(model_emu, cfun)
        step = 1
        horizon = hrz

        # Contraints
        Tmin = np.where((ms.index.hour >= 8) & (ms.index.hour < 17), 21 + 273.15, 19 + 273.15)
        Tmax = np.where((ms.index.hour >= 8) & (ms.index.hour < 17), 22 + 273.15, 24 + 273.15)

        constr = pd.DataFrame(data=np.column_stack((Tmin, Tmax)),
                            columns=['Tmin', 'Tmax'], index=inp_mpc.index)
        constr.to_csv(os.path.join(wdir, 'constr.csv'))

        # Run
        t0 = time.time()
        u, xctr, xemu, yemu, u_hist = mpc.optimize(
            model=model_ctr,
            inp=inp_mpc,
            free=['vpos'],
            ubounds=[(-100., 100.)],
            xbounds=[(Tmin, Tmax)],
            x0=x0,
            maxiter=50,
            ynominal=[2000.],
            step=step,
            horizon=horizon
        )
        cputime = int(time.time() - t0)

        # Save results
        u.to_csv(os.path.join(wdir, 'u.csv'))
        xctr.to_csv(os.path.join(wdir, 'xctr.csv'))
        xemu.to_csv(os.path.join(wdir, 'xemu.csv'))
        yemu.to_csv(os.path.join(wdir, 'yemu.csv'))

        for i in range(len(u_hist)):
            u_hist[i].to_csv(os.path.join(wdir, 'u{}.csv'.format(i)))

        with open(os.path.join(wdir, 'cputime.txt'), 'w') as f:
            f.write("CPU time: {} s".format(cputime))