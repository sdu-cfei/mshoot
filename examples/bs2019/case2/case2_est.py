"""Parameter estimation in all FMUs used in Case 1"""
import os

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import modestpy

# Paths
ms_file = os.path.join('examples', 'bs2019', 'measurements.csv')
fmu_dir = os.path.join('examples', 'bs2019', 'case2', 'models')
res_dir = os.path.join('examples', 'bs2019', 'case2', 'results', 'est')

if not os.path.exists(res_dir):
    os.makedirs(res_dir)

# FMU list
fmus = os.listdir(fmu_dir)

# Training and validation periods
trn_t0 = 0
trn_t1 = trn_t0 + 5 * 86400  # 4
vld_t0 = trn_t0
vld_t1 = trn_t0 + 9 * 86400

# Read measurements
ms = pd.read_csv(ms_file)
ms['datetime'] = pd.to_datetime(ms['datetime'])
ms = ms.set_index('datetime')

# Resample
ms = ms.resample('1h').mean().ffill().bfill()

# Assign model inputs
inp = ms[['solrad', 'Tout', 'occ', 'dpos', 'vpos']]
inp['time'] = (inp.index - inp.index[0]).total_seconds()  # ModestPy needs index in seconds
inp = inp.set_index('time')                               # ModestPy needs index named 'time'
inp.to_csv(os.path.join(res_dir, 'inp.csv'))

ax = inp.loc[trn_t0:trn_t1].plot(subplots=True)
fig = ax[0].get_figure()
fig.savefig(os.path.join(res_dir, 'inp_training.png'), dpi=200)

ax = inp.loc[vld_t0:vld_t1].plot(subplots=True)
fig = ax[0].get_figure()
fig.savefig(os.path.join(res_dir, 'inp_validation.png'), dpi=200)

# Assign model desired outputs
ideal = ms[['T']]
ideal['time'] = (ideal.index - ideal.index[0]).total_seconds()  # ModestPy needs index in seconds
ideal = ideal.set_index('time')                                 # ModestPy needs index named 'time'
ideal.to_csv(os.path.join(res_dir, 'ideal.csv'))

ax = ideal.loc[trn_t0:trn_t1].plot(subplots=True)
fig = ax[0].get_figure()
fig.savefig(os.path.join(res_dir, 'ideal_training.png'), dpi=200)

ax = ideal.loc[vld_t0:vld_t1].plot(subplots=True)
fig = ax[0].get_figure()
fig.savefig(os.path.join(res_dir, 'ideal_validation.png'), dpi=200)

# Parameters
known = {
    'Vi': 139. * 3.5,
    'maxHeat': 2689.,
    'maxVent': 4800.,
    'Tve': 21.
}

est = dict()
est['shgc'] = (1.0, 0.0, 10.0)
est['tmass'] = (5., 1., 50.)
est['RExt'] = (1., 0.5, 4.)
est['occheff'] = (1., 0.5, 3.0)

# Initial condition parameters:
ic_param = dict()  # Empty, because MShoot needs to manipulate states directly

# Estimation
ga_opts = {'maxiter': 50, 'tol': 1e-7, 'lhs': True, 'pop_size': 40}
scipy_opts = {
    'solver': 'L-BFGS-B',
    'options': {'maxiter': 50, 'tol': 1e-12}
    }

# Iterate over all FMUs
for fmu in fmus:

    wdir = os.path.join(res_dir, fmu.split('.')[0])
    fmu_file = os.path.join(fmu_dir, fmu)

    if not os.path.exists(wdir):
        os.makedirs(wdir)

    if 'dpos0' in fmu:
        inp = inp.drop('dpos', axis=1)

    session = modestpy.Estimation(wdir, fmu_file, inp, known, est, ideal,
        lp_n = 1,
        lp_len = trn_t1 - trn_t0,
        lp_frame = (trn_t0, trn_t1),
        vp = (vld_t0, vld_t1),
        methods = ('GA', 'SCIPY'),
        ga_opts = ga_opts,
        scipy_opts = scipy_opts,
        ic_param=ic_param,
        ftype = 'RMSE',
        seed = 12345)

    estimates = session.estimate()

    # Validation
    vld = session.validate()
    vld_err = vld[0]
    vld_res = vld[1]

    with open(os.path.join(wdir, 'vld_err.txt'), 'w') as f:
        for k in vld_err:
            f.write("{}: {:.5f}\n".format(k, vld_err[k]))

    vld_res.to_csv(os.path.join(wdir, 'vld_res.csv'))

    # Save all parameters (except IC parameters)
    parameters = pd.DataFrame(index=[0])
    for p in estimates:
        parameters[p] = estimates[p]
    for p in known:
        parameters[p] = known[p]

    for p in ic_param:
        parameters = parameters.drop(p, axis=1)

    parameters.to_csv(os.path.join(wdir, 'parameters.csv'), index=False)

    # Check how the estimates are far from the bounds (relative estimates -> esrel)
    esrel = pd.DataFrame(index=[0])
    for p in estimates:
        lb = est[p][1]  # Lower bound
        ub = est[p][2]  # Upper bound
        esrel[p] = (estimates[p] - lb) / (ub - lb)

    esrel.to_csv(os.path.join(wdir, 'parameters_rel.csv'), index=False)
