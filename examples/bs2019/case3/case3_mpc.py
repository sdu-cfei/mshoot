import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logging
import os
import time

import mshoot

# Set up logging
logging.basicConfig(filename='mpc_case3.log', filemode='w', level='DEBUG')

# Random seed
np.random.seed(12345)

# Paths
ms_file = os.path.join('examples', 'bs2019', 'measurements.csv')
fmu = os.path.join('examples', 'bs2019', 'case3', 'models', 'r1c1co2_dymola_1e-11.fmu')

# Simulation period
t0 = '2018-04-05 00:00:00'
t1 = '2018-04-08 00:00:00'

# Read measurements
ms = pd.read_csv(ms_file)
ms['datetime'] = pd.to_datetime(ms['datetime'])
ms = ms.set_index('datetime')
ms = ms.loc[t0:t1]

# Resample
ms = ms.resample('20min').mean().ffill().bfill()

# Assign model inputs
inp = ms[['solrad', 'Tout', 'occ', 'dpos', 'vpos']]
inp['time'] = (inp.index - inp.index[0]).total_seconds()
inp = inp.set_index('time')

# Initial state
airT0 = 20. + 273.15
CO20 = 500.
x0 = [airT0, CO20]

# Cost function
def cfun(xdf, ydf):
    """
    :param ydf: DataFrame, model states
    :param ydf: DataFrame, model outputs
    :return: float
    """
    Qr = ydf['Qr'].iloc[-1]
    Ve = ydf['vetot'].iloc[-1]
    return Qr ** 2 + Ve ** 2

# Iterate over horizons
horizons = [3, 6, 9]
for hrz in horizons:

    par_file = os.path.join('examples', 'bs2019', 'case3', 'results', 'est',
                            'parameters.csv')
    wdir = os.path.join('examples', 'bs2019', 'case3', 'results', 'mpc',
                        'h{}'.format(hrz))

    # Skip the case, if results already there
    if os.path.exists(wdir):
        pass

    else:
        os.makedirs(wdir)

        # Read parameters and modify Tve to allow cooling
        pm = pd.read_csv(par_file)
        parameters = dict()
        for p in pm:
            parameters[p] = pm.iloc[0][p]

        # parameters['Tve'] = 18.  # In real building this temperature is higher.
        #                          # Here, we assume that the air passes through the
        #                          # rotary wheel HX (air preheated to around 18 degC),
        #                          # but not through the heating coil.

        # Instantiate emulation and control models
        model_emu = mshoot.SimFMU(
            fmu,
            outputs=['T', 'Qr', 'vetot'],
            states=['cair.T', 'co2.balance.CO2ppmv_i'],
            parameters=parameters,
            verbose=False)

        model_ctr = mshoot.SimFMU(
            fmu,
            outputs=['T', 'Qr', 'vetot'],
            states=['cair.T', 'co2.balance.CO2ppmv_i'],  # States should be initialized with fixed=False
            parameters=parameters,
            verbose=False)

        # Instantiate MPCEmulation
        mpc = mshoot.MPCEmulation(model_emu, cfun)
        step = 1
        horizon = hrz

        # Contraints
        Tmin = np.where((ms.index.hour >= 8) & (ms.index.hour < 17), 21 + 273.15, 19 + 273.15)
        Tmax = np.where((ms.index.hour >= 8) & (ms.index.hour < 17), 22 + 273.15, 24 + 273.15)
        CO2min = 400.
        CO2max = 800.

        constr = pd.DataFrame(data=np.column_stack((Tmin, Tmax)),
                            columns=['Tmin', 'Tmax'], index=inp.index)
        constr['CO2min'] = CO2min
        constr['CO2max'] = CO2max
        constr.to_csv(os.path.join(wdir, 'constr.csv'))

        # Run
        t0 = time.time()
        u, xctr, xemu, yemu, u_hist = mpc.optimize(
            model=model_ctr,
            inp_ctr=inp.copy(),
            inp_emu=inp.copy(),
            free=['vpos', 'dpos'],
            ubounds=[(0., 100.), (0., 100.)],
            xbounds=[(Tmin, Tmax), (CO2min, CO2max)],
            x0=x0,
            maxiter=150,
            ynominal=[300., 1e7, 2e6],
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
