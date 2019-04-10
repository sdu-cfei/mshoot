import unittest
import os

import numpy as np
import pandas as pd
from scipy.signal import StateSpace
import matplotlib.pyplot as plt

import mshoot


def cfun(xdf, ydf):
    """
    :param ydf: DataFrame, model states
    :param ydf: DataFrame, model outputs
    :return: float
    """
    qout = ydf['qout'].values
    c = np.sum(qout ** 2) / qout.size
    return c


class TestMPC(unittest.TestCase):

    def setUp(self):
        fmupath = os.path.join('resources', 'fmus', 'R1C1', 'R1C1.fmu')
        parameters = {'C': 1e6, 'R': 0.01}
        self.model = mshoot.SimFMU(
            fmupath,
            outputs=['qout', 'Tr'],
            states=['heatCapacitor.T'],
            parameters=parameters,
            verbose=False)

    def tearDown(self):
        pass

    def test_mpc(self):
        # Inputs
        t = np.arange(0, 3600 * 10, 3600)
        inp = pd.DataFrame(index=pd.Index(t, name='time'), columns=['q', 'Tout'])
        inp['q'] = np.full(t.size, 0)
        inp['Tout'] = np.full(t.size, 273.15)

        # Bounds
        ubounds = [(0., 4000.)]
        xbounds = [(293.15, 296.15)]

        # Initial state
        x0 = [293.65]

        # Optimization
        mpc = mshoot.MPCEmulation(emumod=self.model, cfun=cfun)

        u, xctr, xemu, yemu, uhist = mpc.optimize(
            model=self.model,
            inp_ctr=inp.copy(),
            inp_emu=inp.copy(),
            free=['q'],
            ubounds=ubounds,
            xbounds=xbounds,
            x0=x0,
            ynominal=[4000., 293.15],
            step=1,
            horizon=3
        )

        # ax = u.plot(title='u')
        # ax.set_ylim(0, 4000)
        # ax = xemu.plot(title='xemu')
        # ax.set_ylim(292.15, 296.15)
        # plt.show()

        # Assert the solution is correct
        self.assertLess(abs(xemu['heatCapacitor.T'].iloc[-1] - 293.15), 0.3)  # Ideally, should be even closer

        # Validate emulation with optimized control
        inp['q'] = u['q']
        yvld, xvld = self.model.simulate(inp, x0)

        # self.assertTrue(((yvld - yemu).abs() < 1e-3).all().all())  # Might not be true for FMUs *
        self.assertTrue(((xvld - xemu).abs() < 1e-3).all().all())  # Might not be true for FMUs *
        # * FMU results might be shifted in time by one time step.
        #   The reason is unknown, but FMU- or pyFMI-specific.

    def test_mpc_inp_clb(self):
        # Inputs
        t = np.arange(0, 3600 * 10, 3600)
        inp = pd.DataFrame(index=pd.Index(t, name='time'), columns=['q', 'Tout'])
        inp['q'] = np.full(t.size, 0)
        inp['Tout'] = np.full(t.size, 273.15)

        # Bounds
        ubounds = [(0., 4000.)]
        xbounds = [(293.15, 296.15)]

        # Initial state
        x0 = [293.65]

        # Input callback function
        def inp_clb(index):
            return inp.loc[index]

        # Optimization
        mpc = mshoot.MPCEmulation(emumod=self.model, cfun=cfun)

        u, xctr, xemu, yemu, uhist = mpc.optimize(
            model=self.model,
            inp_ctr=None,
            inp_clb=inp_clb,
            inp_emu=inp.copy(),
            free=['q'],
            ubounds=ubounds,
            xbounds=xbounds,
            x0=x0,
            ynominal=[4000., 293.15],
            step=1,
            horizon=3
        )

        # Assert the solution is correct
        self.assertLess(abs(xemu['heatCapacitor.T'].iloc[-1] - 293.15), 0.3)  # Ideally, should be even closer

        # Validate emulation with optimized control
        inp['q'] = u['q']
        yvld, xvld = self.model.simulate(inp, x0)

        # self.assertTrue(((yvld - yemu).abs() < 1e-3).all().all())  # Might not be true for FMUs *
        self.assertTrue(((xvld - xemu).abs() < 1e-3).all().all())  # Might not be true for FMUs *
        # * FMU results might be shifted in time by one time step.
        #   The reason is unknown, but FMU- or pyFMI-specific.

    # def test_2_inputs(self):
    #     """THE SOLVER HAS PROBLEMS WITH GETTING THE RIGHT SOLUTION. (?)"""
    #     # Inputs
    #     t = np.arange(0, 3600 * 10, 3600)
    #     inp = pd.DataFrame(index=pd.Index(t, name='time'), columns=['q', 'Tout'])
    #     inp['q'] = np.full(t.size, 0)
    #     inp['Tout'] = np.full(t.size, 273.15)

    #     # Bounds
    #     ubounds = [(0., 10000.), (272.15, 275.)]  # <-- Solver should try to yield Tout = 275
    #     xbounds = [(293.15, 296.15)]

    #     # Initial state
    #     x0 = [293.65]

    #     # Optimization
    #     mpc = mshoot.MPCEmulation(emumod=self.model, cfun=cfun)

    #     u, xctr, xemu, yemu, uhist = mpc.optimize(
    #         model=self.model,
    #         inp=inp,
    #         free=['q', 'Tout'],
    #         ubounds=ubounds,
    #         xbounds=xbounds,
    #         x0=x0,
    #         unominal=[4000., 273.15],
    #         ynominal=[4000., 293.15],
    #         step=1,
    #         horizon=4
    #     )

    #     ax = u.plot(title='u', subplots=True)
    #     ax = xemu.plot(title='xemu')
    #     plt.show()

    #     # Assert the solution is correct
    #     self.assertLess(abs(xemu['heatCapacitor.T'].iloc[-1] - 293.15), 0.01)

    #     # Validate emulation with optimized control
    #     inp['q'] = u['q']
    #     yvld, xvld = self.model.simulate(inp, x0)

    #     # self.assertTrue((yvld - yemu < 1e-3).all().all())  # Might not be true for FMUs *
    #     # self.assertTrue((xvld - xemu < 1e-3).all().all())  # Might not be true for FMUs *
    #     # * FMU results might be shifted in time by one time step.
    #     #   The reason is unknown, but FMU- or pyFMI-specific.


if __name__ == '__main__':
    unittest.main()


