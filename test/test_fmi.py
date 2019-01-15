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


class TestFMI(unittest.TestCase):

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

    def test_simulation(self):
        print('** Test FMU simulation **')
        # Inputs
        t = np.arange(0, 86401, 3600)
        udf = pd.DataFrame(index=pd.Index(t, name='time'), columns=['q', 'Tout'])
        udf['q'] = np.full(t.size, 0)
        udf['Tout'] = np.full(t.size, 273.15)

        # Initial state
        x0 = [273.15 + 20]

        # Simulation
        ydf, xdf = self.model.simulate(udf, x0)
        print('Outputs:')
        print(ydf)
        print('States')
        print(xdf)

        self.assertEqual(xdf.loc[0, 'heatCapacitor.T'], x0[0])
        self.assertLess(xdf.iloc[-1]['heatCapacitor.T'], 273.16)
        self.assertTrue((xdf['heatCapacitor.T'] == ydf['Tr']).all())

    def test_set_parameters(self):
        fmupath = os.path.join('resources', 'fmus', 'R1C1', 'R1C1.fmu')
        parameters = {'C': 2e6, 'R': 0.05}
        model2 = mshoot.SimFMU(
            fmupath,
            outputs=['qout', 'Tr'],
            states=['heatCapacitor.T'],
            parameters=parameters,
            verbose=False)
        
        # Inputs
        t = np.arange(0, 86401, 3600)
        udf = pd.DataFrame(index=pd.Index(t, name='time'), columns=['q', 'Tout'])
        udf['q'] = np.full(t.size, 0)
        udf['Tout'] = np.full(t.size, 273.15)

        # Initial state
        x0 = [273.15 + 20]

        # Simulation
        ydf1, xdf1 = self.model.simulate(udf, x0)
        ydf2, xdf2 = model2.simulate(udf, x0)

        self.assertFalse((ydf1 == ydf2).all().all())
        self.assertFalse((xdf1 == xdf2).all().all())

    def test_optimal_control(self):
        print('** Test FMU optimal control **')
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
        mpc = mshoot.MShoot(cfun=cfun)
        udf, xdf = mpc.optimize(
            model=self.model,
            inp=inp,
            free=['q'],
            ubounds=ubounds,
            xbounds=xbounds,
            x0=x0,
            uguess=None,
            unominal=[4000.],
            ynominal=[4000., 293.15],
        )

        # ax = udf.plot(title='udf')
        # ax.set_ylim(0, 4000)
        # ax = xdf.plot(title='xdf')
        # ax.set_ylim(293.15, 296.15)
        # plt.show()

        self.assertLess(abs(xdf['x0'].iloc[-1] - 293.15), 0.01)
        self.assertLess(abs(udf['q'].iloc[-1] - 2000.), 0.01)


if __name__ == '__main__':
    unittest.main()


