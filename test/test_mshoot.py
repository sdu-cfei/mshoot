import unittest

import numpy as np
import pandas as pd
from scipy.signal import StateSpace
import matplotlib.pyplot as plt

import mshoot


# Model
class SSModel(mshoot.SimModel):
    def __init__(self):
        self.ss = self.get_model()

    def get_model(self):
        # State-space R3C3
        # =================
        # Constants
        Vi = 400.  # Indoor volume [m3]

        # Parameters
        Ri = 0.25 / (3. * Vi ** (2./3.))  # Int. wall resistance [K/W]
        Re1 = 1.00 / (3. * Vi ** (2./3.))  # Ext. wall resistance [K/W]
        Re2 = 0.25 / (3. * Vi ** (2./3.))  # Ext. wall resistance [K/W]
        Ci = 30.0 * Vi * 1.2 * 1005.  # Int. wall capacitance [J/K]
        Ce = 20.0 * Vi * 1.2 * 1005.  # Ext. wall capacitance [J/K]
        Cr = 10.0 * Vi * 1.2 * 1005.  # Zone capacitance [J/K]
        fsol = 1.8   # Solar coeff.
        qmet = 100.  # Metabolic heat gain [W/person]

        # States, inputs and outputs:
        # n (states)  : [Tr, Ti, Te]
        # p (inputs)  : [Tout, Hglo, qhvac, nocc]
        # q (outputs) : [Tr, Ti, Te, qhvac]

        # State matrix (n x n)
        A = np.array([
            [-1/Ri/Cr - 1/Re1/Cr, 1/Ri/Cr,             1/Re1/Cr],
            [1/Ri/Ci,            -1/Ri/Ci,                    0],
            [1/Re1/Ce,                  0, -1/Re1/Ce - 1/Re2/Ce]
        ])  # [Tr, Ti, Te]^T

        # Input matrix (n x p)
        B = np.array([
            [0.,       fsol/Cr, 1/Cr, qmet/Cr],
            [0.,            0.,   0.,      0.],
            [1/Re2/Ce,      0.,   0.,      0.]
        ])

        # Output matrix (q x n)
        C = np.array([
            [1., 0., 0.],  # Tr
            [0., 1., 0.],  # Ti
            [0., 0., 1.],  # Te
            [0., 0., 0.]   # qhvac
        ])

        # Feedthrough matrix (q x p)
        D = np.array([
            [0., 0., 0., 0.],  # Tr
            [0., 0., 0., 0.],  # Ti
            [0., 0., 0., 0.],  # Te
            [0., 0., 1., 0.]   # qhvac
        ])

        return StateSpace(A, B, C, D)

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

# Cost function
def cfun(xdf, ydf):
    """
    :param ydf: DataFrame, model states
    :param ydf: DataFrame, model outputs
    :return: float
    """
    qhvac = ydf['qhvac']
    n = qhvac.size
    c = (np.square(qhvac).sum() / n)
    return c

# Unit tests
class TestMShoot(unittest.TestCase):

    def setUp(self):
        self.model = SSModel()

        # Prepare inputs
        index = np.arange(0, 3600*6 + 1, 3600)
        data = np.zeros((index.size, 4))
        data[:, 0] = 0.  # Tout
        data[:, 1] = 0.  # Hglo
        data[:, 2] = 0.  # qhvac
        data[:, 3] = 0.  # nocc
        self.inp = pd.DataFrame(
            data=data,
            columns=['Tout', 'Hglo', 'qhvac', 'nocc'],
            index=index
        )

    def tearDown(self):
        pass

    def test_optimize_const_bounds(self):
        """Test with constant u and x bounds"""
        # Bounds
        ubounds = [(0., 4000.)]
        xbounds = [(20., 23.), (0., 50.), (0., 50.)]
        # Initial state
        x0 = [21., 21., 21.]
        # Instantiate multiple shooting
        optc = mshoot.MShoot(cfun=cfun)
        # Optimize
        udf, xdf = optc.optimize(
            model=self.model,
            inp=self.inp,
            free=['qhvac'],
            ubounds=ubounds,
            xbounds=xbounds,
            ynominal=[20., 20., 20., 4000.],
            x0=x0,
            uguess=None,
            join=1
        )

        # Final room temperature should be equal to lower bound
        self.assertLess(np.abs(xdf.iloc[-1][0] - xbounds[0][0]), 1e-4,
                         "Final state different than lower bound")
        # Heating should be provided
        self.assertGreater(udf['qhvac'].sum(), 0, "qhvac should be positive")
        # No heating at the beginning is needed

    def test_optimize_const_bounds_joined(self):
        """Test if good solution is achieved also with joined intervals"""
        # Bounds
        ubounds = [(0., 4000.)]
        xbounds = [(21., 23.), (0., 50.), (0., 50.)]
        # Initial state
        x0 = [21., 21., 21.]
        # Instantiate multiple shooting
        optc = mshoot.MShoot(cfun=cfun)
        # Optimize
        udf, xdf = optc.optimize(
            model=self.model,
            inp=self.inp,
            free=['qhvac'],
            ubounds=ubounds,
            xbounds=xbounds,
            ynominal=[20., 20., 20., 4000.],
            x0=x0,
            uguess=None,
            join=3
        )
        # Final room temperature should be equal to lower bound
        self.assertLess(np.abs(xdf.iloc[-1][0] - xbounds[0][0]), 1e-4,
                         "Final state different than lower bound")
        # Heating should be provided
        self.assertGreater(udf['qhvac'].sum(), 0, "qhvac should be positive")

    def test_optimize_step_bounds(self):
        """Test preheating (step increase in state lower bound)"""
        # Bounds
        Tlo = np.full(self.inp.index.size, 16.)
        Tlo[-2:] = 19.
        plt.plot(Tlo)
        ubounds = [(0., 5000.)]
        xbounds = [(Tlo, 30.), (0., 50.), (0., 50.)]
        # Initial state
        x0 = [16., 16., 16.]
        # Instantiate multiple shooting
        optc = mshoot.MShoot(cfun=cfun)
        # Optimize
        udf, xdf = optc.optimize(
            model=self.model,
            inp=self.inp,
            free=['qhvac'],
            ubounds=ubounds,
            xbounds=xbounds,
            ynominal=[20., 20., 20., 4000.],
            x0=x0,
            uguess=None,
            join=1
        )
        # Final room temperature should be equal to lower bound
        self.assertLess(np.abs(xdf.iloc[-1][0] - xbounds[0][0][-1]), 1e-4,
                         "Final state different than lower bound")
        # Room temperature should be always higher or equal to the lower bound
        self.assertTrue((xdf['x0'].values >= Tlo).all(), "Temperature bound violated")
        # Heating should be provided
        self.assertGreater(udf['qhvac'].sum(), 0, "qhvac should be positive")

    def test_infeasible_x0(self):
        """Test if good solution is found even with infeasible x0"""
        # Bounds
        ubounds = [(0., 4000.)]
        xbounds = [(21., 23.), (0., 50.), (0., 50.)]
        # Initial state
        x0 = [15., 21., 21.]  # <- infeasible x0
        # Instantiate multiple shooting
        optc = mshoot.MShoot(cfun=cfun)
        # Optimize
        udf, xdf = optc.optimize(
            model=self.model,
            inp=self.inp,
            free=['qhvac'],
            ubounds=ubounds,
            xbounds=xbounds,
            ynominal=[20., 20., 20., 4000.],
            x0=x0,
            uguess=None,
            join=3
        )

        # Final room temperature should be equal to lower bound
        self.assertLess(np.abs(xdf.iloc[-1][0] - xbounds[0][0]), 1e-4,
                         "Final state different than lower bound")
        # Heating should be provided
        self.assertGreater(udf['qhvac'].sum(), 0, "qhvac should be positive")

    def test_normalized_uy(self):
        """Test the right solution is achieved with nominal u and y"""
        # Bounds
        ubounds = [(0., 4000.)]
        xbounds = [(21., 23.), (0., 50.), (0., 50.)]
        # Initial state
        x0 = [21., 21., 21.]
        # Instantiate multiple shooting
        optc = mshoot.MShoot(cfun=cfun)
        # Optimize
        udf, xdf = optc.optimize(
            model=self.model,
            inp=self.inp,
            free=['qhvac'],
            ubounds=ubounds,
            xbounds=xbounds,
            x0=x0,
            uguess=None,
            unominal=[4000.],
            ynominal=[20., 20., 20., 4000.],
            join=3
        )
        # Final room temperature should be equal to lower bound
        self.assertLess(np.abs(xdf.iloc[-1][0] - xbounds[0][0]), 1e-4,
                         "Final state different than lower bound")
        # Heating should be provided
        self.assertGreater(udf['qhvac'].sum(), 0, "qhvac should be positive")

    def test_uguess_outside_bounds(self):
        """Test initial guess for u outsidethe feasible region
        (should be adjusted)"""
        # Bounds
        ubounds = [(0., 4000.)]
        xbounds = [(21., 23.), (0., 50.), (0., 50.)]
        # Initial u
        uguess = np.full((self.inp.index.size - 1, 1), 9999.)
        # Initial state
        x0 = [21., 21., 21.]
        # Instantiate multiple shooting
        optc = mshoot.MShoot(cfun=cfun)
        # Optimize
        udf, xdf = optc.optimize(
            model=self.model,
            inp=self.inp,
            free=['qhvac'],
            ubounds=ubounds,
            xbounds=xbounds,
            x0=x0,
            uguess=uguess,
            ynominal=[20., 20., 20., 4000.],
            join=1,
            maxiter=2
        )
        # Heating power should not exceed 4000 W
        self.assertTrue((udf['qhvac'] < 4000.).all())

if __name__ == "__main__":
    unittest.main()