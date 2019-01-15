from scipy.signal import StateSpace
import numpy as np


def get_model():

    # Build model: R3C3 with a feedthrough price signal
    # =================================================
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
    # p (inputs)  : [Tout, Hglo, qhvac, nocc, price]
    # q (outputs) : [Tr, Ti, Te, qhvac, price]

    # State matrix (n x n)
    A = np.array([
        [-1/Ri/Cr - 1/Re1/Cr, 1/Ri/Cr,             1/Re1/Cr],
        [1/Ri/Ci,            -1/Ri/Ci,                    0],
        [1/Re1/Ce,                  0, -1/Re1/Ce - 1/Re2/Ce]
    ])  # [Tr, Ti, Te]^T

    # Input matrix (n x p)
    B = np.array([
        [0.,       fsol/Cr, 1/Cr, qmet/Cr, 0.],
        [0.,            0.,   0.,      0., 0.],
        [1/Re2/Ce,      0.,   0.,      0., 0.]
    ])

    # Output matrix (q x n)
    C = np.array([
        [1., 0., 0.],  # Tr
        [0., 1., 0.],  # Ti
        [0., 0., 1.],  # Te
        [0., 0., 0.],  # qhvac
        [0., 0., 0.]   # price
    ])

    # Feedthrough matrix (q x p)
    D = np.array([
        [0., 0., 0., 0., 0.],  # Tr
        [0., 0., 0., 0., 0.],  # Ti
        [0., 0., 0., 0., 0.],  # Te
        [0., 0., 1., 0., 0.],  # qhvac
        [0., 0., 0., 0., 1.]   # price
    ])

    return StateSpace(A, B, C, D)
