import logging
import abc
import time
import copy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mshoot.mshoot import MShoot


class MPCEmulation():
    """
    MPC emulation is used to test the MPC strategy in a simulated
    environment. Two models are required: (1) control model,
    (2) emulation model.

    The control model is used in the optimization.
    The emulation model is the virtual mock-up of the real system.

    The cost function should return a scalar and have
    the following signature:

    ```
        def cost(xdf, ydf):
            # user code here
            return cost
    ```

    where xdf is the state data frame and ydf is the output
    data frame.

    :param emumod: SimModel, emulation model
    :param cfun: user cost function
    """
    def __init__(self, emumod, cfun):
        self.log = logging.getLogger('MPCEmulation')
        self.log.info('Instantiate MPCEmulation')
        # Emulation model
        self.emumod = emumod
        # User cost function
        self.cost_user = cfun

        pd.set_option('display.max_columns', 50)

    def optimize(self, model, inp, free, ubounds, xbounds,
                 x0, ynominal=None, maxiter=50,
                 step=1, horizon=10):
        """
        Optimize problem using progressive multiple shooting optimizations.

        Return:
        - u - DataFrame, optimal free inputs
        - xctr - DataFrame, control states
        - xemu - DataFrame, emulation states
        - yemu - DataFrame, emulation outputs
        - uhist - list of DataFrames, optimal free inputs per MPC interval

        :param model: SimModel, control model
        :param inp: DataFrame, fixed inputs, index with time
        :param free: list, names of free inputs
        :param ubounds: list of tuples of floats, free input bounds
        :param xbounds: list of tuples of vectors, state bounds
        :param x0: list, initial state
        :param ynominal: list, nominal values of outputs (for regularization)
        :param maxiter: int, maximum number of iterations (default 50)
        :param step: int, MPC re-run step - number of `inp` rows (default 1)
        :param horizon: int, opt. horizon - number of `inp` rows (default 10)
        :return: u, xctr, xemu, yemu, uhist
        """
        self.log.info("**Start MPC optimization**")

        # Sanity checks
        assert len(ubounds) == len(free)
        assert len(xbounds) == len(x0)
        if ynominal is not None:
            pass  # No information about model outputs passed to this method

        # Make sure x0 is ndarray
        x0 = np.array(x0).astype(float)

        # Assert index type is int
        inp.index = inp.index.astype(int)

        t = inp.index[0]  # Current time
        u = pd.DataFrame(
            index=inp.index,
            columns=free,
            data=np.zeros((inp.index.size, len(free))))  # Optimal free inputs
        xemu = pd.DataFrame(index=inp.index)  # Emulation states
        yemu = pd.DataFrame(index=inp.index)  # Emulation outputs
        xctr = pd.DataFrame(index=inp.index)  # Control states
        uhist = list()  # List of optimal solutions from all intervals

        # Extend bounds if given as floats
        xbounds = list(xbounds)
        for i in range(len(xbounds)):
            xb = xbounds[i]
            newxb = list()
            for b in xb:
                if isinstance(b, float) or isinstance(b, int):
                    newxb.append(np.full(inp.index.size, b))
                else:
                    newxb.append(b)
            xbounds[i] = newxb

        # Instantiate optimizer
        ms = MShoot(self.cost_user)

        # Start MPC loops
        i = 0
        dt = 0
        while (i + horizon <= inp.index.size):
            self.log.debug("Current MPC time: {} s".format(t))
            print("Progress: {:.1f}%".format(float(i) / inp.index.size * 100.))  # TODO: To log
            # Calculate MPC re-run time step
            dt = inp.index[i+step] - inp.index[i]

            # Define inputs for next period
            nxt_inp = inp.iloc[i:i+horizon+1].copy()
            nxt_xbounds = [
                (b[0][i:i+horizon+1], b[1][i:i+horizon+1]) for b in xbounds
            ]
            if i == 0:
                nxt_x0 = x0
            else:
                nxt_x0 = xemu.iloc[i].values
            uguess = u.loc[nxt_inp.index].iloc[:-1].values  # One element shorter than inp

            self.log.debug("Next inputs:\n{}".format(nxt_inp))
            self.log.debug("Next xbounds:\n{}".format(nxt_xbounds))
            self.log.debug("Next x0:\n{}".format(nxt_x0))
            self.log.debug("Next uguess:\n{}".format(uguess))

            # Optimize
            udf, xdf = ms.optimize(
                model=model,
                inp=nxt_inp,
                free=free,
                ubounds=ubounds,
                xbounds=nxt_xbounds,
                x0=nxt_x0,
                uguess=uguess,
                ynominal=ynominal,
                join=1,  # TODO: add to arguments
                maxiter=maxiter
            )

            # Assert index type is int
            udf.index = udf.index.astype(int)
            xdf.index = xdf.index.astype(int)

            # Add udf to `uhist`
            uhist.append(udf.copy())

            # Update `u`
            u.loc[udf.index] = udf.copy()

            # Save control states to `xctr`
            if len(xctr.columns) == 0:
                # Add columns
                for c in xdf.columns:
                    xctr[c] = np.nan
            xctr.loc[xdf.index] = xdf.copy()

            # Progress emulation by `step`
            nxt_inp.loc[udf.index, free] = udf[free].copy()
            nxt_inp = nxt_inp.dropna()  # Last row of free inp is NaN
            self.log.debug(
                "Progress emulation with...\n"
                "nxt_inp=\n{}\n"
                "nxt_x0=\n{}\n".format(nxt_inp, nxt_x0)
            )
            ey, ex = self.emumod.simulate(nxt_inp, nxt_x0, save_state=True)

            # ...or re-simulate entire period each time
            # inp.loc[udf.index, free] = udf[free].copy()
            # ey, ex = self.emumod.simulate(inp, x0)

            ey.index = ey.index.astype(int)  # Assert index type is int
            ex.index = ex.index.astype(int)  # Assert index type is int
            if len(xemu.columns) == 0:
                # Add columns
                for c in ex.columns:
                    xemu[c] = np.nan
            if len(yemu.columns) == 0:
                # Add columns
                for c in ey.columns:
                    yemu[c] = np.nan

            # Save emulation results to `xemu` and `yemu`
            # Comment:
            #   For unknown reason, ex.index sometimes contains elements
            #   not present in xemu.index... Therefore, select intersection
            intersect = [ti for ti in ex.index if ti in xemu.index]
            xemu.loc[intersect] = ex.loc[intersect].copy()

            intersect = [ti for ti in ey.index if ti in yemu.index]
            yemu.loc[intersect] = ey.loc[intersect].copy()

            self.log.info("Updated xemu:\n{}".format(xemu))
            self.log.info("Updated xctr:\n{}".format(xctr))

            # Increase `t` and `i`
            t += dt
            i += step

        return u, xctr, xemu, yemu, uhist
