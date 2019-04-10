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
    MPC emulation is used to test an MPC strategy in a simulated
    environment. Two models are required: (1) control model,
    (2) emulation model.

    The control model is used in the optimization.
    The emulation model is the virtual mock-up of the real system.

    The cost function should return a scalar and have
    the following signature:

    .. code::

        def cost(xdf, ydf):
            # user code here
            return cost


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

    def optimize(self, model, inp_ctr, inp_emu, free, ubounds, xbounds,
                 x0, ynominal=None, maxiter=50,
                 step=1, horizon=10, inp_clb=None):
        """
        Optimize problem using progressive multiple shooting optimizations.

        `inp_ctr` and `inp_emu` have to have the same index (control
        and emulation inputs need to be aligned).

        If `inp_clb` is provided with a callback function, then it overrides
        `inp_ctr`. The function can be used to provide new inputs for each
        receding horizon, emulating in example new weather/occupancy forecasts.
        The callback function takes only one argument `index`, which is a numpy
        1D array with time steps for which the inputs have to be provided,
        e.g. if `index = [300, 600, 900, 1200]` is passed to the function,
        it should return a data frame with index `[300, 600, 900, 1200]`
        and columns for all input variables. Since `inp_ctr` is a required argument,
        for clarity it is advised to set it to `None` in cases when `inp_clb` is provided.
        See the example `examples/mpc_fmi_6.py` for how to use `inp_clb`.

        Return:

        - u - DataFrame, optimal free inputs
        - xctr - DataFrame, control states
        - xemu - DataFrame, emulation states
        - yemu - DataFrame, emulation outputs
        - uhist - list of DataFrames, optimal free inputs per MPC interval

        :param model: SimModel, control model
        :param inp_ctr: DataFrame, fixed inputs for control model, index with time
        :param inp_emu: DataFrame, fixed inputs for emulation model, index with time
        :param free: list, names of free inputs
        :param ubounds: list of tuples of floats, free input bounds
        :param xbounds: list of tuples of vectors, state bounds
        :param x0: list, initial state
        :param ynominal: list, nominal values of outputs (for regularization)
        :param maxiter: int, maximum number of iterations (default 50)
        :param step: int, MPC re-run step - number of `inp` rows (default 1)
        :param horizon: int, opt. horizon - number of `inp` rows (default 10)
        :param inp_clb: callback function, see description
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

        # Initialize inp_ctr from callback, if provided.
        # Initialization based on inp_emu. In each MPC loop
        # pass, a portion of the inp_ctr will be
        # overwritten by the callback function.
        if inp_clb is not None:
            inp_ctr = inp_emu.copy()

        # Assert index type is int
        inp_ctr.index = inp_ctr.index.astype(int)
        inp_emu.index = inp_emu.index.astype(int)

        # Assert control and emulation input indexes as the same
        assert np.isclose(inp_ctr.index.values, inp_emu.index.values).all(), \
            "Control and emulation input indexes are not equal!"

        # Initialize optimal free inputs u (to be returned)
        t = inp_emu.index[0]  # Current time
        u = pd.DataFrame(
            index=inp_emu.index,
            columns=free,
            data=np.zeros((inp_emu.index.size, len(free))))  # Optimal free inputs
        xemu = pd.DataFrame(index=inp_emu.index)  # Emulation states
        yemu = pd.DataFrame(index=inp_emu.index)  # Emulation outputs
        xctr = pd.DataFrame(index=inp_emu.index)  # Control states
        uhist = list()  # List of optimal solutions from all intervals

        # Extend bounds if given as floats
        xbounds = list(xbounds)
        for i in range(len(xbounds)):
            xb = xbounds[i]
            newxb = list()
            for b in xb:
                if isinstance(b, float) or isinstance(b, int):
                    newxb.append(np.full(inp_emu.index.size, b))
                else:
                    newxb.append(b)
            xbounds[i] = newxb

        # Instantiate optimizer
        ms = MShoot(self.cost_user)

        # Start MPC loops
        i = 0
        dt = 0
        while (i + horizon <= inp_emu.index.size):
            self.log.debug("Current MPC time: {} s".format(t))
            print("Progress: {:.1f}%".format(float(i) / inp_emu.index.size * 100.))  # TODO: To log
            # Calculate MPC re-run time step
            dt = inp_emu.index[i+step] - inp_emu.index[i]

            # Define inputs for next period
            nxt_inp_ctr = inp_ctr.iloc[i:i+horizon+1].copy()
            nxt_inp_emu = inp_emu.iloc[i:i+horizon+1].copy()

            # Overwrite control inputs using the callback function (if provided)
            if inp_clb is not None:
                index = inp_ctr.index[i:i+horizon+1]

                self.log.debug(
                    'Using inp_clb for getting new inputs. index = {}'.format(index)
                )

                nxt_inp_ctr = inp_clb(index=index)

                assert np.isclose(nxt_inp_ctr.index.values, index.values).all(), \
                    'Index of the dataframe returned by' + \
                    ' inp_clb not consistent with emulation input dataframe'

            nxt_xbounds = [
                (b[0][i:i+horizon+1], b[1][i:i+horizon+1]) for b in xbounds
            ]

            if i == 0:
                nxt_x0 = x0
            else:
                nxt_x0 = xemu.iloc[i].values

            uguess = u.loc[nxt_inp_ctr.index].iloc[:-1].values  # One element shorter than inp

            self.log.debug("nxt_inp_ctr:\n{}".format(nxt_inp_ctr))
            self.log.debug("nxt_xbounds:\n{}".format(nxt_xbounds))
            self.log.debug("next_x0:\n{}".format(nxt_x0))
            self.log.debug("uguess:\n{}".format(uguess))

            # Optimize
            udf, xdf = ms.optimize(
                model=model,
                inp=nxt_inp_ctr,
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
            nxt_inp_emu.loc[udf.index, free] = udf[free].copy()
            nxt_inp_emu = nxt_inp_emu.dropna()  # Last row of free inp is NaN

            self.log.debug(
                "Progress emulation with...\n"
                "nxt_inp_emu=\n{}\n"
                "nxt_x0=\n{}\n".format(nxt_inp_emu, nxt_x0)
            )
            ey, ex = self.emumod.simulate(nxt_inp_emu, nxt_x0, save_state=True)

            # Assert index type is int
            # (mshoot doesn't have control over the model implementation)
            ey.index = ey.index.astype(int)  # Assert index type is int
            ex.index = ex.index.astype(int)  # Assert index type is int

            # TODO: Figure out why the following is needed
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
