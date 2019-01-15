import logging
import abc
import time
import copy

import numpy as np
import pandas as pd
from multiprocessing.pool import ThreadPool
from scipy.optimize import fmin_slsqp
import matplotlib.pyplot as plt

from mshoot.optimize import optimize


class SimModel(abc.ABC):
    """MShoot statefull model interface."""

    @abc.abstractmethod
    def simulate(self, udf, x0, **kwargs):
        """
        Simulate the model using the provided inputs `udf`
        and initial state `x0`.

        The DataFrame should have the following content:
        - index - time in seconds and equal steps,
        - columns - input data,
        - column names - input variable names.

        Return two DataFrames, `ydf` and `xdf`, with
        outputs and states, respectively, and with the same
        structure as `udf`.

        :param udf: DataFrame, shape (n_steps, n_variables)
        :param x0: vector, size (n_states, )
        :param **kwargs: Additional arguments required by some interfaces
        :return: ydf, xdf
        """
        pass


class MShoot():

    def __init__(self, cfun):
        self.log = logging.getLogger('MShoot')
        self.log.info('Instantiate MShoot')
        # User cost function
        self.cost_user = cfun

    def test_input_const(self, u, unames, model, inp, x0):
        """
        Test input `u`. `u` is assumed to be constant, so must
        be given as a 1D array (vector) of length ``len(unames)``.

        :param u: vector, constant inputs
        :param unames: list(str)
        :param model: SimModel
        :param inp: DataFrame
        :param x0: vector
        :return: ydf (DataFrame), xdf (DataFrame)
        """
        inp = inp.copy()
        nfree = len(unames)
        for n, i in zip(unames, np.arange(nfree)):
            inp[n] = u[i]
        ydf, xdf = model.simulate(inp, x0)
        return ydf, xdf

    def test_input(self, u, unames, model, inp, x0):
        """
        Test input `u`. `u` is assumed to be time-varying, but is
        given as a flattened array (vector) of length:
        ``len(unames) * inp.index.size``.

        :param u: vector, flattened inputs
        :param unames: list(str)
        :param model: SimModel
        :param inp: DataFrame
        :param x0: vector
        :return: ydf (DataFrame), xdf (DataFrame)
        """
        assert u.shape[0] == inp.shape[0], \
            'Sizes of u ({}) and inp ({}) are incompatible' \
            .format(u.shape[0], inp.shape[0])
        inp = inp.copy()
        # Reshape `u` and add to `inp`
        nfree = len(unames)
        u2D = u.reshape((-1, nfree))
        for n, i in zip(unames, np.arange(nfree)):
            inp[n] = u2D[:, i]
        # Calculate and return ydf and xdf
        ydf, xdf = model.simulate(inp, x0)
        return ydf, xdf

    def extend_bounds(self, bounds, n):
        """
        Unfold scalar bounds to constant vectors.
        All scalar bounds are copied to n-element arrays.

        :param bounds: list(tuple)
        :paran n: int
        :return: array, indexing: [variable][lo/up][interval]
        """
        for i in range(len(bounds)):
            lb = bounds[i][0]  # Lower bound
            ub = bounds[i][1]  # Upper bound

            # If lb or ub is float or int, change to array with n-elements
            if isinstance(lb, (int, float)):
                lb = np.full(n, lb)
            if isinstance(ub, (int, float)):
                ub = np.full(n, ub)
            bounds[i] = (lb, ub)

        return np.array(bounds)

    def slice_inputs(self, inp, nrows):
        """
        Slice `inp` into even adjacent parts with `nrows` each.

        :param inp: DataFrame
        :param nrows: int
        :return: list of DataFrames
        """
        slinp = list()  # Inputs for each subperiod
        i = 0
        rows_left = inp.index.size
        while rows_left - 1 > 0:
            if rows_left >= nrows:
                slinp.append(inp.iloc[i:i+nrows])
                # "-1" to get adjacent subperiods
                i += nrows - 1
                rows_left -= nrows - 1
            else:
                self.log.warning(
                    "{} rows of `inp` not included in optimization"
                    .format(rows_left - 1)
                )
                break
        msg = "Sliced inputs ({}):\n".format(len(slinp))
        for si, i in zip(slinp, np.arange(len(slinp))):
            msg += "[{}]\n{}\n".format(i, si)
        self.log.info(msg)

        return slinp

    def slice_bounds(self, bounds, nrows):
        """
        Slice `bounds` into ``n = steps // (nrows-1)`` parts.

        Each bound must be given as a tuple of vectors. Each vector element
        represents the bound for a specific time.

        Return:

            list(list(tuple(float, float)))
            |    |    |     ^ lo.  ^ up. bound
            |    |    ^ variable
            |    ^ period
            ^ list of periods

        :param bounds: list(tuple(vector, vector))
        :param nrows: int
        :return: list(list(tuple(float, float)))
        """
        slb = list()  # Sliced bounds
        steps = bounds[0][0].size
        i = 0

        while i < steps:
            # "-1" to get adjacent subperiods
            if i % (nrows - 1) == 0:
                row = [[b[0][i], b[1][i]] for b in bounds]
                slb.append(row)
            i += 1

        return slb

    def trim_to_bounds(self, xt, bounds):
        """
        Returns xt if lies within bounds, otherwise returns
        the nearest bound.

        :param xt: vector, state at time t
        :param bounds: list(tuple), bounds for each state
        :return: vector
        """
        rsv = 1e-3

        xt = np.copy(xt)
        for i in range(len(bounds)):
            lb = bounds[i][0]
            ub = bounds[i][1]
            if xt[i] < lb:
                xt[i] = lb + lb * rsv
            elif xt[i] > ub:
                xt[i] = ub - ub * rsv
        return xt

    def x0_within_bounds(self, x0, bounds):
        """
        Return a tuple of two tuples, each containing
        one bool per state. The booleans are True
        is respective bounds are not violated.

        The first tuple refers to the lower bound.
        The second tuple refers to the upper bound.

        :param xt: vector, state at time t
        :param bounds: list(tuple), bounds for each state
        :return: tuple(tuple(bool, ...), tuple(bool, ...))
        """
        lo_ok = [True for x in x0]
        hi_ok = [True for x in x0]

        for i in range(len(x0)):
            lb = bounds[i][0]
            ub = bounds[i][1]
            if x0[i] < lb:
                lo_ok[i] = False
            elif x0[i] > ub:
                hi_ok[i] = False

        lo_ok = tuple(lo_ok)
        hi_ok = tuple(hi_ok)

        return (lo_ok, hi_ok)

    def norm_ubounds(self, ubounds, unominal):
        """
        Normalize bounds formatted as:
        [(lo1, hi1), (lo2, hi2), ..., (lon, hin)]

        Return in the same format.
        """
        nub = list()
        i = 0
        for b in ubounds:
            nub.append((b[0] / unominal[i], b[1] / unominal[i]))
            i += 1
        return nub

    def optimize(self, model, inp, free, ubounds, xbounds,
                 x0, uguess=None, unominal=None, ynominal=None,
                 join=1, maxiter=50):
        """
        Multiple shooting in which the problem is transcribed
        into a single sparse NLP with adjoint variable continuity
        conditions.

        `ubounds` has to be a list of 2-tuples of floats (one tuple
        per free input).

        `xbounds` has to be a list of 2-tuples of floats or vectors
        (one tuple per state). Floats are used for constant bounds.
        Vectors can be used to define time-varying bounds. The size
        of a vector must be N+1, where N is the number of rows in `inp`.

        `join` can be used to solve on a coarser mesh. By default
        `join` equals to 1, meaning each interval is treated separately.
        On larger problems it may be usefull to merge intervals though
        (less adjoint variables).

        Return a tuple with 2 DataFrames: `udf` with optimized free
        inputs and `xdf` with resulting adjoint states. Note that
        `xdf` has one more row than `udf` (for final state).

        If `x0` doesn't lie in the feasible region, it will be automatically
        adjusted.

        The shape of `uguess` should be (n_intervals, n_free_inputs).
        The number of intervals equals to ``inp.index.size - 1``.

        :param model: SimModel, control model
        :param inp: DataFrame, fixed inputs, index with time
        :param free: list, names of free inputs
        :param ubounds: list of tuples, free input bounds
        :param xbounds: list of tuples, state bounds
        :param x0: 1D array, initial state
        :param uguess: 2D array, initial guess for free inputs
        :param ynominal: list, nominal values of outputs (for regularization)
        :param join: int, number of intervals to join (default 1)
        :param maxiter: int, maximum number of iterations (default 50)
        :return: udf, xdf (two DataFrames with optimized inputs and states)
        """
        self.log.info("Start multiple shooting optimization ({}-{} s)"
                      .format(inp.index[0], inp.index[-1]))

        time0 = time.time()

        # Sanity checks
        assert len(ubounds) == len(free)
        assert len(xbounds) == len(x0)
        if unominal is not None:
            assert len(unominal) == len(free)
        if ynominal is not None:
            pass  # No information about model outputs passed to this method

        # Assert index type is int
        inp.index = inp.index.astype(int)

        # Copy lists
        xbounds = xbounds.copy()
        ubounds = ubounds.copy()

        # Make sure all inputs are floats
        inp = inp.astype(float)  # TODO: validate

        # Make sure x0 is ndarray
        x0 = np.array(x0).astype(float)

        # Parameters
        n_states = x0.size
        n_free = len(free)

        # Number of inp rows in each interval:
        # the lower, the easier to solve, but more iterations
        n_rows = join + 1

        # Input and output time
        tin = inp.index.values
        tout = tin[np.arange(0, tin.size, n_rows - 1)]

        self.log.debug('inp ({} rows):\n{}'.format(inp.index.size, inp))
        self.log.debug('ubounds:\n{}'.format(ubounds))
        self.log.debug('xbounds:\n{}'.format(xbounds))

        # If state bounds are given as scalars, unfold to constant vectors
        # indexing: [variable][lo/up][interval]
        xbounds = self.extend_bounds(xbounds, n=inp.index.size)

        # NOT USED ==================================================
        # # State nominal values based on upper bounds
        # xnominal = np.ones(n_states)

        # for i in range(n_states):
        #     xnominal[i] = xbounds[i][1].max()

        # # Normalized state bounds
        # xbnorm = np.ones(xbounds.shape)

        # for i in range(n_states):
        #     xbnorm[i] = xbounds[i] / xnominal[i]
        # ===========================================================

        # Input u nominal values based on max absolute bounds
        unominal = np.ones(n_free)

        for i in range(n_free):
            unominal[i] = max(abs(ubounds[i][0]), abs(ubounds[i][1]))

        # Set y nominal values to 1 if not provided
        if ynominal is None:
            ynominal = 1.
        else:
            ynominal = np.array(ynominal)

        # Normalize u bounds
        self.log.debug("Normalizing ubounds and xbounds")
        ubnorm = self.norm_ubounds(ubounds, unominal)
        self.log.debug("Normalized ubounds:\n{}".format(ubounds))

        # Slice inp, xbounds
        slinp = self.slice_inputs(inp, n_rows)
        slxb = self.slice_bounds(xbounds, n_rows)

        # Number of subperiods
        n_interv = len(slinp)

        # Unfold ubnorm
        ubnorm_unfold = [b for b in ubnorm for i in range(n_interv)]

        # Merge bounds (ux)
        xb = list()  # xbounds reshaped for the needs of uxbounds
        for b in slxb:
            xb.extend(b)
        uxbounds = ubnorm_unfold + xb
        self.log.debug("Merged bounds (u+x): {}".format(uxbounds))

        # Assert x0 is within xbounds and if not, enforce it
        x0_ok = self.x0_within_bounds(x0, slxb[0])

        all_x0_ok = np.array(x0_ok).all()

        if not all_x0_ok:
            msg = "x0 outside bounds ({} not in {}).".format(x0, slxb[0])
            self.log.warning(msg)

            # Relax state bounds
            for i in range(n_states):
                if x0_ok[0][i] is False:
                    slxb[0][i][0] = x0[i] - 1e-3 * x0[i]  # Single-interval relaxation
                    slxb[1][i][0] = (slxb[1][i][0] + slxb[0][i][0]) / 2.
                elif x0_ok[1][i] is False:
                    slxb[0][i][1] = x0[i] + 1e-3 * x0[i]  # Single-interval relaxation
                    slxb[1][i][1] = (slxb[1][i][1] + slxb[0][i][1]) / 2.

            msg = "RELAXING STATE CONSTRAINTS: {}".format(slxb)
            self.log.warning(msg)

        # Initial guess for adjoint states
        shape_xadj = (n_interv + 1, n_states)
        xadj = np.full(shape_xadj, 0.)
        xadj[0] = x0

        # Assert xadj within bounds
        for i in range(n_interv + 1):
            xadj[i] = self.trim_to_bounds(xadj[i], slxb[i])

        self.log.debug("Initial guess for adjoint states:\n{}".format(xadj))

        # Flatten xadj
        xadj = xadj.flatten()

        # Initial guess for u
        shape_u = (n_interv, n_free)
        size_u = shape_u[0] * shape_u[1]

        if uguess is None:
            self.log.debug("uguess is None -> assume zeros")
            uguess = np.zeros(shape_u)
        elif uguess.shape != shape_u:
            msg = "Incorrect shape of uguess: "
            msg += "is {}, should be {}".format(uguess.shape, shape_u)
            if join > 1:
                msg += " (note that join={})".format(join)
            self.log.error(msg)
            raise ValueError(msg)

        self.log.debug("uguess=\n{}".format(uguess))

        # Normalize u
        ugnorm = uguess / unominal
        self.log.debug("Normalized uguess=\n{}".format(ugnorm))

        # Assert normalized uguess within normalized bounds
        for i in range(ugnorm.shape[1]):
            ugnorm[:, i] = np.where(ugnorm[:, i] < ubnorm[i][0],
                                    ubnorm[i][0], ugnorm[:, i])
            ugnorm[:, i] = np.where(ugnorm[:, i] > ubnorm[i][1],
                                    ubnorm[i][1], ugnorm[:, i])

        self.log.debug("Normalied uguess adjusted to bounds:\n{}".format(ugnorm))

        # Flatten u
        u = ugnorm.flatten()

        assert u.size == size_u, "Incorrect size of uguess " \
            "({} instead of {})".format(u.size, size_u)
        self.log.debug("Flattened uguess: {}".format(u))

        # Merge xadj and u into ux
        ux = np.hstack((u, xadj))
        self.log.debug("Initial guess for u+x: {}".format(ux))

        # Cost and constraint functions arguments
        args = (self, n_interv, free, n_free, model, slinp,
                shape_u, shape_xadj, x0, n_states,
                unominal, ynominal)
        MShoot.args = args

        # Optimize using multiple shooting
        # with constraints on u and x and user cost function
        opt, fx, its, imode, smode = optimize(
            func=MShoot.cfun,
            x0=ux,
            jac=MShoot.jac,
            bounds=uxbounds,
            f_eqcons=MShoot.f_eqcons,
            args=args,
            iter=maxiter
        )

        msg = "==================================\n"
        msg = "Solution complete...\n"
        msg += "- iterations:   {}\n".format(its)
        msg += "- final cost:   {}\n".format(fx)
        msg += "- exit mode:    {}\n".format(imode)
        msg += "- exit message: {}\n".format(smode)
        self.log.info(msg)

        # Get u and xadj from the result array
        u = opt[:size_u].reshape(shape_u)
        xadj = opt[size_u:].reshape(shape_xadj)

        # Filling strategy
        # udf has one row less than xdf -> how to align u with the index?
        # 'prepend' - the first row is duplicated
        # 'append'  - the last row is duplicated
        # 'linear'  - the profile is stretched out
        fill_stg = 'linear'

        if fill_stg == 'prepend':
            ures = np.concatenate((u[[0]], u), axis=0)
        elif fill_stg == 'append':
            ures = np.concatenate((u, u[[-1]]), axis=0)
        elif fill_stg == 'linear':
            ures = np.concatenate((u[[0]], u), axis=0)
            for i in range(1, u.shape[0]):
                ures[i] = np.mean((u[[i - 1]], u[[i]]), axis=0)

        # Put outputs into DataFrames
        udf = pd.DataFrame(ures, columns=free,
                           index=pd.Index(tout, name='time'))
        xdf = pd.DataFrame(xadj,
                           columns=['x{}'.format(i) for i in range(n_states)],
                           index=pd.Index(tout, name='time'))

        # Denormalize
        udf = udf * unominal

        self.log.info("u=\n{}".format(udf))
        self.log.info("x=\n{}".format(xdf))
        self.log.info('Computational time: {}s'.format(time.time() - time0))

        return udf, xdf

    @staticmethod
    def cfun(ux, *args):
        """
        Cost function applied to all intervals.

        `ux` contains free inputs and adjoint states.
        """
        # Unpack arguments
        self, n_interv, free, n_free, model, slinp, shape_u, \
            shape_xadj, x0, n_states, unominal, ynominal = args

        # Extract `u` and `xadj` from `ux`
        size_u = shape_u[0] * shape_u[1]
        u = ux[:size_u].reshape(shape_u)
        u = u.reshape((n_interv, n_free))
        xadj = ux[size_u:].reshape(shape_xadj)

        # Denormalize inputs with unominal
        u = u * unominal

        # Initialize vector for interval costs
        costs = np.zeros(n_interv)

        # Sequential simulations
        # ======================
        for i in range(n_interv):
            ydf, xdf = self.test_input_const(
                u=u[i],
                unames=free,
                model=model,
                inp=slinp[i],
                x0=xadj[i])
            # Normalize output with ynominal
            ydf = ydf / ynominal
            # Calculate cost using user-defined function
            # and based on output DataFrame `ydf`
            costs[i] = self.cost_user(xdf, ydf)

        # Return the mean of interval costs
        return costs.mean()

    @staticmethod
    def jac(ux, *args):
        """
        Calculate Jacobian of the cost function.

        `ux` contains free inputs and adjoint states.
        """
        delta = 1e-9   # Step size
        j = np.zeros(ux.size)
        f0 = MShoot.cfun(ux, *args)
        fi = np.zeros(ux.size)

        # Unpack arguments
        self, n_interv, free, n_free, model, slinp, shape_u, \
            shape_xadj, x0, n_states, unominal, ynominal = args

        # Input size
        size_u = shape_u[0] * shape_u[1]

        # Calculate only with respect to inputs (and not adjoint states)
        for i in range(size_u):
            uxi = ux.copy()
            uxi[i] = uxi[i] + delta
            fi[i] = MShoot.cfun(uxi, *args)
            j[i] = (fi[i] - f0) / delta

        return j

    @staticmethod
    def f_eqcons(ux, *args):
        """
        Equality constraints for adjoint variables.
        SLSQP tries to keep all elements of `discont` equal to 0.

        `ux` contains free inputs and adjoint states.
        """
        # Unpack inputs
        self, n_interv, free, n_free, model, slinp, shape_u, \
            shape_xadj, x0, n_states, unominal, ynominal = args

        # Extract `u` and `xadj` from `ux`
        size_u = shape_u[0] * shape_u[1]
        u = ux[:size_u].reshape(shape_u)
        u = u.reshape((n_interv, n_free))
        xadj = ux[size_u:].reshape(shape_xadj)

        # Denormalize inputs with unominal
        u = u * unominal

        # Initialize vectors for adjoint states
        xi = np.full((n_interv, n_states), np.nan)  # Initial states
        xf = np.full((n_interv, n_states), np.nan)  # Final states

        # Go through all intervals
        # Sequential simulations
        # ======================
        for i in range(n_interv):
            ydf, xdf = self.test_input_const(u=u[i],
                                                unames=free,
                                                model=model,
                                                inp=slinp[i],
                                                x0=xadj[i])
            # Extract initial and final state for interval `i`
            xi[i, :] = xdf.values[0, :]
            xf[i, :] = xdf.values[-1, :]

        # Calculate discontinuities
        discont = np.full((n_interv+1, n_states), np.nan)
        discont[0] = xi[0] - x0             # Initial state
        discont[-1] = xf[-1] - xadj[-1]     # Final state
        for i in range(1, n_interv):        # Stop iteration one before final
            discont[i] = xi[i] - xf[i-1]    # Adjoint states

        # Return flattened array
        return discont.flatten()
