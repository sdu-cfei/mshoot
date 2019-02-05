import numpy as np
import scipy.optimize

def optimize(func, x0, jac, bounds, f_eqcons, args, iter, solver='SLSQP'):
    """
    Interface to various optimization solvers. The API similar to
    scipy.optimize.fmin_slsqp, which initially was the only available
    solver.

    Supported solvers: 'SLSQP'

    The equality constraints `f_eqcons` must be given as a callable function
    f(x,*args), returning a 1-D array in which each element must equal 0.0
    in a successfully optimized problem. If f_eqcons is specified, eqcons is ignored.

    Return:
    out : ndarray of float, the final minimizer of func.
    fx : ndarray of float, if full_output is true, the final value of the objective function.
    its : int, if full_output is true, the number of iterations.
    imode : int, if full_output is true, the exit mode from the optimizer (see below).
    smode : string, if full_output is true

    :param func: callable f(x,*args), objective function
    :param x0: 1-D ndarray of float, initial guess for the independent variable(s)
    :param jac: callable f(x, *args), obj. function Jacobian
    :param bounds: A list of tuples specifying the lower and upper bound for each independent variable [(xl0, xu0),(xl1, xu1),...]
    :param f_eqcons: callable f(x,*args), returns a 1-D array in which each element must equal 0.0
    :param args: sequence, additional arguments passed to func and f_eqcons
    :param iter: int, maximum number of iterations
    :param solver: str, solver name
    :return: out, fx, its, imode, smode (description above)
    """

    # SLSQP (dedicated interface) ===================================
    if solver == 'fmin_slsqp':
        out, fx, its, imode, smode = scipy.optimize.fmin_slsqp(
            func=func,
            x0=x0,
            bounds=bounds,
            fprime=jac,
            f_eqcons=f_eqcons,
            args=args,
            iprint=2,
            acc=1e-6,
            iter=iter,
            full_output=True,
            epsilon=1e-9
        )

    # SLSQP =========================================================
    elif solver == 'SLSQP':

        res = scipy.optimize.minimize(
            fun=func,
            jac=jac,
            x0=x0,
            args=args,
            method=solver,
            bounds=bounds,
            constraints={'type': 'eq', 'fun': f_eqcons, 'args': args},
            options={'iprint': 2, 'disp': True, 'maxiter': iter}
        )

        # Unpack optimization results
        out = res.x
        fx = res.fun
        its = res.nit
        imode = res.status
        smode = res.message

        return out, fx, its, imode, smode

    else:
        print("Unknown solver")  # TODO: log
        raise KeyError("Unknown solver")

    # Return ========================================================
    return out, fx, its, imode, smode