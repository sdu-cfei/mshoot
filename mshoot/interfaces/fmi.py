import logging
import os
import sys
from collections import OrderedDict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyfmi

from mshoot import SimModel


class SimFMU(SimModel):

    def __init__(self, fmupath, outputs=[], states=[], parameters={},
                 verbose=False):
        """
        :param fmupath: str, path to FMU
        :param outputs: list(str), monitored outputs names
        :param states: list(str), monitored states names
        :param parameters: dict, parameters names and values
        :param verbose: bool, whether to suppress pyfmi prints
        """
        print(fmupath)
        self.fmu = pyfmi.load_fmu(fmupath)
        self.outputs = outputs
        self.states = states
        self.verbose = verbose

        # Get initial state
        # Comment:
        #   The model has to be initialized to read the state variables.
        #   The easiest way to initialize is to run a short simulation.
        dummy_result = self.fmu.simulate(start_time=0, final_time=1)
        self.x0 = self._get_state()

        # Reset the FMU
        self.fmu.reset()

        # Set parameters
        for n in parameters:
            self.fmu.set(n, parameters[n])

    def _get_state(self):
        """
        Return an ordered dictionary with state names as keys
        and state values as values.
        """
        # Return dictionary, keys - state names, values - state values
        x = OrderedDict()

        # Get dictionary, keys - state names,
        # values - ScalarViariable instances (svi)
        x_svi = self.fmu.get_states_list()

        for s in x_svi:
            vr = x_svi[s]._get_value_reference()
            x[s] = self.fmu.get_real(vr)[0]  # [0] because 1-element array

        return x

    def simulate(self, udf, x0, save_state=False):
        """
        Simulate the model using the provided inputs `udf`
        and initial state `x0`.

        The DataFrame should have the following content:
        - index - time in seconds and equal steps, named 'time',
        - columns - input data,
        - column names - input variable names.

        The order of `x0` should reflect the one used in `states`.

        Return two DataFrames, `ydf` and `xdf`, with
        outputs and states, respectively, and with the same
        structure as `udf`.

        :param udf: DataFrame, shape (n_steps, n_variables)
        :param x0: vector, size (n_states, )
        :return: ydf, xdf
        """
        assert udf.index.name == 'time'

        start = udf.index[0]  # Start time
        stop = udf.index[-1]  # Final time
        ncp = udf.index.size - 1  # Number of communication points

        # Prepare inputs for pyfmi:
        # From pyfmi documentation:
        # "The input should be a 2-tuple consisting of first the names
        # of the input variable(s) and then the data matrix"
        # df = udf.reset_index()
        inp = (udf.columns, udf.reset_index().values)

        # FMI options
        opts = self.fmu.simulate_options()
        opts['ncp'] = ncp
        opts['result_handling'] = 'memory'              # Prevents saving result file
        opts['result_handler'] = 'ResultHandlerMemory'  # Prevents saving result file
        # if 'solver' in opts:
        #     # Model Exchange
        #     opts['solver'] = 'CVode'
        #     opts['CVode_options'] = {'rtol': 1e-6, 'atol': 1e-6}

        # Initial states from previous FMU simulation
        for n in self.x0:
            self.fmu.set(n, self.x0[n])

        # Initial states overriden by the user
        i = 0
        for n in self.states:
            self.fmu.set(n, x0[i])
            i += 1

        # Simulate
        if not self.verbose:
            nullf = open(os.devnull, 'w')
            sys.stdout = nullf

        res = self.fmu.simulate(start_time=start, final_time=stop,
                                input=inp, options=opts)

        if not self.verbose:
            sys.stdout = sys.__stdout__
            nullf.close()

        # Update state (use only in emulation)
        if save_state:
            self.x0 = self._get_state()

        # Outputs
        t = res['time']

        ydf = pd.DataFrame(index=pd.Index(t, name='time'))
        xdf = pd.DataFrame(index=pd.Index(t, name='time'))

        for n in self.outputs:
            ydf[n] = res[n]

        for n in self.states:
            xdf[n] = res[n]

        # Align time with udf
        # BUG: round-off errors for large numbers! is this code even needed?
        # ydf = ydf.loc[[i for i in t if i in udf.index]]
        # xdf = xdf.loc[[i for i in t if i in udf.index]]

        # Reset (note: won't work with E+)
        self.fmu.reset()

        return ydf, xdf


if __name__ == "__main__":
    # DEMO: SIMULATE
    # ==============
    # Load FMU
    fmupath = os.path.join('resources', 'fmus', 'R2C2', 'R2C2.fmu')
    parameters = {'C': 1e6}
    model = SimFMU(
        fmupath,
        outputs=['qout', 'Tr'],
        states=['heatCapacitor1.T'],
        parameters=parameters,
        verbose=True)

    # Inputs
    t = np.arange(0, 86401, 3600)
    udf = pd.DataFrame(index=pd.Index(t, name='time'), columns=['q', 'Tout'])
    udf['q'] = np.full(t.size, 100)
    udf['Tout'] = np.full(t.size, 273.15)

    # Initial state
    x0 = [273.15 + 20]

    ydf, xdf = model.simulate(udf, x0)

    ydf.plot(subplots=True, title='ydf')
    xdf.plot(subplots=True, title='xdf')
    plt.show()
