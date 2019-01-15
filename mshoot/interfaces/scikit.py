import logging
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sklearn
import sklearn.linear_model
import sklearn.neighbors
import sklearn.ensemble
import sklearn.multioutput
import sklearn.neural_network
import sklearn.svm
from sklearn.preprocessing import StandardScaler

from mshoot import SimModel

class SimScikit(SimModel):
    """
    Unlike other interface classes (e.g. SimFMU) which are simple wrappers
    for the existing models, SimScikit works more like an automation layer
    for Scikit-Learn models, covering model instantiation, training,
    and prediction (simulation). In addition, the selected ML model is
    turned into a dynamic model by:
    - predicting the time derivative of the state,
    - running a model in a time loop, at each step predicting one-step ahead,
    - using the time-backshifted output as one of the predictors.

    The way this class should be used is as follows:
    1. Select the model type with __init__(...)
    2. Train the model with train(...)
    3. Simulate with simulate(...)

    TODO: Linear interpolation between time steps -> 0.5 * (feat[i-1] + feat[i])
          (currently x[i] is calculated based on x[i-1], dt[i], and feat[i] only)

    :param model: str, Scikit-Learn regression model type
    :param use_state: bool, If True, state is included in the features
    :param **args: Optional arguments to be passed to a sklearn model
    """
    def __init__(self, model="Linear", use_state=False, **args):
        self.labels = list()     # Model output names
        self.states = list()     # State names
        self.outputs = list()    # Output names
        self.model_type = model  # Name of the model from sklearn
        self.model = None        # Actual model instance
        self.use_state = use_state
        self.scaler_feat = StandardScaler()
        self.scaler_labs = StandardScaler()

        available_models = ['Linear', 'Ridge', 'NearestNeighbors', 'RandomForest',
            'NeuralNetwork', 'SVM']

        if model not in available_models:
            msg = "Unsupported model type: {}\n".format(self.model_type)
            msg += "Available models: {}".format(available_models)
            logging.error(msg)
            raise ValueError(msg)

        if self.model_type == "Linear":
            self.model = sklearn.linear_model.LinearRegression(**args)
        elif self.model_type == "Ridge":
            self.model = sklearn.linear_model.Ridge(**args)
        elif self.model_type == "NearestNeighbors":
            self.model = sklearn.neighbors.KNeighborsRegressor(**args)
        elif self.model_type == "RandomForest":
            rf = sklearn.ensemble.RandomForestRegressor(**args)
            self.model = sklearn.multioutput.MultiOutputRegressor(estimator=rf)
        elif self.model_type == "NeuralNetwork":
            self.model = sklearn.neural_network.MLPRegressor(**args)
        elif self.model_type == "SVM":
            svr = sklearn.svm.SVR(**args)
            self.model = sklearn.multioutput.MultiOutputRegressor(estimator=svr)
        else:
            assert False, "This line should never be reached"

    def train(self, udf, xdf, ydf):
        """
        Train the model. The model can be retrained any time
        between the simulations.

        :param udf: DataFrame, shape (n_steps, n_variables)
        :param xdf: DataFrame, shape (n_steps, n_states)
        :param ydf: DataFrame, shape (n_steps, n_outputs)
        :return: None
        """
        assert isinstance(udf, pd.DataFrame)
        assert isinstance(xdf, pd.DataFrame)

        # Features
        # ===============================
        # Extract features and time steps
        feat = udf.values
        t = udf.index.values
        dt = t - np.roll(t, 1)
        dt[0] = t[1] - t[0]
        dt = dt.reshape(-1, 1)

        # Add previous states to the features
        x = xdf.values
        if self.use_state:
            xprev = np.roll(x, 1, axis=0)
            xprev[0] = np.nan
            feat = np.hstack((feat, xprev))

        # Labels
        # ===============================
        # These labels are used in ydf, xdf returned by simulate()
        # Note that states and labels may have the same names
        states = [s + '(sta)' for s in xdf.columns]
        outputs = [s + '(out)' for s in ydf.columns]
        self.labels = states + outputs
        self.states = states
        self.outputs = outputs

        # Extract state changes
        labs = x - np.roll(x, 1, axis=0)
        labs[0] = x[1] - x[0]

        # Add outputs
        labs = np.hstack((labs, ydf.values))

        # Training
        # ===============================
        # Exclude the first row of features, due to unavailable previous state
        feat = feat[1:, :]
        labs = labs[1:, :]

        # Scale data
        self.scaler_feat.fit(feat)
        self.scaler_labs.fit(labs)

        # Train
        self.model.fit(self.scaler_feat.transform(feat), self.scaler_labs.transform(labs))

    def simulate(self, udf, x0, **kwargs):
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

        :param udf: DataFrame, shape (n_variables, n_steps)
        :param x0: vector, size (n_states, )
        :return: ydf, xdf
        """
        # Features
        # ===============================
        # Extract features and time steps
        feat = udf.values
        t = udf.index.values
        dt = t - np.roll(t, 1)
        dt[0] = t[1] - t[0]
        dt = dt.reshape(-1, 1)

        # Simulation
        # ===============================
        x = np.full((len(dt), len(self.states) + len(self.outputs)), np.nan)
        x[0, :len(self.states)] = x0
        x[0, len(self.states):] = 0
        for i in range(1, len(dt)):
            # Get current observation
            feat_i = feat[i:i+1]

            # Add previous state to features
            if self.use_state:
                feat_i = np.hstack((feat_i, x[i-1:i, :len(self.states)]))

            # Scale features
            feat_i = self.scaler_feat.transform(feat_i)

            # Predict next state
            y = self.model.predict(feat_i)
            y = self.scaler_labs.inverse_transform(y)  # Rescale outputs
            x[i, :len(self.states)] = x[i-1, :len(self.states)] + y[0, :len(self.states)]
            x[i, len(self.states):] = y[0, len(self.states):]

        # Pack results into data frames
        # ===============================
        xydf = pd.DataFrame(index=udf.index, data=x, columns=self.labels)
        ydf = xydf[self.outputs]
        xdf = xydf[self.states]

        # Rollback to original names (remove '(sta)' and '(out)')
        ydf = ydf.rename(columns = {s:s[:-5] for s in self.outputs})
        xdf = xdf.rename(columns = {s:s[:-5] for s in self.states})

        return ydf, xdf


if __name__ == "__main__":
    # Example:
    # ======================================
    # Dynamic model of the process with two outputs (y1, y2):
    def sim_y(x1, x2, y0, dt):
        """
        :param x1: vector, input 1
        :param x2: vector, input 2
        :param y0: vector, initial y
        :param dt: float, time step size
        :return: output vector
        """
        assert len(x1) == len(x2)

        y = np.full((x1.size, 2), np.nan)
        y[0] = y0

        for i in range(1, len(x1)):
            y[i, 0] = y[i-1, 0] + (x1[i] + x2[i]) * dt
            y[i, 1] = y[i-1, 1] + (2 * x1[i] - x2[i]) * dt

        return y

    # Time step
    dt = 1.

    # Timeline
    t = np.arange(0, 100, dt)

    # Inputs
    x1 = np.sin(t * 4 * np.pi / max(t))
    x2 = t / 60.

    # Actual (ground truth) output
    y_actual = sim_y(x1, x2, y0=[0., 1.], dt=dt)

    # Noisy output
    y_noisy = y_actual + (np.random.rand(*y_actual.shape) - 0.5)

    # Labels for training of the ML model
    labels = pd.DataFrame(index=t, data=y_noisy, columns=['y1', 'y2'])

    # Inputs for training and simulation of the ML model
    udf = pd.DataFrame(
        data = np.hstack([x1.reshape(-1, 1), x2.reshape(-1, 1)]),
        index = pd.Index(t, name='time'),
        columns=['x1', 'x2']
    )

    # Use the scikit interface
    simsci = SimScikit(model='Linear', use_state=True)
    simsci.train(udf, labels)
    ydf, xdf = simsci.simulate(udf, y_noisy[0])

    # Print and plot the results
    plt.plot(t, y_actual, label="Actual", ls=':')
    plt.plot(t, y_noisy, label="Noisy", ls='--')
    plt.plot(t, ydf, label="Predicted")
    plt.legend()
    plt.show()
