"""
This script shows how to fit a scikit model to a state-space model.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logging
import mshoot
from mshoot.interfaces.scikit import SimScikit
from examples.models import statespace
from resources.ou44 import ou44_data

# Set up logging
logging.basicConfig(filename='test.log', filemode='w', level='DEBUG')

# Get "ground truth" data using the state space model
# ===================================================

# State space model
class SSModel(mshoot.mshoot.SimModel):
    def __init__(self, ss):
        """
        Constructor is optional.
        """
        self.ss = ss

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

# Read inputs: period [h], dt [min]
data = ou44_data.get_inputs(period=240, dt=60)
data = data.reset_index()
data.index = (data['datetime'] - data['datetime'][0]).dt.total_seconds()

udf_ss = data[['Tout', 'Hglo', 'qhvac', 'nocc']]

# Instantiate control and emulation models
model_ss = SSModel(statespace.get_model())

ydf_ss, xdf_ss = model_ss.simulate(udf_ss, np.array([20., 20., 20.]))

# Train scikit models
# ===================
train_end = 86400 * 7
utrain = udf_ss.loc[:train_end]
xtrain = xdf_ss.loc[:train_end]
ytrain = ydf_ss.loc[:train_end]

# Choose ML model:
ml = SimScikit(model="Linear", use_state=True)  # <<< WORKS BEST
# ml = SimScikit(model="Ridge", use_state=True)
# ml = SimScikit(model="NearestNeighbors", use_state=True, n_neighbors=5, weights='distance')
# ml = SimScikit(model="RandomForest", use_state=True)
# ml = SimScikit(model="SVM", use_state=True, kernel='linear')
# ml = SimScikit(model="NeuralNetwork", use_state=True, hidden_layer_sizes=(20, ),
#     activation='logistic', shuffle=True, max_iter=5000, solver='lbfgs')  # << WORKS FINE


ml.train(utrain, xtrain, ytrain)

# Simulate scikit models
# ======================
ydf_ml, xdf_ml = ml.simulate(udf_ss, xdf_ss.iloc[0])

print(ydf_ml.head())
print(xdf_ml.head())


# Change index from [s] to [h]
# ============================
def s2h(*args):
    for df in args:
        df.index /= 3600.

s2h(xdf_ss, xdf_ml)

# Compare state space with linear regression
# ==========================================
ax = xdf_ss.plot(color=['red', 'blue', 'green'])
ylim = ax.get_ylim()
ax.set_ylim(ylim)
ax.plot(xdf_ml['Tr'], label='Tr(ML)', ls='-.', c='red')
ax.plot(xdf_ml['Ti'], label='Ti(ML)', ls='-.', c='blue')
ax.plot(xdf_ml['Te'], label='Te(ML)', ls='-.', c='green')
ax.vlines(x=train_end/3600., ymin=ylim[0], ymax=ylim[1], linestyles='dashed', color='k')
ax.set_xticks(np.arange(0, xdf_ss.index.max() + 1, 24))
ax.set_xlabel('Time [h]')
ax.set_ylabel(r'Temperature [$^\circ$C]')
ax.legend()
ax.text(x=train_end/3600.-3, y=ylim[1]-0.5, s='Training', horizontalalignment='right')
ax.text(x=train_end/3600.+3, y=ylim[1]-0.5, s='Validation', horizontalalignment='left')
ax.set_title('LTI State-Space Model vs. Machine Learning Model')

# udf_ss.plot(subplots=True)
plt.show()
