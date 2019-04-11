import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

from pyfmi import load_fmu

# The example shows how to calculate directional derivatives of states
# with respect to inputs u1, u2. Analytical directional derivatives
# available in some FMUs are not yet exploited in mshoot :(
#
# Model: R3C3.mo
#
# Included FMUs: 
#   - R3C3_ME.fmu           - derivatives work
#   - R3C3_CS_CVODE.fmu     - derivatives work
#
# Not included FMUs:
#   - R3C3_CS_DYMOLA.fmu    - derivatives not available

model = load_fmu(os.path.join('examples', 'derivatives', 'R3C3_CS_CVODE.fmu'))
model.simulate(final_time=1.)

inputs_names = ['u1', 'u2']
inputs_ref = model.get_model_time_varying_value_references(inputs_names)[0]

states = model.get_states_list()
states_ref = [s.value_reference for s in states.values()]

derivatives = model.get_derivatives_list()
derivatives_ref = [d.value_reference for d in derivatives.values()]

v = np.array([-1000., 1000.])

dd = model.get_directional_derivative(inputs_ref, derivatives_ref, v)

print(states)
print(derivatives)
print(dd)
