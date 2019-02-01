# Examples

The directory contains exemplary scripts showing how to set up a simulation using MShoot.

Most of the scripts are named using the convention `type_interface_x`, where:

- `type` can be `mpc` or `mshoot`,
- `interface` can be `generic`, `fmi`, or `scikit`, 
- `x` is the example number (different settings are used in each).

The `mpc` examples present how to use the class `MPCEmulation` for the MPC loops (consisting of multiple optimization periods).
The `mshoot` examples present how to set up and solve single optimization period using multiple shooting.

The script `train_scikit.py` shows how to train and use a machine learning model (from scikit-learn)
using the MShoot's scikit interface.

Finally, the directory `bs2019` contains scripts used to produce results for the paper presenting MShoot and submitted
to Building Simulation 2019 (Rome, September 2019).
