.. mshoot documentation master file, created by
   sphinx-quickstart on Tue Feb  5 18:28:42 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome
=======

`mshoot` is a Python package for Model Predictive Control (MPC).

The package is based on the multiple shooting method, which is a simulation-based
dynamic optimization method. In multiple shooting the optimization period
is divided into N subperiods, each one simulated separately. The optimization
is solved using nonlinear programing (NLP). The optimization problem consists
of a cost function, which is calculated separately for each subperiod and then summed,
and state continuity constraints (:math:`x_{N-1}^R = x_{N}^L`).

.. image:: gfx/multiple_shooting.png

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   api
   license



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
