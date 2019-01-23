# MShoot - Multiple Shooting Model Predictive Control

MShoot works with both physical and data-driven models. Models can be connected through one of the three available interfaces: the Functional Mock-up Interface, the scikit-learn interface, and the generic Python interface.

![Architecture](/examples/bs2019/figs/architecture.png)

The control model is used inside the MPC loop (figure above) for the control strategy optimization for each optimization horizon. The optimization is performed using the multiple shooting method, in which the optimization horizon is divided into N subintervals (figure below). The state trajectories within each subinterval are obtained through simulation. The states are free to vary within the feasible region (hard constraints), however if the initial state violates the constraints, the constraints are temporarily relaxed in the affected subintervals so that the solver can proceed with the calculations. The NLP solver solves the constrained minimization problem augmented by the state continuity constraints.

![Multiple shooting](/examples/bs2019/figs/multiple_shooting.png)

The documentation is coming soon!
