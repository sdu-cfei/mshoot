import unittest
from test.test_mshoot import TestMShoot
from test.test_fmi import TestFMI
from test.test_mpc import TestMPC
import logging

# Set up logging
logging.basicConfig(filename='test.log', filemode='w', level='DEBUG')

def suite():
    suite = unittest.TestSuite()
    suite.addTest(TestMShoot('test_optimize_const_bounds'))
    suite.addTest(TestMShoot('test_optimize_const_bounds_joined'))
    suite.addTest(TestMShoot('test_optimize_step_bounds'))
    suite.addTest(TestMShoot('test_infeasible_x0'))
    suite.addTest(TestMShoot('test_normalized_uy'))
    suite.addTest(TestFMI('test_simulation'))
    suite.addTest(TestFMI('test_optimal_control'))
    suite.addTest(TestMPC('test_mpc'))
    suite.addTest(TestMPC('test_mpc_inp_clb'))
    return suite


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    runner.run(suite())