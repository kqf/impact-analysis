import random

from impact.constants import k_fm, k_norm
from test.configurable import Configurable
import impact.errors.real as err_real

from impact.parametrization.numeric import Numeric
from impact.parametrization.symbolic import Symbolic


import numpy as np
import cmath
import unittest
import sympy as smp

# from cmath import complex

class TestRealErrors(Configurable):

    def setUp(self):
        super(TestRealErrors, self).setUp()
        cov_size = self.data['pdimension']
        # Create some fake covariance matrix
        cov = [[ 1./(i ** 2 + j + 1) 
            for i in range(cov_size)] for j in range(cov_size)]
            
        self.dataset.covariance = cov
        self.dataset.parameters = self.data['initial_parameters'] + [self.dataset.sigma, self.dataset.rho]

    def tearDown(self):
        self.dataset.covariance = None
        self.dataset.parameters = None


    # @unittest.skip('This one should has the lowest priority,\
        # as there is an error in the formula, this error should be studied later')
    def test_real_part_of_impact_amplitude_error(self):
        evaluator = err_real.Error(Numeric())
        data = map(
            lambda x: evaluator.evaluate_(x, self.dataset), 
            np.linspace(1e-5, 3, 100)
        )

        self.dataset.covariance = None
        self.dataset.parameters = None 

        nominal_real_impact_amplitude_error = self.data['real_impact_amplitude_error']
        mymsg = 'Actual values differ from nominal estimates.\n\nActual values: {}'.format(data)
        for a, b in zip(data, nominal_real_impact_amplitude_error):
                self.assertAlmostEqual(a, b, msg = mymsg)



class TestNumericSymbolicConsistency(Configurable):

    def setUp(self):
        super(TestNumericSymbolicConsistency, self).setUp()
        self.parameters = self.data['initial_parameters'] + [self.dataset.sigma, self.dataset.rho]


        # TODO: Check extra k_norm factor
    def test_sympy_calculates_partial_derivatives(self):
        es, ep = Symbolic(), Numeric()
        numeric  = ep.d_a1, ep.d_a2, ep.d_b1, ep.d_b2, ep.d_b4, ep.d_as, ep.d_rho
        symbolic = es.d_a1, es.d_a2, es.d_b1, es.d_b2, es.d_b4, es.d_as, es.d_rho # Why there is no b3?

        for symbolic, numeric in zip(symbolic, numeric):

            for t in np.linspace(0.1, 10):
                symval = symbolic(t, self.parameters)
                trueval = numeric(t, self.parameters)
                self.assertAlmostEqual(symval, trueval) 