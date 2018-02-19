import unittest

import cmath
import numpy as np
import sympy as smp


from test.configurable import Configurable
from impact.parametrization.numeric import Numeric
from impact.parametrization.symbolic import Symbolic



class TestNumericSymbolicConsistency(Configurable):

    def setUp(self):
        super(TestNumericSymbolicConsistency, self).setUp()
        self.parameters = self.data['initial_parameters'] + [self.dataset.sigma, self.dataset.rho]


    def test_sympy_calculates_partial_derivatives(self):
        es, ep = Symbolic(), Numeric()
        numeric  = ep.d_a1, ep.d_a2, ep.d_b1, ep.d_b2, ep.d_b3, ep.d_b4, ep.d_as, ep.d_rho
        symbolic = es.d_a1, es.d_a2, es.d_b1, es.d_b2, es.d_b3, es.d_b4, es.d_as, es.d_rho 

        for symbolic, numeric in zip(symbolic, numeric):
            for t in np.linspace(0.1, 10):
                symval = symbolic(t, self.parameters)
                trueval = numeric(t, self.parameters)
                self.assertAlmostEqual(symval, trueval) 