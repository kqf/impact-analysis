import unittest

import numpy as np

from impact.parametrization.numeric import Numeric
from impact.parametrization.symbolic import Standard
from test.configurable import Configurable


class TestNumericSymbolicConsistency(Configurable):

    def setUp(self):
        super(TestNumericSymbolicConsistency, self).setUp()
        data = self.data["test_real_errors"]
        self.params = data["params"] + [self.dataset.sigma, self.dataset.rho]

    @unittest.skip("Check the normalisation in numeric amplitude")
    def test_sympy_calculates_partial_derivatives(self):
        es, ep = Standard(), Numeric()
        numeric = (ep.d_a1, ep.d_a2, ep.d_b1, ep.d_b2,
                   ep.d_b3, ep.d_b4, ep.d_as, ep.d_rho)
        standard = (es.d_a1, es.d_a2, es.d_b1, es.d_b2,
                    es.d_b3, es.d_b4, es.d_as, es.d_rho)

        for standard, numeric in zip(standard, numeric):
            # print standard
            for t in np.linspace(0.1, 10):
                symval = standard(t, self.params)
                trueval = numeric(t, self.params)
                self.assertAlmostEqual(symval, trueval)
