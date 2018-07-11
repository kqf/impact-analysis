import numpy as np
from scipy import integrate
from scipy.special import j0
from math import pi
from impact.constants import k_fm

from impact.parametrization.numeric import Numeric
# from impact.utils import hankel_transform
from test.configurable import Configurable


class TestHankelTransformation(Configurable):

    def setUp(self):
        super(TestHankelTransformation, self).setUp()
        data = self.data["test_hankel_transform"]
        self.parameters = data["params"] + \
            [self.dataset.sigma, self.dataset.rho]
        self.cov_size = 6

    def npoints(self):
        return np.linspace(1e-5, 3, 100)

    def testRealAmplitude(self):
        model = Numeric()

        def real_gamma_explicit_form(b, p):
            def f(q):
                return q * j0(b * q / k_fm) * \
                    model.amplitude(q * q, p).real / 8 / pi
            result = integrate.quad(f, 0, np.infty)[0]
            return -result

        def f1(x):
            return model.real_gamma(x, self.parameters)

        def f2(x):
            return real_gamma_explicit_form(x, self.parameters)

        for b in self.npoints():
            self.assertEqual(f1(b), f2(b))
