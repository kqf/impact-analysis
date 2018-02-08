import unittest
from test.configurable import Configurable
from impact.utils import hankel_transform
from impact.parametrization.numeric import Numeric
import numpy as np

class TestHankelTransformation(Configurable):

	def setUp(self):
		super(TestHankelTransformation, self).setUp()
		self.parameters = self.data['initial_parameters'] + [self.dataset.sigma, self.dataset.rho]
		self.cov_size = 6


	def npoints(self):
		return np.linspace(1e-5, 3, 100)


	def testRealAmplitude(self):
		from scipy import integrate
		from scipy.special import j0
		from math import sqrt, pi
		from impact.constants import k_fm, k_norm

		model = Numeric()

		def real_gamma_explicit_form(b, p):
		    f = lambda q :  q * j0(b * q / k_fm) * model.amplitude(q * q, p).real / sqrt(pi * k_norm)
		    result = integrate.quad(f, 0, np.infty)[0]
		    return -result

		f1 = lambda x: model.real_gamma(x, self.parameters)
		f2 = lambda x: real_gamma_explicit_form(x, self.parameters)

		for b in self.npoints():
				self.assertAlmostEqual(f1(b), f2(b))