import unittest
from test.configurable import Configurable
from impact.Formulas import amplitude, hankel_transform, getRealError, getRealGammaError, real_gamma
import numpy as np

class TestHankelTransformation(Configurable):

	def setUp(self):
		super(TestHankelTransformation, self).setUp()
		self.parameters = self.data['initial_parameters'] + [self.data['SIGMA'], self.data['RHO']]
		self.cov_size = 6


	def npoints(self):
		return np.linspace(1e-5, 3, 100)


	def testRealAmplitude(self):
		from scipy import integrate
		from scipy.special import j0
		from impact.Formulas import k_fm, k_norm
		from math import sqrt, pi

		def real_gamma_explicit_form(b, p):
		    f = lambda q :  q * j0(b * q / k_fm) * amplitude(q * q, p).real / sqrt(pi * k_norm)
		    result = integrate.quad(f, 0, np.infty)[0]
		    return -result

		f1 = lambda x: real_gamma(x, self.parameters)
		f2 = lambda x: real_gamma_explicit_form(x, self.parameters)

		for b in self.npoints():
				self.assertAlmostEqual(f1(b), f2(b))

	# TODO: Implement uniform algorithm for all methods
	@unittest.skip('This test will fail as hankel_transform is not fully implemented.')
	def testRealGammaError(self):
		cov = [[ 1./(i ** 2 + j + 1) for i in range(self.cov_size)] for j in range(self.cov_size)]

		@hankel_transform
		def real_gamma_error(x, p, covariance, dsigma, drho):
			return getRealError(x, p, covariance, dsigma, drho)

		# TODO: do we need to pass sigma and rho explicitly ?
		nominal = [getRealGammaError(b, self.parameters, cov, *self.parameters[-2:]) for b in self.npoints()]
		actual  =  [real_gamma_error(b, self.parameters, cov, *self.parameters[-2:]) for b in self.npoints()]

		for a, b in zip(nominal, actual):
				self.assertAlmostEqual(a, b)