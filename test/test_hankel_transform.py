import unittest
from impact.Formulas import getRealGamma, getImagGamma, amplitude, hankel_transform, getRealError, getRealGammaError
from test.configurable import Configurable
import numpy as np

class TestHankelTransformation(Configurable):

	def setUp(self):
		super(TestHankelTransformation, self).setUp()
		self.parameters = self.data['initial_parameters'] + [self.data['SIGMA'], self.data['RHO']]
		self.cov_size = 6


	def npoints(self):
		return np.linspace(1e-5, 3, 100)

	# TODO: Leave only one test case, for the real/imag part.
	#       Swap definitions: Use decorators in code and hardocded version here for test purposes.

	def testRealAmplitude(self):

		@hankel_transform
		def real_gamma(x, p):
			return -amplitude(x, p).real

		f1 = lambda x: real_gamma(x, self.parameters)
		f2 = lambda x: getRealGamma(x, self.parameters)

		for b in self.npoints():
				self.assertAlmostEqual(f1(b), f2(b))


	def testRealAmplitude(self):

		@hankel_transform
		def imag_gamma(x, p):
			return -amplitude(x, p).imag

		f1 = lambda x: imag_gamma(x, self.parameters)
		f2 = lambda x: getImagGamma(x, self.parameters)

		for b in self.npoints():
				self.assertAlmostEqual(f1(b), f2(b))


	# TODO: Implement uniform algorithm for all methods
	@unittest.skip('This test will fail as hankel_transform is not fully implemented.')
	def testRealGammaError(self):
		cov = [[ 1./(i ** 2 + j + 1) for i in range(self.cov_size)] for j in range(self.cov_size)]

		@hankel_transform
		def real_gamma_error(x, p, covariance, dsigma, drho):
			return getRealError(x, p, covariance, dsigma, drho)

		nominal = [getRealGammaError(b, self.parameters, cov, *self.parameters[-2:]) for b in self.npoints()]
		actual  =  [real_gamma_error(b, self.parameters, cov, *self.parameters[-2:]) for b in self.npoints()]

		for a, b in zip(nominal, actual):
				self.assertAlmostEqual(a, b)