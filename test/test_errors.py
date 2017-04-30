import random

from test.configurable import Configurable
from impact.errors import RealPartErrorEvaluator
from impact.constants import k_fm, k_norm

import numpy as np
import unittest

class TestRealErrors(Configurable):

	def setUp(self):
		super(TestRealErrors, self).setUp()
		self.parameters = self.data['initial_parameters'] + [self.data['SIGMA'], self.data['RHO']]
		uncertanities = self.data['DSIGMA'], self.data['DRHO']

		self.real_impact = self.data['real_impact_amplitude_error']
		cov_size = self.data['pdimension']

		# Create some fake covariance matrix
		cov = [[ 1./(i ** 2 + j + 1) 
			for i in range(cov_size)] for j in range(cov_size)]

		self.evaluator = RealPartErrorEvaluator(cov, *uncertanities)
		self.analytic_amplitude = self.analytic_formula()

	def npoints(self):
		return np.linspace(1e-5, 3, 100)


	def analytic_formula(self):
		try:
			import sympy 
			from sympy import exp as Exp
			from sympy import sqrt, pi
			allpars = "t a1 a2 a4 b1 b2 b3 b4 a5 b5 b6 a_s rho"
			t, a1, a2, a4, b1, b2, b3, b4, a5, b5, b6, a_s, rho =  sympy.symbols(allpars)
			a_s = a_s / (sqrt(pi *  k_norm) * 4)
			alpha = (1 - 1j*rho)*(a_s + a4)
			return 1j*alpha*( a1*Exp(-0.5*alpha*b1*t) + (1 - a1)*Exp(-0.5*alpha*b2*t) ) - 1j*a4*Exp(-0.5*b4*t) - a4*rho/((1 + t/b5)**4)

		except ImportError:
			return None


	def testRealImpactError(self):
		self.real_gamma_error = lambda x: self.evaluator.breal_error(x, self.parameters)
		data = map(self.real_gamma_error, self.npoints())

		for a, b in zip(data, self.real_impact):
				self.assertAlmostEqual(a, b)


	@unittest.skip('t')
	def testParameterA1Contribution(self):
		b1 = [i for i in self.analytic_amplitude.free_symbols if i.name == 'b1'][0]
		partial_derivative = self.analytic_amplitude.diff(b1)
		# print partial_derivative
		# TODO: Compare this partial derivative to evaluator.b1 function 


