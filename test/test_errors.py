import random

from test.configurable import Configurable
from impact.errors import RealPartErrorEvaluator
from impact.constants import k_fm, k_norm

import numpy as np
import unittest
# from cmath import complex

class TestRealErrors(Configurable):

	def setUp(self):
		super(TestRealErrors, self).setUp()
		self.parameters = self.data['initial_parameters'] + [self.data['SIGMA'], self.data['RHO']]
		uncertanities = self.data['DSIGMA'], self.data['DRHO']

		self.real_impact = self.data['real_impact_amplitude_error']
		cov_size = self.data['pdimension']

		self.variable_names = 'a1 a2 a4 b1 b2 b3 b4 a5 b5 b6 a_s rho'.split()
		self.variables = {v: p for v, p in zip(self.variable_names, self.parameters)}

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
			t, a1, a2, a4, b1, b2, b3, b4, a5, b5, b6, a_s, rho =  sympy.symbols(['t'] + self.variable_names)
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


	# @unittest.skip('This test fails due to problems with math')
	def testParameterA1Contribution(self):
		b1 = [i for i in self.analytic_amplitude.free_symbols if i.name == 'a1'][0]
		partial_derivative = self.analytic_amplitude.diff(b1)

		# TODO: These functions give different results. Check the formula/values by hand
		#      
		# TODO: Check extra k_norm factor

		print self.parameters
		for t in np.linspace(0.2, 10):
			self.variables['t'] = t
			# print t
			analytic = complex(partial_derivative.evalf(subs = self.variables))
			# print analytic
			trueval = self.evaluator.d_a1(t, self.parameters)
			# print trueval


			self.assertAlmostEqual(analytic.real , trueval.real)
			# self.assertAlmostEqual(analytic.imag , trueval.imag)



