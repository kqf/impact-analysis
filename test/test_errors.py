import random

from test.configurable import Configurable
from impact.errors import RealPartErrorEvaluator
from impact.constants import k_fm, k_norm

import numpy as np
import cmath
import unittest
import sympy as smp

# from cmath import complex

class TestRealErrors(Configurable):

	def setUp(self):
		super(TestRealErrors, self).setUp()
		self.parameters = self.data['initial_parameters'] + [self.data['SIGMA'], self.data['RHO']]
		uncertanities = self.data['DSIGMA'], self.data['DRHO']

		self.real_impact = self.data['real_impact_amplitude_error']
		cov_size = self.data['pdimension']

		self.variable_names = 'a1 a2 a4 b1 b2 b3 b4 a5 b5 b6 a_s rho'.split()
		self.variables = smp.symbols(self.variable_names)
		self.t = smp.Symbol('t')
		self.longMessage = True

		# self.variables = {v: p for v, p in zip(self.variable_names, self.parameters)}

		# Create some fake covariance matrix
		cov = [[ 1./(i ** 2 + j + 1) 
			for i in range(cov_size)] for j in range(cov_size)]

		self.evaluator = RealPartErrorEvaluator(cov, *uncertanities)
		self.analytic_amplitude = self.analytic_formula()



	def analytic_formula(self):
		a1, a2, a4, b1, b2, b3, b4, a5, b5, b6, a_s, rho =  self.variables
		t = self.t
		a_s = a_s / (smp.sqrt(smp.pi *  k_norm) * 4)
		alpha = (1 - 1j*rho)*(a_s + a4)
		return 1j*alpha*( a1*smp.exp(-0.5*alpha*b1*t) + (1 - a1)*smp.exp(-0.5*alpha*b2*t) ) - 1j*a4*smp.exp(-0.5*b4*t) - a4*rho/((1 + t/b5)**4)



	@unittest.skip('This one should has the lowest priority,\
		as there is an error in the formula, this error should be studied later')
	def testRealImpactError(self):
		rgamma_error = lambda x: self.evaluator.breal_error(x, self.parameters)
		data = map(rgamma_error, np.linspace(1e-5, 3, 100))

		mymsg = 'Actual values differ from nominal estimates.\n\nActual values: {}'.format(rgamma_error)
		for a, b in zip(data, self.real_impact):
				self.assertAlmostEqual(a, b, msg = mymsg)


	# TODO: Check extra k_norm factor
	def testCheckExplicitFormula(self):
		fvariables = 'a1', 'a4', 'b1', 'b2'
		fmethods = self.evaluator.d_a1, self.evaluator.d_a4, self.evaluator.d_b1, self.evaluator.d_b2

		for arg, method in zip(fvariables, fmethods):
			fpar = next((i for i in self.analytic_amplitude.free_symbols if i.name == arg), None)
			partial_derivative = smp.lambdify((self.t, self.variables), self.analytic_amplitude.diff(fpar), 'numpy')

			mymsg = '\nThere is an error in formulas for partial derivative of A(s, t) over {0}'.format(arg)
			for t in np.linspace(0.2, 10):
				analytic = complex(partial_derivative(t + 0.1, self.parameters))
				trueval = method(t, self.parameters)

				self.assertAlmostEqual(analytic.real , trueval, msg = mymsg)



