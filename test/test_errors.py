import random

from test.configurable import Configurable
from impact.errors import RealPartErrorEvaluator

import numpy as np

class TestRealErrors(Configurable):

	def setUp(self):
		super(TestRealErrors, self).setUp()
		parameters = self.data['initial_parameters'] + [self.data['SIGMA'], self.data['RHO']]
		self.real_impact = self.data['real_impact_amplitude_error']

		# TODO: Move it to the config file?
		cov_size = 6
		cov = [[ 1./(i ** 2 + j + 1) 
			for i in range(cov_size)] for j in range(cov_size)]

		# TODO: Use here dsigma, drho
		self.real_gamma_error = lambda x: RealPartErrorEvaluator(cov, *parameters[-2:]).breal_error(x, parameters)

	def npoints(self):
		return np.linspace(1e-5, 3, 100)


	def testRealImpactError(self):
		data = map(self.real_gamma_error, self.npoints())

		for a, b in zip(data, self.real_impact):
				self.assertAlmostEqual(a, b)
