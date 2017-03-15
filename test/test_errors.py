import random
from impact.Formulas import getRealGammaError
from test.configurable import Configurable
import numpy as np

class TestRealErrors(Configurable):

	def setUp(self):
		super(TestRealErrors, self).setUp()
		self.parameters = self.data['initial_parameters'] + [self.data['SIGMA'], self.data['RHO']]
		self.real_impact = self.data['real_impact_amplitude_error']
		# TODO: Move it to the config file
		self.cov_size = 6

	def npoints(self):
		return np.linspace(1e-5, 3, 100)


	def testRealImpactError(self):
		cov = [[ 1./(i ** 2 + j + 1) 
            for i in range(self.cov_size)] for j in range(self.cov_size)]

		data = [getRealGammaError([b], self.parameters, cov, *self.parameters[-2:]) for b in self.npoints()]

		for a, b in zip(data, self.real_impact):
				self.assertAlmostEqual(a, b)
