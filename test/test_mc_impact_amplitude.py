from test.configurable import Configurable
from impact.ComputeGamma import ComputeGamma
import random

class TestMCImpactAmplitude(Configurable):
	"""
		This test checks if the procedure of mc impact amplitude generation 
		remains the same. 
	"""

	def setUp(self):
		super(TestMCImpactAmplitude, self).setUp()
		self.nominal_value = self.data['mc_impact_amplitude']
		random.seed(1234)


	def testValues(self):
		# TODO: Try to add quiet/dead mode
		c = ComputeGamma(self.infile, self.PROCESS, self.ENERGY, self.SIGMA, self.RHO) 
		result = c.generate_mc_gamma(100, 1, self.DSIGMA)

		for a, b in zip(result, self.nominal_value):
				self.assertAlmostEqual(a, b)
