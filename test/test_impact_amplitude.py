from impact.ComputeGamma import ComputeGamma
from test.configurable import Configurable

class TestImpactAmplitude(Configurable):

	def setUp(self):
		super(TestImpactAmplitude, self).setUp()
		self.nominal_value = self.data['impact_amplitude']


	def testValues(self):
		c = ComputeGamma(self.infile, self.PROCESS, self.ENERGY, self.SIGMA, self.RHO) 
		result = c.compute()

		for a, b in zip(result, self.nominal_value):
				self.assertAlmostEqual(a, b)
