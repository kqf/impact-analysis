from impact.ComputeGamma import ComputeGamma
from test.configurable import Configurable

class TestImpactAmplitude(Configurable):

	def setUp(self):
		super(TestImpactAmplitude, self).setUp()
		self.nominal_value = self.data['impact_amplitude']


	def testValues(self):
		# TODO: Try to add quiet/dead mode
		c = ComputeGamma(self.infile, self.PROCESS, self.ENERGY, self.SIGMA, self.RHO) 
		result = c.performComputations()

		for a, b in zip(result, self.nominal_value):
				self.assertAlmostEqual(a, b)
