import random
from impact.impactanalysis import ImpactAnalysis
from test.configurable import Configurable

class TestMCAveragesOFImpactAmplitude(Configurable):

	def setUp(self):
		super(TestMCAveragesOFImpactAmplitude, self).setUp()
		random.seed(1234)
		self.nom_values, self.nom_errors = self.data['mc_averages']


	def testValues(self):
		analysis = ImpactAnalysis(self.infile, self.PROCESS, self.ENERGY, self.SIGMA, self.RHO, self.DSIGMA, self.DRHO, 100)
		values, errors = analysis.run()

		# Test values
		for a, b in zip(values, self.nom_values):
				self.assertAlmostEqual(a, b)

		# Test errors
		for a, b in zip(errors, self.nom_errors):
				self.assertAlmostEqual(a, b)
