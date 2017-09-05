import random
from impact.impactanalysis import ImpactAnalysis
from test.configurable import Configurable

class TestFinalResult(Configurable):

	def setUp(self):
		super(TestFinalResult, self).setUp()
		random.seed(1234)
		#TODO: Urgent! Adde tests for MC averages
		self.nom_values, self.nom_errors = self.data['final_result']
		self.longMessage = True


	def testValues(self):
		analysis = ImpactAnalysis(self.infile, self.PROCESS, self.ENERGY, self.SIGMA, self.RHO, self.DSIGMA, self.DRHO)
		values, errors = analysis.run()

		# Test values
		msg = 'Actual values:\n {}'.format(values)
		for a, b in zip(values, self.nom_values):
				self.assertAlmostEqual(a, b, msg = msg)

		# Test errors
		msg = 'Actual errors:\n {}'.format(errors)
		for a, b in zip(errors, self.nom_errors):
				self.assertAlmostEqual(a, b, msg = msg)
