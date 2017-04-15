from test.configurable import Configurable
from impact.ComputeGamma import ComputeGamma
from impact.model import GammaApproximation

# TODO: try to avoid numpy ?!
import numpy as np

class TestGammaApproximation(Configurable):

	def setUp(self):
		super(TestGammaApproximation, self).setUp()
		self.low_t_extrapolation = self.data['low_t_extrapolation']
		self.parameters = self.data['initial_parameters'] + [self.SIGMA, self.RHO]


	def testLowTExrapolation(self):
		c = ComputeGamma(self.infile, self.PROCESS, self.ENERGY, self.SIGMA, self.RHO) 
		approximator = GammaApproximation(c.dataPoints)

		values = [approximator.im_amplitude_low_t(t, self.parameters) for t in np.linspace(0, 0.5, 100)]

		for a, b in zip(values, self.low_t_extrapolation):
				self.assertAlmostEqual(a, b)
