from test.configurable import Configurable
from impact.model import GammaApproximation
from impact.datapoint import DataReader

# TODO: try to avoid numpy ?!
import numpy as np

class TestGammaApproximation(Configurable):

    def setUp(self):
        super(TestGammaApproximation, self).setUp()
        self.low_t_extrapolation = self.data['low_t_extrapolation']
        self.parameters = self.data['initial_parameters'] + [self.SIGMA, self.RHO]


    def testLowTExrapolation(self):
        data = DataReader(self.ENERGY, self.PROCESS).read(self.infile)
        approximator = GammaApproximation(data)

        values = [approximator.im_amplitude_low_t(t, self.parameters) for t in np.linspace(0, 0.5, 100)]

        for a, b in zip(values, self.low_t_extrapolation):
                self.assertAlmostEqual(a, b)
