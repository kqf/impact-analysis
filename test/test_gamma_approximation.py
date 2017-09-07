from test.configurable import Configurable
import impact.model as model
from impact.datapoint import DataPoint

# TODO: try to avoid numpy ?!
import numpy as np

class TestGammaApproximation(Configurable):

    def setUp(self):
        super(TestGammaApproximation, self).setUp()
        self.low_t_extrapolation = self.data['low_t_extrapolation']
        self.parameters = self.data['initial_parameters'] + [self.SIGMA, self.RHO]


    def testLowTExrapolation(self):
        data = DataPoint.read(self.ENERGY, self.PROCESS, self.infile)
        approximator = model.approx(data)

        values = [approximator.im_amplitude_low_t(t, self.parameters) for t in np.linspace(0, 0.5, 100)]

        for a, b in zip(values, self.low_t_extrapolation):
                self.assertAlmostEqual(a, b)
