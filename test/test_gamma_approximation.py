from test.configurable import Configurable
from impact.model import Approx
from impact.datapoint import DataSet
from impact.parametrization.numeric import Numeric
from impact.model import Approx

# TODO: try to avoid numpy ?!
import numpy as np

class TestGammaApproximation(Configurable):

    def setUp(self):
        super(TestGammaApproximation, self).setUp()
        self.low_t_extrapolation = self.data['low_t_extrapolation']
        self.parameters = self.data['initial_parameters'] + [self.dataset.sigma, self.dataset.rho]


    def test_extrapolates_low_t_correctly(self):
        approximator = Approx(Numeric(), self.dataset.data)

        values = [approximator.im_amplitude_low_t(t, self.parameters) for t in np.linspace(0, 0.5, 100)]

        for a, b in zip(values, self.low_t_extrapolation):
                self.assertAlmostEqual(a, b)
