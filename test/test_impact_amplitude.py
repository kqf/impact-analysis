from test.configurable import Configurable
from impact.datapoint import DataPoint
from impact.datafit import DataFit
import impact.model as model

class TestImpactAmplitude(Configurable):

    def setUp(self):
        super(TestImpactAmplitude, self).setUp()
        self.nominal_value = self.data['impact_amplitude']

    def testValues(self):
        data = DataPoint.read(self.ENERGY, self.PROCESS, self.infile)
        gamma_fitter = DataFit(data, self.PROCESS + str(self.ENERGY), self.PROCESS, self.ENERGY, self.SIGMA, self.RHO) 
        parameters, covariance = gamma_fitter.fit()

        result = model.approx.values(data, parameters)
        for a, b in zip(result, self.nominal_value):
                self.assertAlmostEqual(a, b)
