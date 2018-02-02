from test.configurable import Configurable
from impact.datapoint import DataSet
from impact.datafit import DataFit
import impact.model as model

class TestImpactAmplitude(Configurable):

    def setUp(self):
        super(TestImpactAmplitude, self).setUp()
        self.nominal_value = self.data['impact_amplitude']

    def testValues(self):
        data = DataSet(self.infile, self.data).data
        gamma_fitter = DataFit(self.PROCESS + str(self.ENERGY), self.PROCESS, self.ENERGY, self.SIGMA, self.RHO) 
        parameters, covariance = gamma_fitter.fit(data)

        result = model.approx.values(data, parameters)
        for a, b in zip(result, self.nominal_value):
                self.assertAlmostEqual(a, b)
