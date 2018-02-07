from test.configurable import Configurable
from impact.datapoint import DataSet
from impact.datafit import DataFit
import impact.model as model

class TestImpactAmplitude(Configurable):

    def setUp(self):
        super(TestImpactAmplitude, self).setUp()
        self.nominal_value = self.data['impact_amplitude']

    def testValues(self):
        gamma_fitter = DataFit("") 
        parameters, covariance = gamma_fitter.fit(self.dataset)

        result = model.approx.values(self.dataset.data, self.dataset.parameters, model.impact_range())
        for a, b in zip(result, self.nominal_value):
                self.assertAlmostEqual(a, b)
