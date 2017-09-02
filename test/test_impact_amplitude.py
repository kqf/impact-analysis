from test.configurable import Configurable
from impact.datapoint import DataReader
from impact.datafit import DataFit

class TestImpactAmplitude(Configurable):

    def setUp(self):
        super(TestImpactAmplitude, self).setUp()
        self.nominal_value = self.data['impact_amplitude']

        data = DataReader(self.ENERGY, self.PROCESS).read(self.infile)
        self.gamma_fitter = DataFit(data, self.PROCESS + str(self.ENERGY), self.PROCESS, self.ENERGY, self.SIGMA, self.RHO) 

    def testValues(self):
        result, _ = self.gamma_fitter.fit()

        for a, b in zip(result, self.nominal_value):
                self.assertAlmostEqual(a, b)
