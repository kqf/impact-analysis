from test.configurable import Configurable
from impact.errors_image import ImageError
from impact.datapoint import DataReader
from impact.datafit import DataFit
import random

class TestMCImpactAmplitude(Configurable):
    """
        This test checks if the procedure of mc impact amplitude generation 
        remains the same. 
    """

    def setUp(self):
        super(TestMCImpactAmplitude, self).setUp()
        self.nominal_value = self.data['mc_impact_amplitude']
        random.seed(1234)


    def testValues(self):
        nmc = 100
        data = DataReader(self.ENERGY, self.PROCESS).read(self.infile)
        gamma_fitter = DataFit(data, self.PROCESS + str(self.ENERGY), self.PROCESS, self.ENERGY, self.SIGMA, self.RHO)
        imag_errors = ImageError(data, nmc, self.SIGMA, self.DSIGMA)

        _, parameters = gamma_fitter.fit()

        result = imag_errors.generate_mc_gamma(nmc, parameters, self.DSIGMA)

        for a, b in zip(result, self.nominal_value):
                self.assertAlmostEqual(a, b)
