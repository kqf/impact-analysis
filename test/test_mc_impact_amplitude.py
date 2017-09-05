from test.configurable import Configurable
import impact.errors.imag as err_imag
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
        self.longMessage = True


    def testValues(self):
        nmc = 100
        data = DataReader(self.ENERGY, self.PROCESS).read(self.infile)
        gamma_fitter = DataFit(data, self.PROCESS + str(self.ENERGY), self.PROCESS, self.ENERGY, self.SIGMA, self.RHO)
        imag_errors = err_imag.Error(data, self.SIGMA, self.DSIGMA)

        _, parameters = gamma_fitter.fit()

        result = imag_errors.generate_mc_gamma(parameters)

        msg = "Computed values are:\n{0}".format(result)
        for a, b in zip(result, self.nominal_value):
                self.assertAlmostEqual(a, b, msg = msg)
