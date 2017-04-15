#!/usr/bin/python2.7

from DataPoint import DataPoint, DataReader
from DataFit import DataFit
from Formulas import GammaApproximation 
import random as rnd

class ComputeGamma(object):
    observable = {'pp': 310, 'p#bar{p}': 311}
    def __init__(self, infile, ptype, energy, sigma, rho, mode = 's'):
        self.mode = mode
        self.sigma = sigma
        self.dataPoints = DataReader(energy, self.observable[ptype]).read(infile)
        self.gamma_fitter = DataFit(self.dataPoints, ptype + str(energy), ptype, energy, self.sigma, rho)
        self.parameters = None
        

    def generate_mc(self):
        """Generates MC data"""
        return [DataPoint(p.t, rnd.gauss(p.ds, p.err), p.err, p.lower, p.upper) for p in self.dataPoints]


    def compute(self):
        return self.gamma_fitter.fit()


    def read_parameters(self):
        if not self.parameters:
            self.parameters = self.gamma_fitter.get_save_parameters()
        return self.parameters


    def generate_mc_gamma(self, npoints, prefix, dsigma):
        """Writes points from b = 0 to b = 3 to file"""
        ## Read parameters for real data approximation
        parameters, mcPoints = self.read_parameters(), self.generate_mc()

        ## Get gamma approximator using data points, generate random cross section value
        mc_gamma = GammaApproximation(mcPoints, rnd.gauss(self.sigma, dsigma))
        gamma = lambda x: mc_gamma(x, parameters)

        ## Calculate values of Gamma function for the mc
        return map(gamma, self.impact_range(npoints))


    def impact_range(self, npoints):
        # TODO: merge this method with similar one in DataFit
        return (((1e-5) * (i == 0) + i * 3.0 / npoints) for i in range(101))


    def get_gamma(self, x):
        parameters = self.read_parameters() 
        # TODO: Do we need to pass sigma here?
        computor = GammaApproximation(self.dataPoints, self.sigma)
        return computor(x, parameters)


if __name__ == "__main__":
    main()