#!/usr/bin/python2.7

from DataPoint import DataPoint, DataReader
from DataFit import DataFit
from Formulas import GammaApproximation 
import random as rnd

class ComputeGamma(object):
    observable = {'pp': 310, 'p#bar{p}': 311}
    def __init__(self, infile, ptype, energy, sigma, rho):
        self.sigma = sigma
        self.dataPoints = DataReader(energy, self.observable[ptype]).read(infile)
        self.gamma_fitter = DataFit(self.dataPoints, ptype + str(energy), ptype, energy, self.sigma, rho)
        self.parameters = None
        

    def generate_mc(self):
        """Generates MC data"""
        return [DataPoint(p.t, rnd.gauss(p.ds, p.err), p.err, p.lower, p.upper) for p in self.dataPoints]


    def performComputations(self):
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
        mc_gamma, new_sigma = GammaApproximation(mcPoints), rnd.gauss(self.sigma, dsigma)
        gamma = lambda x: mc_gamma.gamma([x], parameters, new_sigma)

        ## Calculate values of Gamma function for the mc
        return [gamma(value) for value in self.impact_range(npoints)] 


    def impact_range(self, npoints):
        return (((1e-5) * (i == 0) + i * 3.0 / npoints) for i in range(101))


    def get_gamma(self, x):
        parameters = self.read_parameters() 
        computor = GammaApproximation(self.dataPoints)
        return computor.gamma([x], parameters, self.sigma)


def main():
    ENERGY = 7000
    RHO    = 0.14
    SIGMA  = 98.3

    DSIGMA = 2.23
    DRHO =   0.007

    PROCESS = 'pp'

    c = ComputeGamma(PROCESS, ENERGY, SIGMA, RHO) 
    print c.performComputations()
    raw_input('press any key ...')

if __name__ == "__main__":
    main()