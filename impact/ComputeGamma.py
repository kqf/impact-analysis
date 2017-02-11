#!/usr/bin/python2.7

from ROOT import *
from DataPoint import DataPoint, DataReader
from DataFit import DataFit
# from readCSData import *
from Formulas import diff_cs, GammaApproximation, ratio, getRealGamma, getRealGammaError, getRealError, getImagGamma, amplitude
import random as rnd

class ComputeGamma(object):
    __canvas = TCanvas('canvas', 'Impact Analysis', 800, 600)
    observable = {'pp': 310, 'p#bar{p}': 311}

    def __init__(self, ptype, energy, sigma, rho):
        self.sigma = sigma
        self.rho = rho
        self.energy = energy
        self.name = ptype + str(energy)
        self.title = ptype
        # Reading data from file
        # TODO: Add try-catch block for check proper ptype
        self.dataPoints = DataReader(self.energy, self.observable[ptype]).read()
        self.parametersFile = ''
        self.gamma_fitter = DataFit(self.dataPoints, self.name, self.title, self.energy, self.sigma, self.rho)
        

    def generate_mc(self):
        """Generates MC data"""
        return [DataPoint(p.t, rnd.gauss(p.ds, p.err), p.err, p.lower, p.upper) for p in self.dataPoints]

    def performComputations(self):
        return self.gamma_fitter.fit()


    def read_parameters(self):
        try:
            with open(self.gamma_fitter.par_file_name, 'r') as file:
                data = file.readline().split()
                return [float(i) for i in data]
        except IOError:
            return self.gamma_fitter.get_save_parameters()

    def performComputationsMC(self, nuber_of_points, prefix, dsigma):
        """Writes points from b = 0 to b = 3 to file"""
        ## Read parameters for real data approximation
        parameters, mcPoints = self.read_parameters(), self.generate_mc()

        ## Get gamma approximator using data points, generate random cross section value
        mc_gamma, new_sigma = GammaApproximation(mcPoints), rnd.gauss(self.sigma, dsigma)
        gamma = lambda x: mc_gamma.gamma([x], parameters, new_sigma)

        ## Calculate values of Gamma function for the mc
        return [gamma((1e-5) * (i == 0) + i * 3.0 / nuber_of_points) for i in range(101)] 

    def t_max(self):
        return self.dataPoints[-1].t

def main():
    ENERGY = 7000
    RHO    = 0.14
    SIGMA  = 98.3

    DSIGMA = 2.23
    DRHO =   0.007

    PROCESS = 'pp'

    c = ComputeGamma(PROCESS, ENERGY, SIGMA, RHO) 
    print c.performComputations()



    # print 'im Gamma', getRealGamma([1e-5], c.parameters)
    # print 'delta im Gamma', getRealGammaError( [1e-5], c.parameters, c.covariance, DSIGMA, DRHO)
    # print 're gamma', c.getGamma([1e-5], c.parameters)
    # print 're gamma', getImagGamma([1e-5], c.parameters)


    raw_input('press any key ...')

if __name__ == "__main__":
    main()
else:
    print 'Module loaded'
