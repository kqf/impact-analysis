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
        self.nparameters = 12
        self.parametersFile = ''
        self.gamma_fitter = DataFit(self.dataPoints, self.name, self.title, self.energy, self.sigma, self.rho)
        

    def __generateMC(self):
        """Generates MC data"""
        return [DataPoint(p.t, rnd.gauss(p.ds, p.err), p.err, p.lower, p.upper) for p in self.dataPoints]

    def performComputations(self):
        return self.gamma_fitter.fit()


    def performComputationsMC(self, nuber_of_points, prefix, dsigma):
        """Writes points from b = 0 to b = 3 to file"""
        ROOT.gROOT.SetBatch(True)
        file_name = 'parameters_' + self.title + str(self.energy) + '.dat'

        # needs to be deleted from the outside of the class
        self.parametersFile = file_name
        try:
            with open(file_name, 'r') as file:
                data = file.readline().split()
                parameters = [float(i) for i in data]
        except IOError:
            parameters = self.gamma_fitter.get_save_parameters()

        mcPoints = self.__generateMC()

        gammaComputorMC = GammaApproximation(mcPoints)
        new_sigma = rnd.gauss(self.sigma, dsigma)
        getGamma = lambda x: gammaComputorMC.gamma([x], parameters, new_sigma)

        b_max = 3.0
        b = 0
        result = []
        while b <= b_max + 0.000001:
            if b == 0:
                g = getGamma(1e-5)
            else:
                g = getGamma(b)

            b += b_max/nuber_of_points
            # print b, b <= b_max

            result.append(g)
        ROOT.gROOT.SetBatch(False)
        return result

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
