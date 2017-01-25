#!/usr/bin/python2.7

from ROOT import *
from DataPoint import DataPoint, DataReader
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
        self.nParemeters = 12
        self.gammaAtZero = 0
        self.parametersFile = ''
        

    def __generateMC(self):
        """Generates MC data"""
        return [DataPoint(p.t, rnd.gauss(p.ds, p.err), p.err, p.lower, p.upper) for p in self.dataPoints]


    def __createDiffCsGraph(self):
        """Creates TGraphErrors, with differential cross section data"""
        # TODO: Check if data red properly
        graph = TGraphErrors( len(self.dataPoints) )
        graph.SetName(self.name)
        graph.SetTitle(self.title)

        [ graph.SetPoint(i, p.t, p.ds) for i, p in enumerate(self.dataPoints) ]
        [ graph.SetPointError(i, 0, p.err) for i, p in enumerate(self.dataPoints) ]

        graph.GetXaxis().SetTitle('-t, GeV/c')
        graph.GetYaxis().SetTitle('#frac{d#sigma}{dt}, mb/GeV^{2}')
        graph.SetMarkerStyle(20)
        graph.SetMarkerColor(46)

        return graph

    def __createDiffCsFunct(self):
        """Creates TF1 function for fitting data"""

        function = TF1('function', diff_cs, 0, 10, self.nParemeters)

        # parameters = [0.2867574, 0, 0.3005147e-01, 0.1318668e+01
                    # , 0.2405179e+01, 0, 0.1821596e+01, 0 ,0.7800825e-01, 0]
                    # # 0, 0.4509552e+01, 0]
        # parameters = [0.1198810, 0, 0.7276397, 0.7088795,
                      # 0.1382631e+01, 0, 0.4637189e+01, 0, 0.5886629e+00, 0]
        parameters = [0.11986832441123918, 0.0, 1.1660221228353649, 0.44233049876624964,
                       0.8627662804403674, 0.0, 4.63711534711051, 0.0, 0.588952821602961, 0.0 ]

        self.parameters = parameters


        [function.SetParameter(i, par) for i, par in enumerate(parameters)]


        
        function.FixParameter(10, self.sigma)
        function.FixParameter(11, self.rho)

#         function.SetParameter(10, self.sigma)
        # function.SetParameter(11, self.rho)

        # function.SetParLimits(10, self.sigma - 2.2, self.sigma + 2.2)
        # function.SetParLimits(11, self.rho - 0.08, self.rho + 0.01)

        # function.SetParLimits(1, -100, 0)
        function.FixParameter(1, 0)
        function.FixParameter(5, 0)
        function.FixParameter(7, 0)
        function.FixParameter(9, 0)


        function.SetParLimits(0, 0, 6)
        function.SetParLimits(2, -50, 250)
        function.SetParLimits(3, 0, 50)
        function.SetParLimits(4, 0.1, 50)
        function.SetParLimits(6, 0, 100)
        function.SetParLimits(8, 0.078, 100)
        # function.FixParameter(8, 0.078)

        function.SetLineColor(38)
        return function

    def __createRatioFunction(self, parameters):
        """Creates TF1 function for fitting data"""

        function = TF1('ratio', ratio, 0, 10, self.nParemeters)
        # ratio.GetXaxis().SetRange(0, 2)
        [function.SetParameter(i, par) for i, par in enumerate(parameters)]

        function.FixParameter(10, self.sigma)
        function.FixParameter(11, self.rho)
        print 'Ratio at zero is : ', ratio([0], parameters)

        function.SetLineColor(38)
        return function

    def __createGammaFunction(self, parameters):
        """Creates \Gamma(b) functor uses dataPoints !"""
        self.gammaComputor = GammaApproximation(self.dataPoints)
        self.getGamma = lambda x, p: self.gammaComputor.gamma(x, p)

        gamma = TF1('#Gamma(b)', self.getGamma,0, 3, self.nParemeters)
        [gamma.SetParameter(i, p) for i, p in enumerate(parameters)]

        gamma.SetLineColor(46)
        gamma.GetXaxis().SetTitle('b\t,fm')
        gamma.GetYaxis().SetTitle('#Gamma')
        return gamma

    def performComputations(self):
        self.__canvas.Divide(2, 1)
        self.__canvas.cd(1)
        gPad.SetLogy()

        self.graph = self.__createDiffCsGraph()
        self.graph.Draw('AP')
        print self.graph.GetN()

        self.function = self.__createDiffCsFunct()
        self.graph.Fit(self.function,'rE')

        self.function.Draw('same')
        parameters = [ self.function.GetParameter(i) for i in range(self.nParemeters) ]
        parametersErrors = [ self.function.GetParError(i) for i in range(self.nParemeters) ]
        # Needed to compute gamma at zero
        self.parameters = parameters
        self.parametersErrors = parametersErrors

        self.chi2 = self.function.GetChisquare()/ (self.function.GetNDF() if self.function.GetNDF() != 0 else 1)

        self.__gamma = self.__createGammaFunction(parameters)
        self.__canvas.cd(2)
        gPad.SetLogy()
        self.__gamma.Draw()

        self.__canvas.cd(1)
        self.legend = self.getLegendForDiffCS()
        self.legend.Draw()


        gamma_0 = self.getGamma([1e-5], parameters)
        # Needs to be called from outside:
        self.getReal_Gamma = lambda b : self.getGamma(b, parameters)

        self.gammaAtZero = gamma_0

        # print '\Gamma(0) = ', self.gammaAtZero

        with open('fit_parameters.txt', 'a') as file:
            file.write(str(self.energy) + ' ' + str(parameters) + '\n')

        fitter = TVirtualFitter.GetFitter()
        # cov = fitter.GetCovarianceMatrix()


        covariance =  [ [0 for i in range(6)] for j in range(6)]
        for i in range(6):
            for j in range(6):
                covariance[i][j]  =  fitter.GetCovarianceMatrixElement(i, j)

        self.covariance = covariance
        self.__canvas.Update()
        return [ self.getReal_Gamma( [(1e-5) * (i == 0) + i / 10. ]) for i in range(30)]

    def performComputationsMC(self, nuber_of_points, prefix, dsigma):
        """Writes points from b = 0 to b = 3 to file"""

        file_name = 'parameters_' + self.title + str(self.energy) + '.dat'
        # needs to be deleted from the outside of the class
        self.parametersFile = file_name
        try:
            with open(file_name, 'r') as file:
                data = file.readline().split()
                parameters = [float(i) for i in data]
        except IOError:
            self.graph    = self.__createDiffCsGraph()
            self.function = self.__createDiffCsFunct()

            self.graph.Fit(self.function, 'rE')
            parameters = [ self.function.GetParameter(i) for i in range(self.nParemeters) ]
            with open(file_name, 'w') as file:
                [ file.write(str(i) + ' ') for i in parameters ]

        mcPoints = self.__generateMC()

        gammaComputorMC = GammaApproximation(self.dataPoints)
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

        print '%d percent is done' % (prefix)
        return result

    def getCanvas(self):
        return self.__canvas

    def getLegendForDiffCS(self):
        """Creates legend for ds/dt graph"""
        legend = TLegend(0.9, 0.7, 0.3, 0.8)
        legend.AddEntry(self.graph ,
                '#sqrt{s} = '         + str(self.energy) +
                'GeV #sigma_{tot} = ' + str(self.sigma   ) +
                'mb #rho = '          + str(self.rho)    )
        legend.AddEntry(self.function , '#chi^{2}/ndf = ' + str(self.chi2))
        return legend
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
