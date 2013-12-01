#!/usr/bin/python2.7

from ROOT import *
from DataPoint import *
# from readCSData import *
from Formulas import diff_cs, GammaApproximation
from random import gauss

class ComputeGamma(object):
    __canvas = TCanvas('canvas', 'Impact Analysis', 800, 600)
    __observable = {'pp': 310, 'p#bar{p}': 311}
    
    def __init__(self, ptype, energy, sigma, rho):
        self.sigma = sigma
        self.rho = rho

        # TODO: Add try-catch block for check proper ptype
        self.__procesType = self.__observable[ptype]
        self.__energy = energy
        self.__dataPoints = []

        self.name = ptype + str(energy)
        self.title = ptype
        # Reading data from file
        self.__readDataCS()

    def __readDataCS(self):
        """Reading data from alldata_v1_4.dat file"""
        toint1 = lambda x: int( 1000 * float(x) )
        toint2 = lambda x: int( float(x) )

        raw_data = []
        with open('alldata_v1_4.dat','r') as file:
            for line in file:
                data = line.lower().split()
                if toint1(data[0]) == 1000*self.__energy and toint2(data[6]) == self.__procesType:
                    raw_data.append([float(data[1]), float(data[2]), float(data[5])])

        self.__dataPoints = [ DataPoint(i[0], i[1], i[2]) for i in raw_data ]
        self.__dataPoints.sort()

        self.__dataPoints[0].setBinRange( self.__dataPoints[0].t - (self.__dataPoints[0].t - self.__dataPoints[1].t)/2. 
                                       ,( self.__dataPoints[1].t + self.__dataPoints[0].t )/2. )

        [ self.__dataPoints[i].setBinRange(
              ( self.__dataPoints[i - 1].t + self.__dataPoints[i    ].t )/2.
            , ( self.__dataPoints[i    ].t + self.__dataPoints[i + 1].t )/2.
            ) for i in range(1,len(self.__dataPoints) - 1) ]

        self.__dataPoints[-1].setBinRange(
                ( self.__dataPoints[-2].t +  self.__dataPoints[-1].t)/2
                , self.__dataPoints[-1].t + ( self.__dataPoints[-1].t - self.__dataPoints[-2].t )/2. )

    def __generateMC(self):
        """Generates MC data"""
        mc_list = list(self.__dataPoints)
        for i, point in enumerate(self.__dataPoints):
            mu = point.ds
            sigma = point.err
            mc_list[i].ds = gauss(mu, sigma)

        return mc_list

    def __createDiffCsGraph(self):
        """Creates TGraphErrors, with differential cross section data"""
        # TODO: Check if data red properly
        graph = TGraphErrors( len(self.__dataPoints) )
        graph.SetName(self.name)
        graph.SetTitle(self.title)

        [ graph.SetPoint(i, p.t, p.ds) for i, p in enumerate(self.__dataPoints) ]
        [ graph.SetPointError(i, 0, p.err) for i, p in enumerate(self.__dataPoints) ]

        graph.GetXaxis().SetTitle('-t, GeV/c')
        graph.GetYaxis().SetTitle('#frac{d#sigma}{dt}, mb/GeV^{2}')
        graph.SetMarkerStyle(20)
        graph.SetMarkerColor(46)

        return graph

    def __createDiffCsFunct(self):
        """Creates TF1 function for fitting data"""

        function = TF1('function', diff_cs, 0, 10, 9)

        function.SetParameter(0, (30.14)**0.5)
        function.SetParameter(1, (11.29)**0.5 )
        function.SetParameter(2, (0.0025)**0.5)
        function.SetParameter(3, (14.3))
        function.SetParameter(4, (7.67))
        function.SetParameter(5, (1.858))
        function.SetParLimits(5, 0, 20)
        function.FixParameter(6, self.sigma) 
        function.FixParameter(7, 38.94)      #  sigma0
        function.FixParameter(8, self.rho)
        # function.SetParameter(6, self.sigma)
        # function.SetParameter(8, self.rho)

        function.SetLineColor(38)
        return function

    def __createGammaFunction(self, parameters):
        """Creates \Gamma(b) functor uses __dataPoints !"""
        self.gammaComputor = GammaApproximation(self.__dataPoints)
        self.getGamma = lambda x, p: self.gammaComputor.gamma(x, p)

        gamma = TF1('#Gamma(b)', self.getGamma,0, 3, 9)
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

        self.function = self.__createDiffCsFunct()
        self.graph.Fit(self.function, 'r')
        self.function.Draw('same')
        parameters = [ self.function.GetParameter(i) for i in range(9) ]

        self.chi2 = self.function.GetChisquare()/ (self.function.GetNDF() if self.function.GetNDF() != 0 else 1)

        print parameters
        self.__gamma = self.__createGammaFunction(parameters)
        self.__canvas.cd(2)
        gPad.SetLogy()
        self.__gamma.Draw()

        self.__canvas.cd(1)
        self.legend = self.getLegendForDiffCS()
        self.legend.Draw()

        gamma_0 = self.getGamma([1e-5], parameters)

        with open('gamma_at_zero.txt', 'a') as file:
            file.write( '%f\t%f' % (self.__energy, gamma_0) )
        self.__canvas.Update()

    def performComputationsMC(self, nuber_of_points, prefix):
        """Writes points from b = 0 to b = 3 to file"""

        file_name = 'parameters_' + self.title + str(self.__energy) + '.dat'
        try:
            with open(file_name, 'r') as file:
                data = file.readline().split()
                parameters = [float(i) for i in data]
        except IOError:
            self.graph    = self.__createDiffCsGraph()
            self.function = self.__createDiffCsFunct()

            self.graph.Fit(self.function, 'r')
            parameters = [ self.function.GetParameter(i) for i in range(9) ]

            with open(file_name, 'w') as file:
                [file.write(str(i) + ' ') for i in parameters]

        mcPoints = self.__generateMC()

        gammaComputorMC = GammaApproximation(self.__dataPoints)
        getGamma = lambda x: gammaComputorMC.gamma([x], parameters)

        b_max = 3.
        b = 0
        result = []
        while b < b_max:
            if b == 0:
                g = getGamma(1e-5)
            else:
                g = getGamma(b)

            b += b_max/nuber_of_points
            result.append(g)
        print '%d percent is done' % (prefix)
        return result

    def getCanvas(self):
        return self.__canvas

    def getLegendForDiffCS(self):
        """Creates legend for ds/dt graph"""
        legend = TLegend(0.9, 0.7, 0.3, 0.8)
        legend.AddEntry(self.graph ,
                '#sqrt{s} = '         + str(self.__energy) +
                'GeV #sigma_{tot} = ' + str(self.sigma   )  +
                'mb #rho = '          + str(self.rho)     )
        legend.AddEntry(self.function , '#chi^{2}/ndf = ' + str(self.chi2))
        return legend


def main():
    ENERGY = 7000
    RHO    = 0.141
    SIGMA  = 98.58
    PROCESS = 'pp'

    c = ComputeGamma(PROCESS, ENERGY, SIGMA, RHO)
    c.performComputations()
    # c.performComputationsMC(100, 'a')
    raw_input('press any key ...')

if __name__ == "__main__":
    main()
else:
    print 'Module loaded'
