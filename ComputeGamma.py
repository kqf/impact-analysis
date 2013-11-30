#!/usr/bin/python2.7

from ROOT import *
from DataPoint import *
# from readCSData import *
from Formulas import diff_cs, GammaApproximation
from random import gauss

class ComputeGamma(object):
    __canvas = TCanvas('canvas', 'Impact Analysis', 800, 600)
    __observable = {'pp': 310, 'p#bar{p}': 311}
    
    def __init__(self, ptype, energy, sigma, rho, is_mc = False):
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
        if is_mc:
            self.__generateMC()

    def __readDataCS(self):
        """Reading data from alldata_v1_4.dat file"""
        toint = lambda x: int( float(x) )
        raw_data = []
        with open('alldata_v1_4.dat','r') as file:
            for line in file:
                data = line.lower().split()
                if toint(data[0]) == self.__energy and toint(data[6]) == self.__procesType:
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
        for i, point in enumerate(self.__dataPoints):
            mu = point.ds
            sigma = point.err
            self.__dataPoints[i].ds = gauss(mu, sigma)

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
        # diff_cs = lambda x,p: x[0]
        function = TF1('function', diff_cs, 0, 10, 9)

        function.SetParameter(0, (30.14)**0.5)
        function.SetParameter(1, (11.29)**0.5 )
        function.SetParameter(2, (0.0025)**0.5)
        function.SetParameter(3, (14.3))
        function.SetParameter(4, (7.67))
        function.SetParameter(5, (1.858))
        function.SetParLimits(5, 0, 20)
        function.FixParameter(6, self.sigma) #  sigma0
        function.FixParameter(7, 38.94)      #  sigma0
        function.FixParameter(8, self.rho)   #  rho

        function.SetLineColor(38)

        return function
    def __createGammaFunction(self, parameters):
        # TODO: Create function
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
        # TODO: Initialize functor with  dataPoints
        parameters = [ self.function.GetParameter(i) for i in range(9) ]

        print parameters
        self.__gamma = self.__createGammaFunction(parameters)
        self.__canvas.cd(2)
        gPad.SetLogy()
        self.__gamma.Draw()


        self.__canvas.Update()

    def performComputationsMC(self, nuber_of_points, prefix):
        """Writes points from b = 0 to b = 3 to file"""
        self.graph = self.__createDiffCsGraph()

        self.function = self.__createDiffCsFunct()
        # self.graph.Fit(self.function, 'r')

        # TODO: correct hardcoded parameters
        # parameters = [7.84799185336423, 1.5012093138338845, 1.1213117749438828,
                      # 8.276553417529822, 4.488295630955916, 4.599604734749989,
                      # 97.98762962023775, 39.22381587849738, 0.15588368317954668 ]
        parameters = [ 7.428437251797944, 1.8425102008870191, 1.1699384365535779,
                       8.307894270633678, 4.630583811454606, 4.657402188541418,
                       98.58, 38.94, 0.141]
        # parameters = [ self.function.GetParameter(i) for i in range(9) ]

        gammaComputorMC = GammaApproximation(self.__dataPoints)
        getGamma = lambda x: gammaComputorMC.gamma([x], parameters)

        # TODO: Rewrite using list comprehention
        with open(self.name + str(prefix) + '.dat', 'w') as file:
            i = 0
            while i < 3:
                if i == 0:
                    g = getGamma(1e-5)
                else:
                    g = getGamma(i)
                i += 3./nuber_of_points
                file.write(str(i) + '\t' + str(g) + '\n')
    def getCanvas(self):
        return self.__canvas


def main():
    ENERGY = 7000
    RHO    = 0.141
    SIGMA  = 98.58
    PROCESS = 'pp'

    c = ComputeGamma(PROCESS, ENERGY, SIGMA, RHO, False)
    c.performComputations()
    # c.performComputationsMC(100, 'a')
    raw_input('press any key ...')
    

if __name__ == "__main__":
    main()
else:
    print 'Module loaded'
