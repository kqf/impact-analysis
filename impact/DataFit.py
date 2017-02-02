#!/usr/bin/python2
import ROOT
from Formulas import diff_cs, GammaApproximation

class DataFit(object):
    def __init__(self, data, name, title, sigma , rho, nparameters = 12):
        super(DataFit, self).__init__()
        self.canvas = ROOT.TCanvas('canvas', 'Impact Analysis', 800, 600)
        self.data = data
        self.name = name
        self.title = title
        self.sigma = sigma
        self.rho = rho  
        self.nparameters = nparameters


    def differential_cs(self):
        """Creates TGraphErrors, with differential cross section data"""
        # TODO: Check if data red properly
        graph = ROOT.TGraphErrors( len(self.data) )
        graph.SetName(self.name)
        graph.SetTitle(self.title)

        [ graph.SetPoint(i, p.t, p.ds) for i, p in enumerate(self.data) ]
        [ graph.SetPointError(i, 0, p.err) for i, p in enumerate(self.data) ]

        graph.GetXaxis().SetTitle('-t, GeV/c')
        graph.GetYaxis().SetTitle('#frac{d#sigma}{dt}, mb/GeV^{2}')
        graph.SetMarkerStyle(20)
        graph.SetMarkerColor(46)
        return graph


    def differential_cs_approx(self):
        function = ROOT.TF1('function', diff_cs, 0, 10, self.nparameters)

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
        # function.SetParameter(10, self.sigma)
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


    def ratio(self, parameters):
        function = ROOT.TF1('ratio', ratio, 0, 10, self.nparameters)
        # ratio.GetXaxis().SetRange(0, 2)
        [function.SetParameter(i, par) for i, par in enumerate(parameters)]
        function.FixParameter(10, self.sigma)
        function.FixParameter(11, self.rho)
        print 'Ratio at zero is : ', ratio([0], parameters)

        function.SetLineColor(38)
        return function


    def gamma_function(self, parameters):
        """Creates \Gamma(b) functor uses data !"""
        g = GammaApproximation(self.data)
        
        # TODO: try to use functor
        gamma = ROOT.TF1('#Gamma(b)', lambda x, p: g.gamma(x, p), 0, 3, self.nparameters)
        [gamma.SetParameter(i, p) for i, p in enumerate(parameters)]

        gamma.SetLineColor(46)
        gamma.GetXaxis().SetTitle('b\t,fm')
        gamma.GetYaxis().SetTitle('#Gamma')
        return gamma 


    def draw(self):
        self.canvas.Divide(2, 1)
        self.canvas.cd(1).SetLogy()

        graph = self.differential_cs()
        graph.Draw('AP')

        cs_function = self.differential_cs_approx()
        graph.Fit(cs_function,'rE')
        cs_function.Draw('same')

        gamma = self.gamma_function([cs_function.GetParameter(i) for i in range(self.nparameters)])
        self.canvas.cd(2).SetLogy()
        gamma.Draw()
        self.canvas.Update()
        # self.canvas.cd(1)
        # self.legend = self.getLegendForDiffCS()
        # self.legend.Draw()
        raw_input()


        # gamma_0 = self.getGamma([1e-5], parameters)
        # # Needs to be called from outside:
        # self.getReal_Gamma = lambda b : self.getGamma(b, parameters)

        # self.gammaAtZero = gamma_0

        # # print '\Gamma(0) = ', self.gammaAtZero

        # with open('fit_parameters.txt', 'a') as file:
        #     file.write(str(self.energy) + ' ' + str(parameters) + '\n')

        # fitter = TVirtualFitter.GetFitter()
        # # cov = fitter.GetCovarianceMatrix()


        # covariance =  [ [0 for i in range(6)] for j in range(6)]
        # for i in range(6):
        #     for j in range(6):
        #         covariance[i][j]  =  fitter.GetCovarianceMatrixElement(i, j)

        # self.covariance = covariance
        # self.canvas.Update()
        # return [ self.getReal_Gamma( [(1e-5) * (i == 0) + i / 10. ]) for i in range(30)]

