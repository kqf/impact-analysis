#!/usr/bin/python2
import ROOT
from Formulas import diff_cs, GammaApproximation

class DataFit(object):
    def __init__(self, data, name, title, energy, sigma , rho):
        super(DataFit, self).__init__()
        self.canvas = ROOT.TCanvas('c1', 'Impact Analysis', 800, 600)
        self.data = data
        self.name = name
        self.title = title
        self.energy = energy
        self.sigma = sigma
        self.rho = rho  
        self.parameters = [0.11986832441123918, 0.0, 1.1660221228353649, 0.44233049876624964, 
                           0.8627662804403674, 0.0, 4.63711534711051, 0.0, 0.588952821602961, 0.0, self.sigma, self.rho]
        self.nparameters = len(self.parameters)
        self.par_file_name = 'parameters_' + self.title + str(self.energy) + '.dat'


    def differential_cs(self):
        """Creates TGraphErrors, with differential cross section data"""
        # TODO: Check if data red properly
        graph = ROOT.TGraphErrors( len(self.data) )
        graph.SetName(self.name)
        graph.SetTitle(self.title)

        [graph.SetPoint(i, p.t, p.ds) for i, p in enumerate(self.data)]
        [graph.SetPointError(i, 0, p.err) for i, p in enumerate(self.data)]

        graph.GetXaxis().SetTitle('-t, GeV/c')
        graph.GetYaxis().SetTitle('#frac{d#sigma}{dt}, mb/GeV^{2}')
        graph.SetMarkerStyle(20)
        graph.SetMarkerColor(46)
        return graph


    def differential_cs_approx(self):
        function = ROOT.TF1('function', diff_cs, 0, 10, self.nparameters)
        [function.SetParameter(i, par) for i, par in enumerate(self.parameters)]
        function.FixParameter(10, self.sigma)
        function.FixParameter(11, self.rho)

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


    def get_legend(self, graph, function):
        legend = ROOT.TLegend(0.9, 0.7, 0.3, 0.8)
        chi2 = function.GetChisquare() / (function.GetNDF() + 1 * (function.GetNDF() == 0))
        legend.AddEntry(graph ,
                '#sqrt{s} = '         + str(self.energy) +
                'GeV #sigma_{tot} = ' + str(self.sigma   ) +
                'mb #rho = '          + str(self.rho)    )
        legend.AddEntry(function , '#chi^{2}/ndf = ' + str(chi2))
        return legend


    def fit(self):
        self.canvas.Divide(2, 1)
        self.canvas.cd(1).SetLogy()

        graph = self.differential_cs()
        graph.Draw('AP')

        cs_function = self.differential_cs_approx()
        graph.Fit(cs_function,'rE')

        self.parameters = [cs_function.GetParameter(i) for i in range(self.nparameters)]
        cs_function.Draw('same')

        self.legend = self.get_legend(graph, cs_function)
        self.legend.Draw()

        gamma = self.gamma_function([cs_function.GetParameter(i) for i in range(self.nparameters)])
        self.canvas.cd(2).SetLogy()
        gamma.Draw()
        self.canvas.Update()

        # real_gamma = lambda b : self.getGamma(b, parameters) 
        return [gamma.Eval((1e-5) * (i == 0) + i / 10.) for i in range(30)]


    def get_save_parameters(self):
        with open(self.par_file_name, 'w') as f:
            [f.write(str(i) + ' ') for i in self.parameters]
        return self.parameters



    def covariance(self):
        fitter = TVirtualFitter.GetFitter()
        cov = fitter.GetCovarianceMatrix()

        covariance =  [ [0 for i in range(6)] for j in range(6)]
        for i in range(6):
            for j in range(6):
                covariance[i][j]  =  fitter.GetCovarianceMatrixElement(i, j)
        return covariance