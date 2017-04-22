#!/usr/bin/python2
import ROOT
import json
from impact.model import diff_cs, ratio, GammaApproximation

class DataFit(object):
    with open('config/datafit.json') as f:
        conf = json.load(f)

    def __init__(self, data, name, title, energy, sigma , rho):
        super(DataFit, self).__init__()
        self.canvas = ROOT.TCanvas('c1', 'Impact Analysis', 800, 600)
        self.data = data
        self.name = name
        self.title = title
        self.energy = energy
        self.sigma = sigma
        self.rho = rho  
        self.par_file_name = 'parameters_' + self.title + str(self.energy) + '.dat'

        # Configure fit
        self.parameters = self.conf['initial_parameters'] + [self.sigma, self.rho]
        self.par_fixed  = self.conf['par_fixed']
        self.par_limits = self.conf['par_limits']
        self.step_nstep = self.conf['step_nstep']
        self.zero       = self.conf['zero']
        self.cov_size   = self.conf['cov_size']
        self.legend     = self.conf['legend']
        self.t_range    = self.conf['t_range']
        self.b_range    = self.conf['b_range']
        self.nparameters = len(self.parameters)

        # Don't initialize these two without purpose
        self.graph, self.gamma = None, None


    def differential_cs(self):
        """
            Creates TGraphErrors, with differential cross section data
        """
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
        tmin, tmax = self.t_range
        function = ROOT.TF1('function', lambda x, p: diff_cs(x[0], p), tmin, tmax, self.nparameters)
        [function.SetParameter(i, par) for i, par in enumerate(self.parameters)]
        function.FixParameter(10, self.sigma)
        function.FixParameter(11, self.rho)

        map(lambda x :function.FixParameter(*x), self.par_fixed)
        map(lambda x :function.SetParLimits(*x), self.par_limits)
        # function.FixParameter(8, 0.078)

        function.SetLineColor(38)
        return function


    def ratio(self, parameters):
        tmin, tmax = self.t_range
        function = ROOT.TF1('ratio', ratio, tmin, tmax, self.nparameters)
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
        bmin, bmax = self.b_range
        gamma = ROOT.TF1('#Gamma(b)', lambda x, p: g(x[0], parameters), bmin, bmax, self.nparameters)
        [gamma.SetParameter(i, p) for i, p in enumerate(parameters)]

        gamma.SetLineColor(46)
        gamma.GetXaxis().SetTitle('b\t,fm')
        gamma.GetYaxis().SetTitle('#Gamma')
        return gamma 


    def get_legend(self, graph, function):
        legend = ROOT.TLegend(*self.legend)
        chi2 = function.GetChisquare() / (function.GetNDF() + 1 * (function.GetNDF() == 0))
        legend.AddEntry(graph , '#sqrt{s} = ' + str(self.energy) + ' GeV')
        legend.AddEntry(0, '#sigma_{tot} = ' + str(self.sigma) + 'mb #rho = ' + str(self.rho), '')
        legend.AddEntry(function , '#chi^{2}/ndf = ' + str(chi2))
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.04) 
        return legend


    def decorate_pad(self, pad):
        pad.SetLogy()
        pad.SetTickx()
        pad.SetTicky()
        pad.SetGridy()
        pad.SetGridx()


    def fit(self):
        self.canvas.Divide(2, 1)
        self.decorate_pad(self.canvas.cd(1))

        self.graph = self.differential_cs()
        self.graph.Draw('AP')

        cs_function = self.differential_cs_approx()
        self.graph.Fit(cs_function,'rE')
        self.covariance = self.get_covariance()

        self.parameters = [cs_function.GetParameter(i) for i in range(self.nparameters)]
        cs_function.Draw('same')

        self.legend = self.get_legend(self.graph, cs_function)
        self.legend.Draw()

        self.gamma = self.gamma_function([cs_function.GetParameter(i) for i in range(self.nparameters)])
        self.decorate_pad(self.canvas.cd(2))
        self.gamma.Draw()
        self.canvas.Update()

        return map(self.gamma.Eval, self.impact_range())


    def impact_range(self, npoints = None, step = None):
        if not npoints:
            step, npoints = self.step_nstep
        return (self.zero * (i == 0) + i / step for i in range(npoints))


    def get_save_parameters(self):
        with open(self.par_file_name, 'w') as f:
            [f.write(str(i) + ' ') for i in self.parameters]
        return self.parameters



    def get_covariance(self):
        fitter = ROOT.TVirtualFitter.GetFitter()
        cov = fitter.GetCovarianceMatrix()
        func = lambda i, j: fitter.GetCovarianceMatrixElement(i, j) 

        covariance = [[func(i, j) 
            for i in range(self.cov_size)] for j in range(self.cov_size)]

        return covariance