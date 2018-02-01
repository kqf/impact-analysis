#!/usr/bin/python2
import ROOT
import json
import model
import utils as ut

def getGraph(lst):
    graph = ROOT.TGraphErrors()
    graph.SetName('graph')
    graph.SetTitle('graph')
    [graph.SetPoint(i, i * 3./100, p[0]) for i, p in enumerate(lst)]
    [graph.SetPointError(i, 0, p[1]) for i, p in enumerate(lst)]
    return graph

class DataFit(object):
    with open('config/datafit.json') as f:
        conf = json.load(f)

    def __init__(self, name, title, energy, sigma , rho, mode = 's'):
        super(DataFit, self).__init__()
        self.name = name
        self.title = title
        self.energy = energy
        self.sigma = sigma
        self.rho = rho  
        self.mode = mode
        self.par_file_name = 'output/parameters_' + self.title + str(self.energy) + '.dat'

        # Configure fit
        self.inpar      = self.conf['initial_parameters'] + [self.sigma, self.rho]
        self.par_fixed  = self.conf['par_fixed']
        self.par_limits = self.conf['par_limits']
        self.step_nstep = self.conf['step_nstep']
        self.zero       = self.conf['zero']
        self.cov_size   = self.conf['cov_size']
        self.legend     = self.conf['legend']
        self.t_range    = self.conf['t_range']
        self.b_range    = self.conf['b_range']
        imgfile         = self.conf['imgfile']
        self.ofilename = imgfile % name

        # NB: Keep all objects that you want to reuse 
        #     othervise those objects will be delted
        self.cache = []


    def differential_cs(self, data):
        """
            Creates TGraphErrors, with differential cross section data
        """
        graph = ROOT.TGraphErrors(len(data))
        graph.SetName(self.name)
        graph.SetTitle(self.title)

        for i, p in enumerate(data):
            graph.SetPoint(i, p.t, p.ds)
            graph.SetPointError(i, 0, p.err)

        graph.GetXaxis().SetTitle('-t, GeV/c')
        graph.GetYaxis().SetTitle('#frac{d#sigma}{dt}, mb/GeV^{2}')
        graph.SetMarkerStyle(20)
        graph.SetMarkerColor(46)
        return graph


    def differential_cs_approx(self):
        tmin, tmax = self.t_range
        function = ROOT.TF1('function', lambda x, p: model.diff_cs(x[0], p), 
            tmin, tmax, len(self.inpar))

        for i, par in enumerate(self.inpar):
            function.SetParameter(i, par)

        function.FixParameter(10, self.sigma)
        function.FixParameter(11, self.rho)

        map(lambda x: function.FixParameter(*x), self.par_fixed)
        map(lambda x: function.SetParLimits(*x), self.par_limits)
        # function.FixParameter(8, 0.078)

        function.SetLineColor(38)
        return function


    def ratio(self, parameters):
        tmin, tmax = self.t_range
        function = ROOT.TF1('ratio', model.ratio, tmin, tmax, len(self.inpar))
        # ratio.GetXaxis().SetRange(0, 2)

        for i, par in enumerate(parameters):
            function.SetParameter(i, par)

        function.FixParameter(10, self.sigma)
        function.FixParameter(11, self.rho)
        print 'Ratio at zero is : ', ratio([0], parameters)

        function.SetLineColor(38)
        return function


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




    def fit(self, data):
        canvas = ut.canvas("data")
        # self.decorate_pad(canvas.cd(1))

        cs_data, cs_func = self.differential_cs(data), self.differential_cs_approx()
        cs_data.Draw('AP')
        cs_data.Fit(cs_func,'rE')
        cs_func.Draw('same')

        self.legend = self.get_legend(cs_data, cs_func)
        self.legend.Draw()

        # TODO: replace self.inpar and out_parameters
        out_parameters = [cs_func.GetParameter(i) for i in range(len(self.inpar))]
        gamma = model.approx.tf1(data, out_parameters, *self.b_range)
        gamma.Draw()

        # self.decorate_pad(canvas.cd(2))
        canvas.Update()

        if 's' not in self.mode:
            raw_input('pease enter any key ...')
           
        self.cache.append([cs_data, cs_func, gamma])
        return out_parameters, self.get_covariance()

    # TODO: How do we read and store parameters
    def get_covariance(self):
        fitter = ROOT.TVirtualFitter.GetFitter()
        cov = fitter.GetCovarianceMatrix()
        func = lambda i, j: fitter.GetCovarianceMatrixElement(i, j) 

        covariance = [[func(i, j) 
            for i in range(self.cov_size)] for j in range(self.cov_size)]
        return covariance

    def compare_results(self, data, mu, sigma, parameters):
        canvas = ut.canvas("gamma")
        ut.decorate_pad(canvas)

        true_gamma = model.approx.values(data, parameters)
        gamma_with_error = zip(true_gamma, sigma)

        final_result = getGraph(gamma_with_error)
        final_result.SetLineColor(37)
        final_result.SetMarkerColor(37)

        mc_average = zip(mu, [0 for i in mu])
        average_mc = getGraph(mc_average)
        average_mc.SetLineColor(46)
        average_mc.SetMarkerColor(46)

        final_result.Draw()
        average_mc.Draw('same')
        canvas.Update()
        canvas.SaveAs(self.ofilename)

        if 's' not in self.mode:
            raw_input('pease enter any key ...')
            
        return true_gamma