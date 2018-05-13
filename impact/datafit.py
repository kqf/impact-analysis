#!/usr/bin/python2

import json

import ROOT


class DataFit(object):
    def __init__(self, model, conffile='config/datafit.json'):
        super(DataFit, self).__init__()
        self.model = model
        # NB: Keep all objects that you want to reuse
        #     othervise those objects will be delted
        # Configure the fitting functions
        self.configure(conffile)

    def configure(self, conffile):
        with open(conffile) as f:
            self.conf = json.load(f)
        # Configure fit
        #
        self.inpar = self.conf['initial_parameters']
        self.par_fixed = self.conf['par_fixed']
        self.par_limits = self.conf['par_limits']
        self.par_names = self.conf['par_names']

    def covariance(self, size):
        fitter = ROOT.TVirtualFitter.GetFitter()
        # cov = fitter.GetCovarianceMatrix()

        def func(i, j):
            return fitter.GetCovarianceMatrixElement(i, j)

        covariance = [
            [func(i, j) for i in range(size)]
            for j in range(size)
        ]
        return covariance

    def fitfunction(self, parameters, ftype='crossection'):
        npar = len(parameters)

        quantity = self.model.diff_cs if ftype == 'crossection' else self.model.ratio
        function = ROOT.TF1(ftype, lambda x, p: quantity(x[0], p), 0, 10, npar)

        for i, par in enumerate(parameters):
            function.FixParameter(i, par)

        function.FixParameter(npar - 2, parameters[-2])
        function.FixParameter(npar - 1, parameters[-1])
        function.SetParNames(*self.par_names)

        map(lambda x: function.FixParameter(*x), self.par_fixed)
        map(lambda x: function.SetParLimits(*x), self.par_limits)

        function.SetLineColor(38)
        return function

    def fit(self, dataset):
        variables = self.inpar[dataset.index] 
        in_parameters = variables + [dataset.sigma, dataset.rho]

        cs_data = dataset.differential_cs()
        cs_func = self.fitfunction(in_parameters)
        cs_data.Fit(cs_func, '0NrE')

        dataset.parameters = [cs_func.GetParameter(
            i) for i, _ in enumerate(in_parameters)]
        dataset.covariance = self.covariance(len(variables))
        return dataset
