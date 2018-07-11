#!/usr/bin/python2

import json

import ROOT
import pandas as pd


class DataFit(object):
    def __init__(self, model, conffile="config/standard.json"):
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
        self.inpar = self.conf["initial_parameters"]
        self.par_fixed = self.conf["par_fixed"]
        self.par_limits = self.conf["par_limits"]
        self.par_names = self.conf["par_names"]

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

    def fitfunction(self, parameters, ftype="crossection"):
        npar = len(parameters)

        quantity = self.model.diff_cs
        if ftype != "crossection":
            quantity = self.model.ratio

        function = ROOT.TF1(ftype, lambda x, p: quantity(x[0], p), 0, 10, npar)

        for i, par in enumerate(parameters):
            function.SetParameter(i, par)

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
        cs_data.Fit(cs_func, "0NrE")

        dataset.parameters = [
            cs_func.GetParameter(i)
            for i, _ in enumerate(in_parameters)
        ]
        dataset.covariance = self.covariance(len(variables))
        self._save_parameters(cs_func)
        self._save_datapoints(cs_func, dataset)
        self._print_chi2(cs_func, dataset)
        return dataset

    def _save_datapoints(self, cs_func, dataset):
        output = pd.DataFrame()
        output["-t"] = [p.t for p in dataset.data]
        output["obs"] = [p.ds for p in dataset.data]
        output["theory"] = [cs_func.Eval(p.t) for p in dataset.data]
        output["total exp. err"] = [p.err for p in dataset.data]

        par = [
            cs_func.GetParameter(i)
            for i in range(cs_func.GetNpar())
        ]
        amplitudes = [self.model.amplitude(p.t, par)
                      for p in dataset.data]

        output["Re A(s, t)"] = [a.real for a in amplitudes]
        output["Im A(s, t)"] = [a.imag for a in amplitudes]
        oname = "dsigma_{}.csv".format(self.model.name)
        output.to_csv(oname, index=False, sep="\t")

    def _save_parameters(self, cs_func):
        data = {}
        data["par. value"] = [
            cs_func.GetParameter(i)
            for i in range(cs_func.GetNpar())
        ]
        data["par. error"] = [
            cs_func.GetParError(i)
            for i in range(cs_func.GetNpar())
        ]
        names = [
            cs_func.GetParName(i)
            for i in range(cs_func.GetNpar())
        ]
        output = pd.DataFrame(data, index=names)
        output = output.reindex(["par. value", "par. error"], axis=1)
        oname = "parameters-{}.tex".format(self.model.name)
        output.to_latex(oname)

    def _print_chi2(self, cs_func, dataset):
        msg = "{}: chi2/ndof = {}, chi2 = {}, ndf = {}, npoints = {}"
        output = msg.format(
            self.model.name,
            cs_func.GetChisquare() / cs_func.GetNDF(),
            cs_func.GetChisquare(),
            cs_func.GetNDF(),
            len(dataset.data)
        )
        print output
