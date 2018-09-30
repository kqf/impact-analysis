#!/usr/bin/python2

import hashlib

import numpy as np
import pandas as pd
import ROOT

import utils as ut


class NotFittedError(IOError):
    pass


class DataSet(object):
    data_fiels = "s", "-t", "obs", "stat", "syst", "total err.", "code"

    def __init__(self, parameters):
        super(DataSet, self).__init__()
        # Instead of copy constructor
        self._hyperparameters = parameters
        self._hsum = parameters["hsum256"]
        self.filename = parameters["infile"]
        self.ptype = parameters["PROCESS"]
        self.energy = parameters["ENERGY"]
        self.dsigma = parameters["DSIGMA"]
        self.sigma = parameters["SIGMA"]
        self.drho = parameters["DRHO"]
        self.rho = parameters["RHO"]
        self.trange = parameters["TRANGE"]
        self.hadron_minimum = parameters["HADRON_MINIMUM"]
        self.fit_sigma_rho = parameters["FIT_SIGMA_RHO"]
        self.trange[0] = self.trange[0] if self.trange[0] > 0 else -float('inf')
        self.trange[1] = self.trange[1] if self.trange[1] > 0 else float('inf')
        self.index = str(parameters["index"])
        self._data = None
        self._hadron_data = None
        self._parameters = None
        self._covariance = None

        self._validate_dataset()

    def copy(self, data, new_sigma):
        new_set = DataSet(self._hyperparameters)
        new_set._data = data
        new_set.sigma = new_sigma
        new_set._parameters = self._parameters
        new_set._covariance = self._covariance
        return new_set

    def report_chi2(self, model):
        report = self.data_table.copy()
        report["name"] = ["ds/dt"] * len(report)
        print report
        report["theor"] = [
            model.diff_cs(t, self.parameters)
            for t in report["-t"]]
        sigma = {"name": ["sigma"], "-t": [0], "obs": [self.sigma],
                 "total err.": [self.dsigma],
                 "theor": [model.sigma(0, self.parameters)]}

        rho = {"name": ["rho"], "-t": [0], "obs": [self.rho],
               "total err.": [self.drho],
               "theor": [model.ratio(0, self.parameters)]}

        report = report.append(pd.DataFrame(sigma), ignore_index=True)
        report = report.append(pd.DataFrame(rho), ignore_index=True)
        chi2 = (report["obs"] - report["theor"]) / report["total err."]
        report["chi2"] = chi2 ** 2
        print report

    def _validate_dataset(self):
        hashsum = hashlib.sha256()
        with open(self.filename, "rb") as f:
                data = f.read()
        hashsum.update(data)

        # Ignore the hash sum if it's equal to -1
        if int(self._hsum, 16) == -1:
            return

        msg = "You are using a wrong file to test your data. "\
              "Your hash sums don't coincide."\
              "\n\nActual:  {}\nNominal: {}".format(
                  hashsum.hexdigest(), self._hsum)

        assert hashsum.hexdigest() == self._hsum, msg

    def differential_cs(self):
        graph = ROOT.TGraphErrors(len(self.data))
        graph.SetTitle(
            "Differential cross section at {0} TeV".format(self.energy))

        for i, p in enumerate(self.data):
            graph.SetPoint(i, p.t, p.ds)
            graph.SetPointError(i, 0, p.err)

        graph.GetXaxis().SetTitle("-t, GeV/c")
        graph.GetYaxis().SetTitle("#frac{d#sigma}{dt}, mb/GeV^{2}")
        graph.SetMarkerStyle(20)
        graph.SetMarkerColor(46)
        return graph

    def draw(self, canvas=None):
        if not canvas:
            canvas = ut.canvas('dataset')
        graph = self.differential_cs()
        graph.Draw()
        canvas.Update()
        raw_input("Press enter to continue")

    def sigma_rho_differential_cs(self):
        graph = ROOT.TGraphErrors(len(self.data) + 2)
        graph.SetTitle(
            "Differential cross section at {0} TeV".format(self.energy))

        # Artifficially add points at t == 0
        graph.SetPoint(0, -1, self.sigma)
        graph.SetPointError(0, 0, self.dsigma)
        graph.SetPoint(1, -100, self.rho)
        graph.SetPointError(1, 0, self.drho)

        offset = 2
        for i, p in enumerate(self.data):
            graph.SetPoint(offset + i, p.t, p.ds)
            graph.SetPointError(offset + i, 0, p.err)

        graph.GetXaxis().SetTitle("-t, GeV/c")
        graph.GetYaxis().SetTitle("#frac{d#sigma}{dt}, mb/GeV^{2}")
        graph.SetMarkerStyle(20)
        graph.SetMarkerColor(46)
        return graph

    def _read_raw(self, ptype):
        data = pd.read_csv(
            self.filename,
            sep="\s+",
            names=self.data_fiels,
            usecols=[i for i, _ in enumerate(self.data_fiels)]
        )
        data = data[data['code'] == ptype]
        data = data[np.isclose(data['s'], self.energy)]
        data['-t'] = data['-t'].astype('float64')
        data = data[data['-t'].between(*self.trange, inclusive=False)]
        self.data_table = data[['-t', 'obs', 'total err.']]
        return self.data_table.values

    def within_trange(self, t):
        lower, upper = self.trange
        return lower < t < upper

    def read(self):
        ptype = DataPoint.observable[self.ptype]

        d = [DataPoint(*raw) for raw in self._read_raw(ptype)]
        d.sort()

        # Recalculate bin widths properly
        d[0].setBinRange( d[0].t - (d[0].t - d[1].t)/2., ( d[1].t + d[0].t )/2. )
        [d[i].setBinRange(( d[i - 1].t + d[i    ].t )/2., ( d[i    ].t + d[i + 1].t )/2. ) for i in range(1,len(d) - 1)]
        d[-1].setBinRange(( d[-2].t +  d[-1].t)/2 , d[-1].t + ( d[-1].t - d[-2].t )/2. )
        return d

    @property
    def data(self):
        if not self._data:
            self._data = self.read()
        return self._data

    @property
    def hadron_data(self):
        if not self._hadron_data:
            self._hadron_data = [p for p in self.data
                                 if p.t > self.hadron_minimum]
        return self._hadron_data

    @property
    def parameters(self):
        if not self._parameters:
            raise NotFittedError(
                "Trying to use **parameters** without fitting")
        return self._parameters

    @parameters.setter
    def parameters(self, parameters):
        self._parameters = parameters

    @property
    def covariance(self):
        if self._covariance is None:
            raise NotFittedError(
                "Trying to use **covariance** without fitting")
        return self._covariance

    @covariance.setter
    def covariance(self, covariance):
        self._covariance = np.asarray(covariance)
        if self.fit_sigma_rho:
            return
        # When rho and sigma are fixed they are not correlated
        # with other parameters
        self._covariance[:, -1] = 0
        self._covariance[:, -2] = 0
        self._covariance[-1, :] = 0
        self._covariance[-2, :] = 0
        self._covariance[-1][-1] = self.drho ** 2
        self._covariance[-2][-2] = self.dsigma ** 2


class DataPoint(object):
    observable = {"pp": 310, "p#bar{p}": 311}

    def __init__(self, t, ds, err, lower=0, upper=0):
        self.t = t
        self.ds = ds
        self.err = err
        self.lower = lower
        self.upper = upper

    def __lt__(self, other):
        return self.t < other.t

    def __repr__(self):
        # return "t = %f ds/dt = %f" % (self.t, self.ds)
        return "t = %0.4g;\tds/dt = %0.4g\n" % (self.t, self.ds)

    def setBinRange(self, l, u):
        self.lower = l
        self.upper = u
