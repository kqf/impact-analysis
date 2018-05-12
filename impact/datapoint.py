#!/usr/bin/python2

import hashlib
import ROOT


class NotFittedError(IOError):
    pass


class DataSet(object):
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
        self.index = str(parameters["index"])
        self._data = None
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

    def _validate_dataset(self):
        hashsum = hashlib.sha256()
        with open(self.filename, 'rb') as f:
                data = f.read()
        hashsum.update(data)

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

        graph.GetXaxis().SetTitle('-t, GeV/c')
        graph.GetYaxis().SetTitle('#frac{d#sigma}{dt}, mb/GeV^{2}')
        graph.SetMarkerStyle(20)
        graph.SetMarkerColor(46)
        return graph

    def _read_raw(self, ptype):
        def toint(x):
            return int(1000 * float(x))

        raw_data = []
        with open(self.filename, 'r') as f:
            for line in f:
                data = line.lower().split()
                # Silly but simpliest way to handle problem
                if 1000 * self.energy == 52818 and toint(data[0]) == 53018:
                    data[0] = 52.818

                energy_ok = (toint(data[0]) == 1000 * self.energy)
                if not energy_ok:
                    continue

                type_ok = (int(float(data[6])) == ptype)
                if not type_ok:
                    continue

                t_ok = self.within_trange(float(data[1]))
                if not t_ok:
                    continue

                raw_data.append(
                    [float(data[1]), float(data[2]), float(data[5])]
                )

        return raw_data

    def within_trange(self, t):
        lower, upper = self.trange
        if lower < 0 and upper > 0:
            return t < upper

        if lower > 0 and upper < 0:
            return lower < t
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
        if not self._covariance:
            raise NotFittedError(
                "Trying to use **covariance** without fitting")
        return self._covariance

    @covariance.setter
    def covariance(self, covariance):
        self._covariance = covariance


class DataPoint(object):
    observable = {'pp': 310, 'p#bar{p}': 311}

    def __init__(self, t, ds, err, lower=0, upper=0):
        self.t = t
        self.ds = ds
        self.err = err
        self.lower = lower
        self.upper = upper

    def __lt__(self, other):
        return self.t < other.t

    def __repr__(self):
        # return 't = %f ds/dt = %f' % (self.t, self.ds)
        return 't = %0.4g;\tds/dt = %0.4g\n' % (self.t, self.ds)

    def setBinRange(self, l, u):
        self.lower = l
        self.upper = u
