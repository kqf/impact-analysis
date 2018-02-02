#!/usr/bin/python2


class DataSet(object):
    def __init__(self, infile, parameters):
        super(DataSet, self).__init__()
        self.filename = infile
        self.ptype = parameters["PROCESS"]
        self.energy = parameters["ENERGY"]
        self.dsigma = parameters["DSIGMA"]
        self.sigma = parameters["SIGMA"]
        self.drho = parameters["DRHO"]
        self.rho =  parameters["RHO"]
        self._data = None

    def _read_raw(self, ptype):
        toint1 = lambda x: int( 1000 * float(x) )
        toint2 = lambda x: int( float(x) )

        raw_data = []
        with open(self.filename,'r') as f:
            for line in f:
                data = line.lower().split()
                # Silly but simpliest way to handle problem
                if 1000 * self.energy == 52818 and toint1(data[0]) == 53018:
                    data[0] = 52.818

                if toint1(data[0]) == 1000 * self.energy and toint2(data[6]) == ptype and float(data[1]) > 0.1:
                    raw_data.append([float(data[1]), float(data[2]), float(data[5])])
        return raw_data

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



class DataPoint(object):
    observable = {'pp': 310, 'p#bar{p}': 311}
    def __init__(self, t, ds, err, l = 0, u = 0):
        self.t = t
        self.ds = ds
        self.err = err
        
        self.lower = l
        self.upper = u

    def __lt__(self, other):
        return self.t < other.t

    def __repr__(self):
        # return 't = %f ds/dt = %f' % (self.t, self.ds)
        return 't = %0.4g;\tds/dt = %0.4g\n' % (self.t, self.ds)

    def setBinRange(self, l, u):
        self.lower = l
        self.upper = u

