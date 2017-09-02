#!/usr/bin/python2

class DataPoint(object):
    def __init__(self, t_, ds_, err_, l = 0, u = 0):
        self.t = t_
        self.ds = ds_
        self.err = err_
        
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


class DataReader(object):
    observable = {'pp': 310, 'p#bar{p}': 311}
    def __init__(self, energy, ptype):
        super(DataReader, self).__init__()
        self.energy = energy
        self.ptype = self.observable[ptype]
        
    def read_raw(self, filename):
        toint1 = lambda x: int( 1000 * float(x) )
        toint2 = lambda x: int( float(x) )

        raw_data = []
        with open(filename,'r') as f:
            for line in f:
                data = line.lower().split()
                # Silly but simpliest way to handle problem
                if 1000 * self.energy == 52818 and toint1(data[0]) == 53018:
                    data[0] = 52.818

                if toint1(data[0]) == 1000 * self.energy and toint2(data[6]) == self.ptype and float(data[1]) > 0.1:
                    raw_data.append([float(data[1]), float(data[2]), float(data[5])])

        return raw_data

    def read(self, filename = 'alldata_v1_4.dat'):
        d = [DataPoint(i[0], i[1], i[2]) for i in self.read_raw(filename)]
        d.sort()

        # Recalculate bin widths properly
        d[0].setBinRange( d[0].t - (d[0].t - d[1].t)/2., ( d[1].t + d[0].t )/2. )
        [d[i].setBinRange(( d[i - 1].t + d[i    ].t )/2., ( d[i    ].t + d[i + 1].t )/2. ) for i in range(1,len(d) - 1)]
        d[-1].setBinRange(( d[-2].t +  d[-1].t)/2 , d[-1].t + ( d[-1].t - d[-2].t )/2. )

        return d