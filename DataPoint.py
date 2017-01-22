#!/usr/bin/python2

class DataPoint(object):
    def __init__(self, t_, ds_, err_):
        self.t = t_
        self.ds = ds_
        self.err = err_

        self.lower = 0
        self.upper = 0

    def __lt__(self, other):
        return self.t < other.t

    def __repr__(self):
        # return 't = %f ds/dt = %f' % (self.t, self.ds)
        return 't = %d ds/dt = %d \n' % (self.t, self.ds)

    def setBinRange(self, l, u):
        self.lower = l
        self.upper = u
