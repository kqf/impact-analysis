#!/usr/bin/python2.7

from math  import exp, pi, fabs, sqrt, log
from scipy import integrate
from scipy.special import j0, j1
from constants import k_fm, k_norm
from parametrization.symbolic import Symbolic, SymbolicUpdated

import ROOT

# The parametrization used in here
# $$A(s, t) = i (1-i \rho ) \left(a_{4}+\frac{a_{s}}{4 \pi  k_{norm}}\right) \left(a_{1} e^{-0.5 b_{1} (1-i \rho ) t \left(a_{4}+\frac{a_{s}}{4 \pi  k_{norm}}\right)}+(1-a_{1}) e^{-0.5 b_{2} (1-i \rho ) t \left(a_{4}+\frac{a_{s}}{4 \pi k_{norm}}\right)}\right)-i a_{4} e^{-0.5 b_{4} t}-\frac{a_{4} \rho }{\left(\frac{t}{b_{5}}+1\right)^4}$$

def hankel_transform(func):
    def impact_version(b, p, limits = (0, float("inf"))):
        f = lambda q : q * j0(b * q / k_fm) *  func(q * q, p) / sqrt(pi * k_norm)
        result = integrate.quad(f, *limits)[0]  # integral from zero to lower bound
        return result

    return impact_version


class RealGammaEstimator(object):

    def __init__(self, model, outname="real_gamma"):
        self.model = model
        self.outname = outname

    def evaluate(self, dataset, output):
        output[self.outname] = map(lambda x: self.model.real_gamma(x, dataset.parameters), output.index)
        return output[self.outname].values


class ImageGammaEstimator(object):

    def __init__(self, model, outname="image_gamma"):
        self.model = model
        self.outname = outname

    def evaluate(self, dataset, output):
        output[self.outname] = Approx.values(self.model, dataset.data, dataset.parameters, output.index)
        return output[self.outname].values


def impact_range(npoints = 30, step = 10.0, zero = 1e-5):
    return (zero * (i == 0) + i / step for i in range(npoints))


class Approx(object):
    
    def __init__(self, model, data, new_sigma = None):
        self.model = model
        self.dataPoints = data
        self.new_sigma = new_sigma 

        # Define these functions in a very strange way to avoid code replication
        @hankel_transform
        def bspace_low_t(b, p):
            return self.im_amplitude_low_t(b, p)
        self.bspace_low_t = bspace_low_t


    def im_amplitude_low_t(self, t, p):
        """Calculates imaginary part of extrapolated amplitude"""
        parameters = [i for i in p]

        sigma = self.new_sigma if self.new_sigma else parameters[-2] 

        a0 = sigma / (4 * sqrt(pi * k_norm) )
        a1 = sqrt(self.dataPoints[0].ds - self.model.amplitude(self.dataPoints[0].t, p).real ** 2)
        b0 = (1./(self.dataPoints[0].t)) * log(a0 / a1)

        result = a0 * exp(-1. * b0 * t)
        return result

    def integral(self, b, p, i):
        q1, q2 = sqrt(i.lower), sqrt(i.upper)
        return sqrt(fabs(i.ds - self.model.amplitude(i.t, p).real ** 2)) * (q2 * j1(b * q2 / k_fm) -  q1 * j1(b * q1 / k_fm))


    def __call__(self ,b, p):
        """Calculates gamma(b) using data points"""

        # Taking into account low-t extrapolation
        extrapolation1 = self.bspace_low_t(b, p, (0, self.dataPoints[0].lower ** 0.5))


        # High t-contribution
        extrapolation2 = self.model.imag_gamma(b, p, (self.dataPoints[-1].upper ** 0.5, float("inf")))


        # Integrated values from data
        gamma_data = sum(self.integral(b, p, i) for i in self.dataPoints) 

        # All contributions
        result = extrapolation1 + extrapolation2 + (gamma_data * k_fm / sqrt(pi * k_norm) / b)
        return result


    @classmethod
    def _func(klass, model, data, parameters, new_sigma = None):
        self = klass(model, data, new_sigma)
        return lambda b: self(b, parameters)


    @classmethod
    def tf1(klass, model, data, parameters, bmin, bmax):
        """Creates \Gamma(b) functor uses data !"""
        g = klass(model, data)
        gamma = ROOT.TF1('#Gamma(b)', lambda x, p: g(x[0], parameters), bmin, bmax, len(parameters))

        for i, p in enumerate(parameters):
            gamma.SetParameter(i, p)

        gamma.SetLineColor(46)
        gamma.GetXaxis().SetTitle('b\t,fm')
        gamma.GetYaxis().SetTitle('#Gamma')
        return gamma 

    @classmethod
    def values(klass, model, data, parameters, index, new_sigma = None):
        gamma_lambda = klass._func(model, data, parameters, new_sigma)
        gamma = map(gamma_lambda, index)
        return gamma