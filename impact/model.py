#!/usr/bin/python2.7

from math  import exp, pi, fabs, sqrt, log
from scipy import integrate
from scipy.special import j0, j1
from constants import k_fm, k_norm
from parametrization.symbolic import Symbolic, SymbolicUpdated

import ROOT

from impact.utils import hankel_transform, impact_range


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
        output[self.outname] = Approx.values(self.model, dataset, output.index)
        return output[self.outname].values


class Approx(object):
    
    def __init__(self, model, dataset):
        self.model = model
        self.dataset = dataset

        # Define these functions in a very strange way to avoid code replication
        @hankel_transform
        def bspace_low_t(b, p):
            return self.im_amplitude_low_t(b, p)
        self.bspace_low_t = bspace_low_t


    def im_amplitude_low_t(self, t, p):
        """Calculates imaginary part of extrapolated amplitude"""
        data = self.dataset.data
        sigma = self.dataset.sigma
        parameters = [i for i in p]
        a0 = sigma / (4 * sqrt(pi * k_norm) )
        a1 = sqrt(data[0].ds - self.model.amplitude(data[0].t, p).real ** 2)
        b0 = (1./(data[0].t)) * log(a0 / a1)

        result = a0 * exp(-1. * b0 * t)
        return result

    def integral(self, b, p, i):
        q1, q2 = sqrt(i.lower), sqrt(i.upper)
        return sqrt((i.ds - self.model.amplitude(i.t, p).real ** 2)) * (q2 * j1(b * q2 / k_fm) -  q1 * j1(b * q1 / k_fm))


    def __call__(self ,b, p):
        """Calculates gamma(b) using data points"""
        data = self.dataset.data

        # Taking into account low-t extrapolation
        extrapolation1 = self.bspace_low_t(b, p, (0, data[0].lower ** 0.5))


        # High t-contribution
        extrapolation2 = self.model.imag_gamma(b, p, (data[-1].upper ** 0.5, float("inf")))


        # Integrated values from data
        gamma_data = sum(self.integral(b, p, i) for i in data) 

        # All contributions
        result = extrapolation1 + extrapolation2 + (gamma_data * k_fm / sqrt(pi * k_norm) / b)
        return result


    @classmethod
    def _func(klass, model, dataset):
        self = klass(model, dataset)
        return lambda b: self(b, dataset.parameters)


    @classmethod
    def tf1(klass, model, dataset, bmin, bmax):
        """Creates \Gamma(b) functor uses data !"""
        g = klass(model, dataset)
        gamma = ROOT.TF1('#Gamma(b)', lambda x, p: g(x[0], dataset.parameters), bmin, bmax, len(dataset.parameters))

        for i, p in enumerate(dataset.parameters):
            gamma.SetParameter(i, p)

        gamma.SetLineColor(46)
        gamma.GetXaxis().SetTitle('b\t,fm')
        gamma.GetYaxis().SetTitle('#Gamma')
        return gamma 

    @classmethod
    def values(klass, model, dataset, index):
        gamma_lambda = klass._func(model, dataset)
        gamma = map(gamma_lambda, index)
        return gamma