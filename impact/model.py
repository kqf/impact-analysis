#!/usr/bin/python2.7

from math  import exp, pi, fabs, sqrt, log
from scipy import integrate
from scipy.special import j0, j1
import pandas as pd
import progressbar
import random as rnd
import copy

import ROOT

from constants import k_fm, k_norm
from impact.utils import hankel_transform, impact_range
from impact.datapoint import DataPoint
from parametrization.symbolic import Symbolic, SymbolicUpdated


class RealGammaEstimator(object):

    def __init__(self, model, outname="real_gamma"):
        self.model = model
        self.outname = outname

    def evaluate(self, dataset, output):
        output[self.outname] = map(lambda x: self.model.real_gamma(x, dataset.parameters), output.index)
        return output[self.outname].values


class RealGammaErrorEstimator(object):
    def __init__(self, model, outname='real_gamma_error'):
        super(RealGammaErrorEstimator, self).__init__()
        self.outname = outname
        self.model = model

        @hankel_transform
        def evaluate_(x, dataset):
            return self.model.treal_error(x, dataset)
        self.evaluate_ = evaluate_


    def evaluate(self, dataset, output):
        output[self.outname] = map(
            lambda x: self.evaluate_(x, dataset),
            output.index
        )
        return output[self.outname]
        
        
class ImageGammaEstimator(object):
    def __init__(self, model, outname="image_gamma"):
        self.model = model
        self.outname = outname
        # Define these functions in a very strange way to avoid code replication
        @hankel_transform
        def bspace_low_t(b, p):
            return self._im_amplitude_low_t(b, p)
        self.bspace_low_t = bspace_low_t


    def evaluate(self, dataset, output):
        f = lambda x: self._gamma(x, dataset)
        output[self.outname] = map(f, output.index)
        return output[self.outname].values


    def _im_amplitude_low_t(self, t, dataset):
        """Calculates imaginary part of extrapolated amplitude"""
        data = dataset.data
        sigma = dataset.sigma
        a0 = sigma / (4 * sqrt(pi * k_norm) )
        a1 = sqrt(data[0].ds - self.model.amplitude(data[0].t, dataset.parameters).real ** 2)
        b0 = (1./(data[0].t)) * log(a0 / a1)

        result = a0 * exp(-1. * b0 * t)
        return result

    def _integral(self, b, p, i):
        q1, q2 = sqrt(i.lower), sqrt(i.upper)
        return sqrt((i.ds - self.model.amplitude(i.t, p).real ** 2)) * (q2 * j1(b * q2 / k_fm) -  q1 * j1(b * q1 / k_fm))


    def _gamma(self, b, dataset):
        data = dataset.data

        # Taking into account low-t extrapolation
        extrapolation1 = self.bspace_low_t(b, dataset, (0, data[0].lower ** 0.5))

        # High t-contribution
        extrapolation2 = self.model.imag_gamma(b, dataset.parameters, (data[-1].upper ** 0.5, float("inf")))

        # Integrated values from data
        gamma_data = sum(self._integral(b, dataset.parameters, i) for i in data) 

        # All contributions
        result = extrapolation1 + extrapolation2 + (gamma_data * k_fm / sqrt(pi * k_norm) / b)
        return result


class ImageGammaErrorEstimator(object):
    def __init__(self, model, outname="image_error", outaverage="average_impact_amplitude"):
        super(ImageGammaErrorEstimator, self).__init__()
        # This is value should be fixed
        self.generator = GammaGeneratorMC(model)
        self.outname = outname
        self.outaverage = outaverage

        self.gausf = ROOT.TF1('gaussFunction','gaus', 0, 1.4)
        self.gausf.SetParameter(0, 1)
        self.gausf.SetParameter(1, 0.1)
        self.gausf.SetParameter(2, 1)       


    def _average_and_deviation(self, data):
        pivot = data[0]  # chosing some random pivot
        hist = ROOT.TH1F('hist', 'gamma distribution', 100, pivot - 0.2, pivot - 0.2)
        map(hist.Fill, data)
        hist.Fit(self.gausf,'0qr')
        mu, sigma = self.gausf.GetParameter(1), self.gausf.GetParameter(2)
        return mu, sigma


    def _estimate_deviations(self, mc):
        print 'Calculating mean, and deviation for all gamma values'
        bar = progressbar.ProgressBar()
        return [self._average_and_deviation(p) for p in bar(zip(*mc))]


    def evaluate(self, dataset, output):
        mc = self.generator.evaluate(dataset, output.index)
        mc_av_and_deviation = self._estimate_deviations(mc)
        average, deviation = zip(*mc_av_and_deviation)

        output[self.outaverage] = average
        output[self.outname] = deviation
        return average, deviation


class GammaGeneratorMC(object):

    def __init__(self, model, mcsize=100):
        super(GammaGeneratorMC, self).__init__()
        self.image_gamma = ImageGammaEstimator(model)
        self.model = model
        self.mcsize = mcsize


    def _generrate_diff_cs(self, p, dataset):
        cs, re2 = -1, self.model.amplitude(p.t, dataset.parameters).real ** 2
        while (cs - re2) < 0: # Generate only positive cs
            cs = rnd.gauss(p.ds, p.err)
        return cs


    def _datapoints(self, dataset):
        return [
            DataPoint(
                p.t,
                self._generrate_diff_cs(p, dataset),
                p.err,
                p.lower,
                p.upper
            ) 
        for p in dataset.data]
        

    def _gamma(self, dataset, index):
        ## Read parameters for real data approximation
        mcPoints = self._datapoints(dataset)

        ## Get gamma approximator using data points, generate random cross section value
        new_sigma = rnd.gauss(dataset.sigma, dataset.dsigma)

        generated_dataset = copy.deepcopy(dataset)
        generated_dataset.sigma = new_sigma
        generated_dataset._data = mcPoints

        ## Calculate values of Gamma function for the mc
        ## Index should be applied here
        return self.image_gamma.evaluate(generated_dataset, pd.DataFrame(index=index))


    def evaluate(self, dataset, index):
        bar = progressbar.ProgressBar()
        print 'Generating the sample of gamma points'
        # NB: First parameter in this list may differ from range
        mc = [self._gamma(dataset, index) for i in bar(range(self.mcsize))]
        return mc





