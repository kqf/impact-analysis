#!/usr/bin/python2.7

import copy
import random as rnd
from math import exp, log, sqrt

# import impact.utils as ut
import numpy as np
import pandas as pd
import progressbar
import ROOT
from constants import k_fm
from impact.datapoint import DataPoint
from impact.utils import hankel_transform
from impact.datafit import DataFit
from scipy.special import j1


class RealGammaEstimator(object):

    def __init__(self, model, outname="real_gamma"):
        self.model = model
        self.outname = outname

    def evaluate(self, dataset, output):
        output[self.outname] = map(lambda x: self.model.real_gamma(
            x, dataset.parameters) * self.model.h_norm(), output.index)
        return output[self.outname].values


class RealGammaErrorEstimator(object):
    def __init__(self, model, outname='real_gamma_error'):
        super(RealGammaErrorEstimator, self).__init__()
        self.outname = outname
        self.model = model

        @hankel_transform
        def evaluate_(x, dataset):
            return self.model.treal_error(x, dataset) * self.model.h_norm()
        self.evaluate_ = evaluate_

    def evaluate(self, dataset, output):
        output[self.outname] = map(
            lambda x: self.evaluate_(x, dataset),
            output.index
        )
        return output[self.outname]


class ImagGammaParametrization(object):

    def __init__(self, model, outname="imag_gamma_param"):
        self.model = model
        self.outname = outname

    def evaluate(self, dataset, output):
        output[self.outname] = map(lambda x: self.model.imag_gamma(
            x, dataset.parameters) * self.model.h_norm(), output.index)
        return output[self.outname].values


class ImagGammaParametrizationErrorEstimator(object):
    def __init__(self, model, outname='real_gamma_error'):
        super(ImagGammaParametrizationErrorEstimator, self).__init__()
        self.outname = outname
        self.model = model

        @hankel_transform
        def evaluate_(x, dataset):
            return self.model.timag_error(x, dataset) * self.model.h_norm()
        self.evaluate_ = evaluate_

    def evaluate(self, dataset, output):
        output[self.outname] = map(
            lambda x: self.evaluate_(x, dataset),
            output.index
        )
        return output[self.outname]


def sign(x):
    result = np.sign(x)
    if result == 0.0:
        return 1.
    return result


class ImagGammaEstimator(object):
    def __init__(self, model, outname="imag_gamma"):
        self.model = model
        self.outname = outname
        # Define these functions in a very strange
        # way to avoid code replication

        @hankel_transform
        def bspace_low_t(b, p):
            return self._im_amplitude_low_t(b, p) * self.model.h_norm()
        self.bspace_low_t = bspace_low_t

    def evaluate(self, dataset, output):
        def f(x):
            return self._gamma(x, dataset)
        output[self.outname] = map(f, output.index)
        return output[self.outname].values

    def _im_amplitude_low_t(self, t, dataset):
        """Calculates imaginary part of extrapolated amplitude"""
        data = dataset.hadron_data
        sigma = dataset.sigma
        a0 = sigma / self.model.sigma_norm()
        a1 = sqrt(data[0].ds / self.model.dsigdt_norm() - self.model.amplitude(
            data[0].t,
            dataset.parameters).real ** 2
        )
        b0 = (1. / (data[0].t)) * log(a0 / a1)
        result = a0 * exp(-1. * b0 * t)
        return result

    def _integral(self, b, p, i):
        q1, q2 = sqrt(i.lower), sqrt(i.upper)
        diff = (
            i.ds / self.model.dsigdt_norm() -
            self.model.amplitude(i.t, p).real ** 2
        )
        weight = (q2 * j1(b * q2 / k_fm) - q1 * j1(b * q1 / k_fm))
        try:
            image_a = self.model.amplitude(i.t, p).imag
            return sqrt(diff) * weight * sign(image_a)
        except ValueError:
            return 0

    def _gamma(self, b, dataset):
        data = dataset.hadron_data

        # Taking into account low-t extrapolation
        extrapolation1 = self.bspace_low_t(
            b, dataset, (0, data[0].lower ** 0.5))

        # High t-contribution
        data_range = data[-1].upper ** 0.5, float("inf")
        extrapolation2 = self.model.imag_gamma(
            b, dataset.parameters, data_range)
        extrapolation2 *= self.model.h_norm()

        # Integrated values from data
        gamma_data = sum(self._integral(b, dataset.parameters, i)
                         for i in data)

        gamma_data_scaled = gamma_data * k_fm / b * self.model.hdata_norm()
        # All contributions
        result = extrapolation1 + extrapolation2 + gamma_data_scaled
        return result


class ImagGammaErrorEstimator(object):
    def __init__(self, model, outname="imag_gamma_error",
                 outaverage="average_impact_amplitude",
                 resolution=100, n_sigma=2.):
        super(ImagGammaErrorEstimator, self).__init__()
        # This is value should be fixed
        self.generator = GammaGeneratorMC(model, n_sigma=n_sigma)
        self.gamma_resolution = resolution
        self.outname = outname
        self.outaverage = outaverage

    def _average_and_deviation(self, data):
        indata = np.array(data)

        hist = ROOT.TH1F(
            'gamma_hist',
            'MC generated im #Gamma(b) points at same b; Im #Gamma(b); counts',
            self.gamma_resolution,
            np.min(indata) * 0.98,
            np.max(indata) * 1.02
        )

        map(hist.Fill, indata)

        # canvas = ut.canvas('test')

        gausf = ROOT.TF1('gaussFunction', 'gaus')
        gausf.SetParameter(0, 1)
        gausf.SetParameter(1, np.std(indata))
        gausf.SetParameter(2, np.mean(indata))
        hist.Fit(gausf, 'q')
        # hist.Draw()
        # canvas.Update()
        # canvas.SaveAs("test-{0}.pdf".format(rnd.randint(0, 100000)))
        mu, sigma = gausf.GetParameter(1), gausf.GetParameter(2)
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
        # output[self.outaverage] = output.index * 0.0
        # output[self.outname] = output.index * 0.0
        # return output.index * 0.0, output.index * 0.0


class GammaGeneratorMC(object):

    def __init__(self, model, mcsize=100, n_sigma=2):
        super(GammaGeneratorMC, self).__init__()
        self.imag_gamma = ImagGammaEstimator(model)
        self.model = model
        self.mcsize = mcsize
        self.n_sigma = n_sigma

    def _generrate_diff_cs(self, p, dataset):
        cs, re2 = -1, self.model.amplitude(p.t, dataset.parameters).real ** 2
        # Generate only positive cs
        while (cs / self.model.dsigdt_norm() - re2) < 0:
            cs = rnd.gauss(p.ds, self.n_sigma * p.err)
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
            for p in dataset.hadron_data]

    def _gamma(self, dataset, index):
        # Read parameters for real data approximation
        mcPoints = self._datapoints(dataset)

        # Get gamma approximator using data points,
        # generate random cross section value
        new_sigma = rnd.gauss(dataset.sigma, dataset.dsigma)

        generated_dataset = copy.deepcopy(dataset)
        generated_dataset.sigma = new_sigma
        generated_dataset._data = mcPoints

        # Calculate values of Gamma function for the mc
        # Index should be applied here
        return self.imag_gamma.evaluate(
            generated_dataset,
            pd.DataFrame(index=index)
        )

    def evaluate(self, dataset, index):
        bar = progressbar.ProgressBar()
        print 'Generating the sample of gamma points'
        # NB: First parameter in this list may differ from range
        mc = [self._gamma(dataset, index) for i in bar(range(self.mcsize))]
        return mc


class GInelEstimator(object):
    """GInelEstimator

        The $G_{inel}$ function is defined as:
        $$G_{inel} = Im H(s, b) - |H(s, b)|^2 $$,
        where $H(s, b) = 2i\Gamma(s, b)$
    """

    def __init__(self, inreal, inimag, outname):
        super(GInelEstimator, self).__init__()
        self.inreal = inreal
        self.inimag = inimag
        self.outname = outname

    def evaluate(self, dataset, output):
        hreal = 0.5 * output[self.inreal].values
        himag = 0.5 * output[self.inimag].values
        output[self.outname] = himag - hreal ** 2 - himag ** 2
        return output[self.outname].values


class GInelErrorEstimator(object):
    """GInelErrorEstimator

        The $G_{inel}$ function is defined as:
        $$\Delta $G_{inel}$  = (1 - 2 * H(s, b)) ^2 \Delta H ^ 2 $$,
        where $H(s, b) = 2i\Gamma(s, b)$

    """

    def __init__(self, inimag, inreal, inimagerr, inrealerr, outname):
        super(GInelErrorEstimator, self).__init__()
        self.inimag = inimag
        self.inreal = inreal
        self.inimagerr = inimagerr
        self.inrealerr = inrealerr
        self.outname = outname

    def evaluate(self, dataset, output):
        hreal = 0.5 * output[self.inreal].values
        himag = 0.5 * output[self.inimag].values
        dhreal = 0.5 * output[self.inrealerr].values
        dhimag = 0.5 * output[self.inimagerr].values

        output[self.outname] = np.sqrt(
            (dhimag ** 2) * (1 - 2 * himag) ** 2 +
            (dhreal ** 2) * (2 * hreal) ** 2
        )
        return output[self.outname].values


class DataGenerator(object):
    def __init__(self, n_iterations, n_sigma):
        super(DataGenerator, self).__init__()
        self.n_iterations = n_iterations
        self.n_sigma = n_sigma

    def evaluate(self, dataset):
        bar = progressbar.ProgressBar()
        datasets = [self._generate(dataset)
                    for _ in bar(range(self.n_iterations))]
        return datasets

    def _generate(self, dataset):
        # Read parameters for real data approximation
        mc_points = [
            DataPoint(
                p.t,
                rnd.gauss(p.ds, self.n_sigma * p.err),
                p.err,
                p.lower,
                p.upper
            )
            for p in dataset.hadron_data]

        new_sigma = rnd.gauss(dataset.sigma, dataset.dsigma)
        generated_dataset = copy.deepcopy(dataset)
        generated_dataset.sigma = new_sigma
        generated_dataset._data = mc_points
        return generated_dataset


class AlnternativeErrorEstimator(object):

    def __init__(self, conffile, model, n_iterations, n_sigma,
                 outname_mean, outname_std):
        self.conffile = conffile
        self.model = model
        self.outname_std = outname_std
        self.outname_mean = outname_mean
        self.generator = DataGenerator(n_iterations, n_sigma)

    def evaluate(self, dataset, output):
        generated = self.generator.evaluate(dataset)
        bar = progressbar.ProgressBar()
        mc_gamma = np.array([
            self._evaluate_signle(fake_dataset, output.index)
            for fake_dataset in bar(generated)
        ])
        pd.DataFrame(mc_gamma).to_csv("mc-gamma-{}.csv".format(
            self.model.name), index=False)
        output[self.outname_mean] = np.mean(mc_gamma, axis=0)
        output[self.outname_std] = np.std(mc_gamma, axis=0)

    def _evaluate_signle(self, dataset, index):
        local_fitter = DataFit(self.model, self.conffile, mute=True)
        local_fitter.fit(dataset)

        output = pd.DataFrame(index=index)
        ImagGammaParametrization(self.model, outname="imag_gamma").evaluate(
            dataset,
            output
        )
        return output["imag_gamma"].values
