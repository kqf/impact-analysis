
import ROOT
import progressbar
import random as rnd

import impact.model as model
from impact.datapoint import DataPoint


class Error(object):
    def __init__(self, data, sigma, dsigma):
        super(Error, self).__init__()
        self.data = data
        self.sigma = sigma
        self.dsigma = dsigma
        # This is value should be fixed
        self.mcsize = 100

        self.gausf = ROOT.TF1('gaussFunction','gaus', 0, 1.4)
        self.gausf.SetParameter(0, 1)
        self.gausf.SetParameter(1, 0.1)
        self.gausf.SetParameter(2, 1)       

    def average_and_deviation(self, data):
        pivot = data[0]  # chosing some random pivot
        hist = ROOT.TH1F('hist', 'gamma distribution', 100, pivot - 0.2, pivot - 0.2)
        map(hist.Fill, data)
        hist.Fit(self.gausf,'0qr')
        mu, sigma = self.gausf.GetParameter(1), self.gausf.GetParameter(2)
        return mu, sigma


    def estimate_deviations(self, mc):
        print 'Calculating mean, and deviation for all gamma values'
        bar = progressbar.ProgressBar()
        return [self.average_and_deviation(p) for p in bar(zip(*mc))]


    def generate_mc_points(self):
        """Generates MC data"""
        return [DataPoint(p.t, rnd.gauss(p.ds, p.err), p.err, p.lower, p.upper) for p in self.data]


    # TODO: Nb introduce additional parameter for nsample points
    def generate_mc_gamma(self, parameters):
        """Calculates points from b = 0 to b = 3"""
        ## Read parameters for real data approximation
        mcPoints = self.generate_mc_points()

        ## Get gamma approximator using data points, generate random cross section value
        mc_gamma = model.GammaApproximation(mcPoints, rnd.gauss(self.sigma, self.dsigma))
        gamma = lambda x: mc_gamma(x, parameters)

        ## Calculate values of Gamma function for the mc
        return map(gamma, model.impact_range())


    def generate_mc_data(self, parameters):
        bar = progressbar.ProgressBar()
        print 'Generating the sample of gamma points'
        # NB: First parameter in this list may differ from range
        mc = [self.generate_mc_gamma(parameters) for i in bar(range(self.mcsize))]
        return mc


    def evaluate(self, parameters):
        mc = self.generate_mc_data(parameters)
        mc_av_and_deviation = self.estimate_deviations(mc)
        return zip(*mc_av_and_deviation)

