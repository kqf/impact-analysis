
import ROOT
import progressbar
import random as rnd
from datapoint import DataPoint
import model
from model import GammaApproximation 


class ImageError(object):
    def __init__(self, data, nmc, sigma, dsigma):
        super(ImageError, self).__init__()
        self.data = data
        self.nmc = nmc 
        self.sigma = sigma
        self.dsigma = dsigma
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

    def evaluate(self, parameters):
        mc = self.generate_mc_data(parameters)
        mc_av_and_deviation = self.estimate_deviations(mc)
        return mc_av_and_deviation


    def generate_mc_data(self, parameters):
        bar = progressbar.ProgressBar()
        print 'Generating the sample of gamma points'
        # TODO: Improve names in this class
        mc = [self.generate_mc_gamma(100, parameters, self.dsigma) for i in bar(range(self.nmc))]
        return mc

    def generate_mc_gamma(self, npoints, parameters, dsigma):
        """Writes points from b = 0 to b = 3 to file"""
        ## Read parameters for real data approximation
        mcPoints = self.generate_mc_points()

        ## Get gamma approximator using data points, generate random cross section value
        mc_gamma = GammaApproximation(mcPoints, rnd.gauss(self.sigma, dsigma))
        gamma = lambda x: mc_gamma(x, parameters)

        ## Calculate values of Gamma function for the mc
        return map(gamma, model.impact_range(npoints, npoints/3.0))

    def generate_mc_points(self):
        """Generates MC data"""
        return [DataPoint(p.t, rnd.gauss(p.ds, p.err), p.err, p.lower, p.upper) for p in self.data]

