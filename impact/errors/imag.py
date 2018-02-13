
import ROOT
import progressbar
import random as rnd

from impact.model import Approx
from impact.datapoint import DataPoint, DataSet
import copy


class Error(object):
    def __init__(self, model, outname="image_error", outaverage="average_impact_amplitude"):
        super(Error, self).__init__()
        # This is value should be fixed
        self.model = model
        self.mcsize = 100
        self.outname = outname
        self.outaverage = outaverage

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


    def generrate_diff_cs(self, p, dataset):
        cs, re2 = -1, self.model.amplitude(p.t, dataset.parameters).real ** 2
        while (cs - re2) < 0: # Generate only positive cs
            cs = rnd.gauss(p.ds, p.err)
        return cs


    def generate_mc_points(self, dataset):
        return [
            DataPoint(
                p.t,
                self.generrate_diff_cs(p, dataset),
                p.err,
                p.lower,
                p.upper
            ) 
        for p in dataset.data]
        

    # TODO: Nb introduce additional parameter for nsample points
    def generate_mc_gamma(self, dataset, index):
        """Calculates points from b = 0 to b = 3"""
        ## Read parameters for real data approximation
        mcPoints = self.generate_mc_points(dataset)

        ## Get gamma approximator using data points, generate random cross section value
        new_sigma = rnd.gauss(dataset.sigma, dataset.dsigma)

        generated_dataset = copy.deepcopy(dataset)
        generated_dataset.sigma = new_sigma
        generated_dataset._data = mcPoints


        ## Calculate values of Gamma function for the mc
        ## Index should be applied here
        return Approx.values(
            self.model,
            generated_dataset,
            index
        )


    def generate_mc_data(self, dataset, index):
        bar = progressbar.ProgressBar()
        print 'Generating the sample of gamma points'
        # NB: First parameter in this list may differ from range
        mc = [self.generate_mc_gamma(dataset, index) for i in bar(range(self.mcsize))]
        return mc


    def evaluate(self, dataset, output):
        mc = self.generate_mc_data(dataset, output.index)
        mc_av_and_deviation = self.estimate_deviations(mc)
        average, deviation = zip(*mc_av_and_deviation)

        output[self.outaverage] = average
        output[self.outname] = deviation

        # output[self.outaverage] = output.index * 0
        # output[self.outname] = output.index * 0
        return average, deviation



