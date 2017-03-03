#!/usr/bin/python2.7

from ROOT import *
from ComputeGamma import *
from Formulas import getRealGammaError, getRealGamma
import os
import progressbar



def getGraph(lst):
    graph = TGraphErrors()
    graph.SetName('graph')
    graph.SetTitle('graph')
    [graph.SetPoint(i, i * 3./100, p[0]) for i, p in enumerate(lst)]
    [graph.SetPointError(i, 0, p[1]) for i, p in enumerate(lst)]
    return graph



class ErrorEstimator(object):
    def __init__(self, ptype, energy, sigma, rho, dsigma, drho, nmc):
        super(ErrorEstimator, self).__init__()
        self.process = ptype
        self.energy = energy
        self.sigma = sigma
        self.rho = rho  
        self.dsigma = dsigma
        self.drho = drho
        self.nmc = nmc

        self.gausf = TF1('gaussFunction','gaus', 0, 1.4)
        self.gausf.SetParameter(0, 1)
        self.gausf.SetParameter(1, 0.1)
        self.gausf.SetParameter(2, 1)
        self.gamma_estimator = ComputeGamma(self.process, self.energy, self.sigma, self.rho)


    def average_and_deviation(self, data):
        pivot = data[0]  # chosing some random pivot
        hist = TH1F('hist', 'gamma distribution', 100, pivot - 0.2, pivot - 0.2)
        map(hist.Fill, data)
        hist.Fit(self.gausf,'0qr')
        mu, sigma = self.gausf.GetParameter(1), self.gausf.GetParameter(2)
        return mu, sigma


    def generate_mc_data(self):
        bar = progressbar.ProgressBar()
        print 'Generating the sample of gamma points'
        mc = [self.gamma_estimator.generate_mc_gamma(100, i, self.dsigma) for i in bar(range(self.nmc))]
        # os.remove(c.gamma_fitter.par_file_name)
        return mc


    def estimate_deviations(self, mc):
        print 'Calculating mean, and deviation for all gamma values'
        bar = progressbar.ProgressBar()
        return [self.average_and_deviation(p) for p in bar(zip(*mc))]


    def save_errors(self, mc_av_and_deviation):
        f = open(self.process + '-' + str(self.energy) + '-errors.dat','w')
        for i in mc_av_and_deviation: f.write('%f\t%f\n' % i)


    def main(self):
        # Calculate experimental values of gamma
        gamma_estimator = ComputeGamma(self.process, self.energy, self.sigma, self.rho)
        gamma_estimator.performComputations()

        # Generate monte-carlo data
        mc = self.generate_mc_data()
        mc_av_and_deviation = self.estimate_deviations(mc)

        self.draw_results(mc_av_and_deviation, gamma_estimator)

        # TODO: why do we need two outputs? 
        # Save the results
        self.save_result(gamma_estimator, mc_av_and_deviation)

        # Save the results
        self.save_errors(mc_av_and_deviation)


        return zip(*mc_av_and_deviation)


    def draw_results(self, mc_av_and_deviation, gamma_estimator):
        canvas_point = gamma_estimator.gamma_fitter.canvas
        canvas_point.cd(2)

        # Rearrange all necessary quantities
        mu, sigma  = zip(*mc_av_and_deviation)
        true_gamma = map(gamma_estimator.get_gamma, gamma_estimator.impact_range(self.nmc))
        fake_sigma = [0 for i in mu]

        # Remove this
        gamma_vs_errors = zip(true_gamma, sigma)
        average_vs_zero = zip(mu, fake_sigma)


        # Drow both graphs
        #

        final_result = getGraph(gamma_vs_errors)
        final_result.SetLineColor(37)
        final_result.SetMarkerColor(37)

        average_mc = getGraph(average_vs_zero)
        average_mc.SetLineColor(46)
        average_mc.SetMarkerColor(46)

        final_result.Draw('same')
        average_mc.Draw('same')
        canvas_point.Update()
        canvas_point.SaveAs(str(self.energy) + self.process + '.eps')
        raw_input('pease enter any key ...')
        return gamma_estimator


    def save_result(self, gamma_estimator, averaged_mc_gammas):
        format = '%f\t%f\t%f\t%f\t%f\n'
        parameters = gamma_estimator.read_parameters()
        data = (
                    gamma_estimator.get_gamma(1e-5), 
                    averaged_mc_gammas[0][1], 
                    getRealGamma([0], parameters),
                    getRealGammaError([0], parameters, gamma_estimator.gamma_fitter.covariance, self.dsigma, self.drho), 
                    self.energy
                )

        with open('gamma_at_zero_errors.txt', 'a') as ofile:
            ofile.write(format % data)
