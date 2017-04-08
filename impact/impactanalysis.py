#!/usr/bin/python2.7

import os
import progressbar
import json

from ROOT import *
from ComputeGamma import *
from Formulas import getRealGammaError, getRealGamma


def getGraph(lst):
    graph = TGraphErrors()
    graph.SetName('graph')
    graph.SetTitle('graph')
    [graph.SetPoint(i, i * 3./100, p[0]) for i, p in enumerate(lst)]
    [graph.SetPointError(i, 0, p[1]) for i, p in enumerate(lst)]
    return graph


class ImpactAnalysis(object):
    with open('config/impactanalysis.json') as f:
        conf = json.load(f)

    def __init__(self, infile, ptype, energy, sigma, rho, dsigma, drho, nmc, mode= 's'):
        super(ImpactAnalysis, self).__init__()
        self.ofilename = ptype + '-' + str(energy)
        self.energy = energy
        self.dsigma = dsigma
        self.drho = drho
        self.nmc = nmc
        self.mode = mode
        self.gamma_estimator = ComputeGamma(infile, ptype, self.energy, sigma, rho)

        self.gausf = TF1('gaussFunction','gaus', 0, 1.4)
        self.gausf.SetParameter(0, 1)
        self.gausf.SetParameter(1, 0.1)
        self.gausf.SetParameter(2, 1)

        self.points_pref = self.conf['points_pref']
        self.ofile       = self.conf['ofile']
        self.imgfile     = self.conf['imgfile']


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
        return mc


    def estimate_deviations(self, mc):
        print 'Calculating mean, and deviation for all gamma values'
        bar = progressbar.ProgressBar()
        return [self.average_and_deviation(p) for p in bar(zip(*mc))]


    def save_points_vs_errors(self, mc_av_and_deviation, pref):
        with open(self.points_pref % (pref + '-' + self.ofilename), 'w') as f:
            for i in mc_av_and_deviation:
                f.write('%f\t%f\n' % i)


    def run(self):
        # Calculate experimental values of gamma
        self.gamma_estimator.compute()

        # Generate monte-carlo data
        mc = self.generate_mc_data()
        mc_av_and_deviation = self.estimate_deviations(mc)

        result = self.draw_results(mc_av_and_deviation)

        # Save the results at avery gamma point
        self.save_result(result[0][0],
                         result[0][1])

        # Save the results only at zero but add more parameters
        self.save_points_vs_errors(mc_av_and_deviation, 'mc')
        self.save_points_vs_errors(result, 'result')

        return zip(*mc_av_and_deviation)


    def draw_results(self, mc_av_and_deviation):
        canvas_point = self.gamma_estimator.gamma_fitter.canvas
        canvas_point.cd(2)

        # Rearrange all necessary quantities
        mu, sigma  = zip(*mc_av_and_deviation)
        true_gamma = map(self.gamma_estimator.get_gamma, self.gamma_estimator.impact_range(self.nmc))
        fake_sigma = [0 for i in mu]

        # Remove this
        gamma_vs_errors = zip(true_gamma, sigma)
        average_vs_zero = zip(mu, fake_sigma)


        final_result = getGraph(gamma_vs_errors)
        final_result.SetLineColor(37)
        final_result.SetMarkerColor(37)

        average_mc = getGraph(average_vs_zero)
        average_mc.SetLineColor(46)
        average_mc.SetMarkerColor(46)

        final_result.Draw('same')
        average_mc.Draw('same')
        canvas_point.Update()
        canvas_point.SaveAs(self.imgfile % self.ofilename)

        if not 's' in self.mode:
            raw_input('pease enter any key ...')

        return gamma_vs_errors


    def save_result(self, gammazero, sigmazero):
        format = '%f\t%f\t%f\t%f\t%f\n'
        parameters = self.gamma_estimator.read_parameters()
        data = (
                    gammazero,
                    sigmazero, 
                    getRealGamma(0, parameters),
                    getRealGammaError(0, parameters, self.gamma_estimator.gamma_fitter.covariance, self.dsigma, self.drho), 
                    self.energy
                )

        with open(self.ofile, 'a') as ofile:
            ofile.write(format % data)
