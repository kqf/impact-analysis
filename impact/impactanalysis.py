#!/usr/bin/python2.7

import os
import progressbar
import json

from ROOT import *
from model import real_gamma
from errors import RealPartErrorEvaluator
from errors_image import ImageError
from gammacomputor import ComputeGamma


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
        
        # TODO: Read the data here
        self.ofilename = ptype + '-' + str(energy)
        self.energy = energy
        self.dsigma = dsigma
        self.drho = drho

        # TODO: Why does this parameter is not important
        self.nmc = nmc
        self.mode = mode
        self.gamma_estimator = ComputeGamma(infile, ptype, self.energy, sigma, rho)
        self.imag_errors = ImageError(self.gamma_estimator, nmc, dsigma)

        self.points_pref = self.conf['points_pref']
        self.ofile       = self.conf['ofile']
        self.imgfile     = self.conf['imgfile']


    def save_points_vs_errors(self, mc_av_and_deviation, pref):
        with open(self.points_pref % (pref + '-' + self.ofilename), 'w') as f:
            for i in mc_av_and_deviation:
                f.write('%f\t%f\n' % i)


    def run(self):
        # Calculate experimental values of gamma
        self.gamma_estimator.compute()

        # Generate monte-carlo data
        mc_av_and_deviation = self.imag_errors.evaluate()

        result = self.gamma_estimator.\
            gamma_fitter.draw_results(mc_av_and_deviation, self.gamma_estimator,
            self.nmc,self.imgfile % self.ofile, self.mode)

        # Save the results at avery gamma point
        self.save_result(result[0][0],
                         result[0][1])

        # Save the results only at zero but add more parameters
        self.save_points_vs_errors(mc_av_and_deviation, 'mc')
        self.save_points_vs_errors(result, 'result')

        return zip(*mc_av_and_deviation)




    def save_result(self, gammazero, sigmazero):
        format = '%f\t%f\t%f\t%f\t%f\n'
        parameters = self.gamma_estimator.read_parameters()
        real_gamma_error = RealPartErrorEvaluator(self.gamma_estimator.gamma_fitter.covariance, self.dsigma, self.drho)
        data = (
                    gammazero,
                    sigmazero, 
                    real_gamma(0, parameters),
                    real_gamma_error.breal_error(0, parameters), 
                    self.energy
                )

        with open(self.ofile, 'a') as ofile:
            ofile.write(format % data)
