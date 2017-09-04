#!/usr/bin/python2.7

import os
import progressbar
import json

from ROOT import *
from model import real_gamma

import errors.real as err_real 
import errors.imag as err_imag

from datapoint import DataReader
from datafit import DataFit
import model 


class ImpactAnalysis(object):
    with open('config/impactanalysis.json') as f:
        conf = json.load(f)

    def __init__(self, infile, ptype, energy, sigma, rho, dsigma, drho, nmc, mode= 's'):
        super(ImpactAnalysis, self).__init__()

        self.name = ptype + '-' + str(energy)
        self.energy = energy
        self.dsigma = dsigma
        self.drho = drho

        self.nmc = nmc
        self.mode = mode
        self.points_pref = self.conf['points_pref']
        self.ofile = self.conf['ofile']

        self.data = DataReader(energy, ptype).read(infile)
        self.gamma_fitter = DataFit(self.data, self.name, ptype, energy, sigma, rho, mode)
        self.imag_gamma_error = err_imag.Error(self.data, nmc, sigma, dsigma)


    def save_points_vs_errors(self, gamma, gamma_error, pref):
        with open(self.points_pref % (pref + '-' + self.name), 'w') as f:
            for i in zip(gamma, gamma_error):
                f.write('%f\t%f\n' % i)


    def run(self):
        # Calculate experimental values of gamma
        _, parameters = self.gamma_fitter.fit()

        # Estimate average values and  errors from monte-carlo data
        average, sigma = self.imag_gamma_error.evaluate(parameters)

        # TODO: Clean the names of the variables
        # TODO: remove self.nmc
        gamma_lambda = model.GammaApproximation.function_for_parameters(self.data, parameters)
        gamma = map(gamma_lambda, model.impact_range(self.nmc, self.nmc / 3.0))

        self.gamma_fitter.draw_results(average, sigma, gamma)

        # Save gamma at zero
        self.save_result(parameters, gamma[0], sigma[0])

        # Save the results not only at zero but add more parameters
        self.save_points_vs_errors(average, sigma, 'mc')
        self.save_points_vs_errors(gamma, sigma,'result')
        return gamma, sigma




    def save_result(self, parameters, gammazero, sigmazero):
        format = '%f\t%f\t%f\t%f\t%f\n'
        real_gamma_error = err_real.Error(self.gamma_fitter.covariance, self.dsigma, self.drho)
        data = (
                    gammazero,
                    sigmazero, 
                    real_gamma(0, parameters),
                    real_gamma_error.evaluate(0, parameters), 
                    self.energy
                )

        with open(self.ofile, 'a') as ofile:
            ofile.write(format % data)
