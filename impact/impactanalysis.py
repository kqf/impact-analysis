#!/usr/bin/python2.7

import os
import progressbar
import json

from ROOT import *
from model import real_gamma

import errors.real as err_real 
import errors.imag as err_imag

from datapoint import DataPoint
from datafit import DataFit
import model 


class ImpactAnalysis(object):
    with open('config/impactanalysis.json') as f:
        conf = json.load(f)

    def __init__(self, infile, ptype, energy, sigma, rho, dsigma, drho, mode= 's'):
        super(ImpactAnalysis, self).__init__()

        self.name = ptype + '-' + str(energy)
        self.energy = energy
        self.dsigma = dsigma
        self.drho = drho

        self.mode = mode
        self.points_pref = self.conf['points_pref']
        self.ofile = self.conf['ofile']

        # TODO: Move the real gamma error here
        self.data = DataPoint.read(energy, ptype, infile)
        self.gamma_fitter = DataFit(self.data, self.name, ptype, energy, sigma, rho, mode)
        self.imag_gamma_error = err_imag.Error(self.data, sigma, dsigma)


    def save_points_vs_errors(self, gamma, gamma_error, pref):
        with open(self.points_pref % (pref + '-' + self.name), 'w') as f:
            for i in zip(gamma, gamma_error):
                f.write('%f\t%f\n' % i)


    def run(self):
        # Calculate experimental values of gamma
        parameters, covariance = self.gamma_fitter.fit()

        # Estimate average values and  errors from monte-carlo data
        average, sigma = self.imag_gamma_error.evaluate(parameters)

        # Compare mc value with real and return the real values
        gamma = self.gamma_fitter.compare_results(average, sigma, parameters)

        # Save gamma at zero
        self.save_result(parameters, covariance, gamma[0], sigma[0])

        # Save the results not only at zero but add more parameters
        self.save_points_vs_errors(average, sigma, 'mc')
        self.save_points_vs_errors(gamma, sigma,'result')
        return gamma, sigma




    def save_result(self, parameters, covariance, gammazero, sigmazero):
        format = '%f\t%f\t%f\t%f\t%f\n'
        # TODO: move covariance and parameters to the same place in 
        #       in the error estimation
        real_gamma_error = err_real.Error(covariance, self.dsigma, self.drho)
        data = (
                    gammazero,
                    sigmazero, 
                    real_gamma(0, parameters),
                    real_gamma_error.evaluate(0, parameters), 
                    self.energy
                )

        with open(self.ofile, 'a') as ofile:
            ofile.write(format % data)
