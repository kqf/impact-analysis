#!/usr/bin/python2.7

import os
import progressbar
import json

from ROOT import *
from model import RealGammaEstimator, impact_range, ImageGammaEstimator

import errors.real as err_real 
import errors.imag as err_imag

from datapoint import DataSet
from datafit import DataFit
import model 
import pandas as pd


class ImpactAnalysis(object):
    with open('config/impactanalysis.json') as f:
        conf = json.load(f)

    def __init__(self, mode= 's'):
        super(ImpactAnalysis, self).__init__()
        self.ofile = self.conf['ofile']
        self.points_pref = self.conf['points_pref']
        self.mode = mode


    def run(self, dataset):
        gamma_fitter = DataFit(self.mode)
        # Calculate experimental values of gamma

        # Setup the parameters and covariance
        gamma_fitter.fit(dataset)
        pipeline = [
            RealGammaEstimator(),
            ImageGammaEstimator(),
            err_imag.Error(),
            err_real.Error()
        ]

        output = pd.DataFrame(index=impact_range())
        for estimator in pipeline:
            estimator.evaluate(dataset, output)

        print output[:20]
        return output['image_gamma'].values, output['image_error'] 
