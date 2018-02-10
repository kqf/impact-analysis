#!/usr/bin/python2.7

import os
import progressbar
import json

from ROOT import *
from model import RealGammaEstimator, ImageGammaEstimator
from utils import impact_range

import errors.real as err_real 
import errors.imag as err_imag

from datafit import DataFit
import pandas as pd

from impact.parametrization.symbolic import Symbolic


class ImpactAnalysis(object):
    with open('config/impactanalysis.json') as f:
        conf = json.load(f)

    def __init__(self, model=Symbolic()):
        super(ImpactAnalysis, self).__init__()
        self.model = model
        self.ofile = self.conf['ofile']
        self.points_pref = self.conf['points_pref']


    def run(self, dataset):
       # Calculate experimental values of gamma
        gamma_fitter = DataFit(self.model)

        # Setup the parameters and covariance
        gamma_fitter.fit(dataset)
        pipeline = [
            RealGammaEstimator(self.model),
            ImageGammaEstimator(self.model),
            err_imag.Error(self.model),
            err_real.Error(self.model)
        ]

        output = pd.DataFrame(index=impact_range())
        output.index.name = 'b'

        for estimator in pipeline:
            estimator.evaluate(dataset, output)
            
        return output
