#!/usr/bin/python2.7

import os
import progressbar
import json

from ROOT import *
from model import RealGammaEstimator, ImageGammaEstimator, RealGammaErrorEstimator
from utils import impact_range
import errors.imag as err_imag

from datafit import DataFit
import pandas as pd

from impact.parametrization.symbolic import Symbolic
from impact.parametrization.numeric import Numeric


class ImpactAnalysis(object):

    def __init__(self, model=Symbolic(), conffile='config/datafit.json'):
        super(ImpactAnalysis, self).__init__()
        self.model = model
        self.conffile = conffile


    def run(self, dataset):
       # Calculate experimental values of gamma
        gamma_fitter = DataFit(self.model, self.conffile)

        # Setup the parameters and covariance
        gamma_fitter.fit(dataset)
        pipeline = [
            RealGammaEstimator(self.model),
            ImageGammaEstimator(self.model),
            err_imag.Error(self.model),
            RealGammaErrorEstimator(self.model)
        ]

        output = pd.DataFrame(index=impact_range())
        output.index.name = 'b'

        for estimator in pipeline:
            estimator.evaluate(dataset, output)
            
        return output
