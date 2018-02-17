#!/usr/bin/python2.7

import os
import progressbar
import json
import pandas as pd

from impact.estimators import RealGammaEstimator, ImagGammaEstimator, RealGammaErrorEstimator, ImagGammaErrorEstimator
from impact.utils import impact_range
from impact.datafit import DataFit

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
            RealGammaEstimator(self.model, outname="real_gamma"),
            ImagGammaEstimator(self.model, outname="imag_gamma"),
            RealGammaErrorEstimator(self.model, outname="real_gamma_error"),
            ImagGammaErrorEstimator(self.model, outname="imag_gamma_error")
        ]

        output = pd.DataFrame(index=impact_range())
        output.index.name = 'b'

        for estimator in pipeline:
            estimator.evaluate(dataset, output)
            
        return output
