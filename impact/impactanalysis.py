#!/usr/bin/python2.7

import pandas as pd

from impact.utils import impact_range
from impact.datafit import DataFit

from impact.parametrization.symbolic import Standard
from impact.estimators import (
    RealGammaEstimator,
    ImagGammaEstimator,
    RealGammaErrorEstimator,
    ImagGammaErrorEstimator,
    GInelEstimator,
    GInelErrorEstimator,
    ImagGammaParametrization,
    ImagGammaParametrizationErrorEstimator,
    AlnternativeErrorEstimator
)


class ImpactAnalysis(object):

    def __init__(self, model=Standard(),
                 conffile='config/standard.json', n_sigma=2.):
        super(ImpactAnalysis, self).__init__()
        self.model = model
        self.conffile = conffile
        self.n_sigma = n_sigma

    def run(self, dataset):
        # Calculate experimental values of gamma
        gamma_fitter = DataFit(self.model, self.conffile)

        # Setup the parameters and covariance
        gamma_fitter.fit(dataset)
        pipeline = [
            RealGammaEstimator(self.model, outname="real_gamma"),
            ImagGammaEstimator(self.model, outname="imag_gamma"),
            RealGammaErrorEstimator(self.model, outname="real_gamma_error"),
            ImagGammaErrorEstimator(self.model, outname="imag_gamma_error",
                                    n_sigma=self.n_sigma),
            GInelEstimator(inreal="real_gamma",
                           inimag="imag_gamma", outname="g_inel"),
            GInelErrorEstimator(
                inreal="real_gamma",
                inimag="imag_gamma",
                inimagerr="imag_gamma_error",
                inrealerr="real_gamma_error",
                outname="g_inel_error"
            ),
            ImagGammaParametrization(self.model, outname="imag_gamma_param"),
            ImagGammaParametrizationErrorEstimator(
                self.model,
                outname="imag_gamma_error_param"),
            GInelEstimator(inreal="real_gamma",
                           inimag="imag_gamma_param", outname="g_inel_param"),
            GInelErrorEstimator(
                inreal="real_gamma",
                inimag="imag_gamma_param",
                inimagerr="imag_gamma_error_param",
                inrealerr="real_gamma_error",
                outname="g_inel_error_param"
            ),
            AlnternativeErrorEstimator(
                conffile=self.conffile,
                model=self.model,
                n_iterations=30,
                n_sigma=2,
                outname_mean="imag_gamma_gen",
                outname_std="imag_gamma_error_gen",
            )
        ]

        output = pd.DataFrame(index=impact_range())
        output.index.name = 'b'

        for estimator in pipeline:
            estimator.evaluate(dataset, output)

        output["im_h"] = output["imag_gamma"].values * 0.5
        output["re_h"] = -output["real_gamma"].values * 0.5
        output["im_h_error"] = output["imag_gamma_error"].values * 0.5
        output["re_h_error"] = output["real_gamma_error"].values * 0.5

        output["im_h_param"] = output["imag_gamma_param"].values * 0.5
        output["im_h_error_param"] = output["imag_gamma_error_param"].values * 0.5

        output["im_h_gen"] = output["imag_gamma_gen"].values * 0.5
        output["im_h_error_gen"] = output["imag_gamma_error_gen"].values * 0.5
        return output
