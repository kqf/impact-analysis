#!/usr/bin/python
import json
import unittest

import pandas as pd
from impact.datapoint import DataSet
from impact.estimators import DataGenerator
from impact.utils import impact_range

from impact.parametrization.symbolic import FullStandard
from impact.parametrization.symbolic import FullThreePlusOne
from impact.parametrization.symbolic import FullThreePlusTwo
from impact.estimators import AlnternativeErrorEstimator


class RunTheAnalysis(unittest.TestCase):

    @unittest.skip("")
    def test_visualize_the_results(self):
        with open("config/input.json") as f:
            data = json.load(f)
        dataset = DataSet(data["data"][-1])

        generator = DataGenerator(n_iterations=10, n_sigma=2)
        generated = generator.evaluate(dataset)
        for dataset in generated:
            dataset.draw()

    def test_elaborated_estimation(self):
        with open("config/input.json") as f:
            data = json.load(f)

        models = {
            "config/full-standard.json": FullStandard,
            "config/full-triple-exponent.json": FullThreePlusOne,
            "config/full-triple-exponent-general.json": FullThreePlusTwo,
        }

        output = pd.DataFrame(index=impact_range())
        for config, algo in models.iteritems():
            estimator = AlnternativeErrorEstimator(
                config, algo(), 10, 2,
                "generated_im_gamma",
                "generated_im_gamma_error"
            )
            estimator.evaluate(DataSet(data["data"][-1]), output)
            print output
