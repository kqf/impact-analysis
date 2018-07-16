#!/usr/bin/python
import json
import unittest

from impact.datapoint import DataSet
from impact.parametrization.symbolic import FullStandard
from impact.parametrization.symbolic import FullTripleExponent
from impact.parametrization.symbolic import FullTripleExponentGeneral
from impact.vis import Plots


class RunTheSolutin(unittest.TestCase):

    # @unittest.skip("")
    def test_visualize_the_results(self):
        with open("config/input.json") as f:
            data = json.load(f)

        models = {
            # "config/full-standard.json": FullStandard,
            "config/full-triple-exponent.json": FullTripleExponent,
            # "config/full-triple-exponent-general.json": FullTripleExponentGeneral,
        }

        # TODO: Move conigs to the amplitude definitions
        #
        for config, algo in models.iteritems():
            dataset = DataSet(data["data"][-1])
            print algo.name
            visualisator = Plots()
            visualisator.draw_results(algo(), dataset, config)
