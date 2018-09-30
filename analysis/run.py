#!/usr/bin/python
import json
import unittest

from impact.datapoint import DataSet
from impact.parametrization.symbolic import FullStandard
from impact.parametrization.symbolic import FullThreePlusOne
from impact.parametrization.symbolic import FullThreePlusTwo
from impact.vis import Plots
from analysis.utils import pack_the_dataset


class RunTheSolutin(unittest.TestCase):

    # @unittest.skip("")
    def test_visualize_the_results(self):
        with open("config/input.json") as f:
            data = json.load(f)

        models = {
            # "config/full-standard.json": FullStandard,
            "config/full-triple-exponent.json": FullThreePlusOne,
            # "config/full-triple-exponent-general.json": FullThreePlusTwo,
        }

        # TODO: Move conigs to the amplitude definitions
        #
        for config, algo in models.iteritems():
            dataset = DataSet(data["data"][-1])
            print algo.name
            visualisator = Plots()
            visualisator.draw_results(algo(), dataset, config)
        pack_the_dataset(str(dataset.energy) + "GeV")
