#!/usr/bin/python
import json
import unittest

from impact.datapoint import DataSet
from impact.parametrization.symbolic import FullStandard
from impact.parametrization.symbolic import FullTripleExponent
from impact.parametrization.symbolic import FullTripleExponentGeneral
from impact.vis import Plots
from analysis.utils import pack_the_dataset


class RunTheSolutin(unittest.TestCase):

    # @unittest.skip("")
    def test_visualize_the_results(self):
        with open("config/input.json") as f:
            data = json.load(f)

        pref = "config/low-energy/full-"
        models = {
            pref + "standard.json": FullStandard,
            pref + "triple-exponent.json": FullTripleExponent,
            pref + "triple-exponent-general.json": FullTripleExponentGeneral,
        }

        # TODO: Move conigs to the amplitude definitions
        #
        datasets = [0, 2, 4, -5]
        for data_index in datasets:
            dataset = DataSet(data["data"][data_index])
            for config, algo in models.iteritems():
                print algo.name
                visualisator = Plots()
                visualisator.draw_results(algo(), dataset, config)
            pack_the_dataset("result/" + str(dataset.energy) + "GeV")
