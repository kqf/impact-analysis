#!/usr/bin/python
import json
import unittest

from impact.datapoint import DataSet
from impact.parametrization.symbolic import Standard
from impact.parametrization.symbolic import ThreePlusOne
from impact.parametrization.symbolic import ThreePlusTwo
from impact.vis import Plots
from analysis.utils import pack_the_dataset


class RunTheSolutin(unittest.TestCase):

    # @unittest.skip("")
    def test_visualize_the_results(self):
        with open("config/input.json") as f:
            data = json.load(f)

        pref = "config/low-energy/full-"
        models = {
            # pref + "standard.json": Standard,
            pref + "triple-exponent.json": ThreePlusOne,
            # pref + "triple-exponent-general.json": ThreePlusTwo,
        }

        # TODO: Move conigs to the amplitude definitions
        #
        datasets = [-5]
        # datasets = [0, 2, 4, -5]
        for data_index in datasets:
            dataset = DataSet(data["data"][data_index])
            for config, algo in models.iteritems():
                print algo.name
                visualisator = Plots()
                visualisator.draw_results(algo(), dataset, config)
            pack_the_dataset("result/" + str(dataset.energy) + "GeV")
