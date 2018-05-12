#!/usr/bin/python
import json
import unittest

from impact.datapoint import DataSet
from impact.parametrization.symbolic import SymbolicUpdated
from impact.vis import Plots


class RunTheSolutin(unittest.TestCase):

    @unittest.skip('')
    def test_visualize_the_results(self):
        with open('config/input.json') as f:
            data = json.load(f)

        models = {
            # 'config/datafit.json': Symbolic,
            'config/triple-exponent.json': SymbolicUpdated
        }

        # TODO: Move conigs to the amplitude definitions
        #

        for config, algo in models.iteritems():
            dataset = DataSet(data['data'][-1])
            visualisator = Plots()
            visualisator.draw_results(algo(), dataset, config)
