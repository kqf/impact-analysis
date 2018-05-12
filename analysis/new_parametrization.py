#!/usr/bin/python
import json
import unittest

from impact.datapoint import DataSet
from impact.parametrization.symbolic import Symbolic
# from impact.parametrization.symbolic import TripleExponent
from impact.vis import Plots


class ValidateNewSolution(unittest.TestCase):

    def test_parametrization_fits_the_data(self):
        with open('config/input.json') as f:
            data = json.load(f)

        model = Symbolic()

        totem7tev = DataSet(data['data'][-2])

        visualisator = Plots()
        # visualisator.draw_results(model,
        # totem7tev, 'config/triple-exponent.json')
        visualisator.draw_results(model, totem7tev, 'config/datafit.json')

    @unittest.skip('')
    def test_the_results(self):
        with open('config/input.json') as f:
            data = json.load(f)

        model = Symbolic()

        totem7tev = DataSet(data['data'][0])

        visualisator = Plots()
        # visualisator.draw_results(model,
        # totem7tev, 'config/triple-exponent.json')
        visualisator.draw_results(model, totem7tev, 'config/datafit.json')
