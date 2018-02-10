#!/usr/bin/python
import unittest
import json
import sys

from impact.datapoint import DataSet
from impact.impactanalysis import ImpactAnalysis


from test.configurable import Configurable
from impact.parametrization.symbolic import Symbolic
from impact.vis import Plots




class ValidateNewSolution(unittest.TestCase):
        
    def test_parametrization_fits_the_data(self):
        with open('config/input.json') as f:
            data = json.load(f)

        totem7tev = DataSet(data['data'][-2])

        visualisator = Plots()
        visualisator.fit(Symbolic(), totem7tev)   