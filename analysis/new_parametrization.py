#!/usr/bin/python
import unittest
import json
import sys


from impact.parametrization.symbolic import Symbolic, SymbolicUpdated
from impact.datapoint import DataSet
from impact.vis import Plots
from impact.impactanalysis import *




class ValidateNewSolution(unittest.TestCase):
        
    def test_parametrization_fits_the_data(self):
        with open('config/input.json') as f:
            data = json.load(f)

        model = SymbolicUpdated()

        totem7tev = DataSet(data['data'][-2])

        visualisator = Plots()
        visualisator.fit(model, totem7tev, 'config/triple-exponent.json')   