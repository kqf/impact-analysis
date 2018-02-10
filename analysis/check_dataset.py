#!/usr/bin/python
import unittest
import json
import sys

from impact.datapoint import DataSet
from impact.impactanalysis import ImpactAnalysis

class AnalyzeSingleDataset(unittest.TestCase):
        

    def test_selected_dataset(self):
        with open('config/input.json') as f:
            data = json.load(f)

        p = data['data'][-1]
        analysis = ImpactAnalysis()
        dataset = DataSet(p)
        analysis.run(dataset)

