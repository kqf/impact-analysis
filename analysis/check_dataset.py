#!/usr/bin/python
import unittest
import json
import sys

from impact.impactanalysis import ImpactAnalysis

class AnalyzeSingleDataset(unittest.TestCase):
        

    def test_selected_dataset(self):
        with open('config/input.json') as f:
            data = json.load(f)

        infile = data['infile']
        
        p = data['data'][-1]
        analysis = ImpactAnalysis(infile, p, mode='v')
        values, errors = analysis.run()
