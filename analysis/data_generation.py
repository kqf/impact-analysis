#!/usr/bin/python
import json
import unittest

from impact.datapoint import DataSet
from impact.estimators import DataGenerator


class RunTheAnalysis(unittest.TestCase):

    # @unittest.skip("")
    def test_visualize_the_results(self):
        with open("config/input.json") as f:
            data = json.load(f)
        dataset = DataSet(data["data"][-1])

        generator = DataGenerator(n_iterations=10, n_sigma=2)
        generated = generator.evaluate(dataset, None)
        for dataset in generated:
            dataset.draw()
