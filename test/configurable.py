import json
import unittest

from impact.datapoint import DataSet


class Configurable(unittest.TestCase):

    def setUp(self):
        with open('config/test.json') as f:
            self.data = json.load(f)

        self.dataset = DataSet(self.data['dataset'])
        self.longMessage = True
