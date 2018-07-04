import json
import unittest
from impact.datapoint import DataSet


class Configurable(unittest.TestCase):

    def setUp(self):
        with open('config/test.json') as f:
            self.data = json.load(f)

        index = self.data['test_dataset_index']
        with open('config/input.json') as f:
            data_list = json.load(f)["data"]
            params = data_list[index]
        self.dataset = DataSet(params)
        self.longMessage = True
