import json
import unittest


class Configurable(unittest.TestCase):

    def setUp(self):
        with open('config/test.json') as f:
            self.data = json.load(f)

        index = self.data['test_dataset_index']
        with open('config/input.json') as f:
            data_list = json.load(f)["data"]
            params = data_list[index]

        self.sigma_rho = [98, 0.14]
        self.parms = params
        self.longMessage = True
