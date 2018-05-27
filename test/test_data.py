import json
from test.configurable import Configurable
from impact.datapoint import DataSet


class TestDataset(Configurable):

    def test_input_datasets(self):
        with open('config/input.json') as f:
            input_data = json.load(f)

        output = [
            len(DataSet(dset).data)
            for dset in input_data['data']
        ]
        self.assertEqual(output, self.data["dataset_sizes"])
