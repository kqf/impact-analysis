import json

import pytest
from test.configurable import Configurable
from impact.datapoint import DataSet


class TestDataset(Configurable):

    @pytest.mark.onlylocal
    def test_input_datasets(self):
        with open('config/input.json') as f:
            input_data = json.load(f)

        output = [
            len(DataSet(dset).data)
            for dset in input_data['data']
        ]
        data = self.data["test_data"]
        self.assertEqual(output, data["dataset_sizes"])
