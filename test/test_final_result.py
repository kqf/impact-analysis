import random
import unittest

import pandas as pd
from impact.impactanalysis import ImpactAnalysis
from test.configurable import Configurable
from impact.datapoint import DataSet


class TestFinalResult(Configurable):

    def setUp(self):
        super(TestFinalResult, self).setUp()
        random.seed(1234)
        data = self.data["test_final_result"]
        self.nominal = pd.DataFrame(data).infer_objects().set_index('b')
        self.longMessage = True

    @unittest.skip("The procedure is not well established yet")
    def test_calculates_default_param(self):
        analysis = ImpactAnalysis(n_sigma=1.0)
        output = analysis.run(DataSet(self.parameters))

        message = "\n"
        for column in self.nominal.columns:
            msg = "\"{}\": {},\n"
            message += msg.format(column, list(output[column].values))

        for column in self.nominal.columns:
            for pair in zip(self.nominal[column].values,
                            output[column].values):
                self.assertAlmostEqual(*pair, places=5, msg=message)
