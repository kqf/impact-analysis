import random
from impact.impactanalysis import ImpactAnalysis
from test.configurable import Configurable

import pandas as pd


class TestFinalResult(Configurable):

    def setUp(self):
        super(TestFinalResult, self).setUp()
        random.seed(1234)
        self.nominal = pd.DataFrame(
            self.data['final_result']).infer_objects().set_index('b')

    def testValues(self):
        analysis = ImpactAnalysis()
        output = analysis.run(self.dataset)

        for column in self.nominal.columns:
            print column
            for pair in zip(self.nominal[column].values,
                            output[column].values):
                self.assertAlmostEqual(*pair, places=5)
