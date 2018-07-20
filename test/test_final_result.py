import random

import pandas as pd
from impact.impactanalysis import ImpactAnalysis
from test.configurable import Configurable


class TestFinalResult(Configurable):

    def setUp(self):
        super(TestFinalResult, self).setUp()
        random.seed(1234)
        data = self.data["test_final_result"]
        self.nominal = pd.DataFrame(data).infer_objects().set_index('b')
        self.longMessage = True

    def test_calculates_default_param(self):
        analysis = ImpactAnalysis(n_sigma=1.0)
        output = analysis.run(self.dataset)

        message = "\n"
        for column in self.nominal.columns:
            msg = "\"{}\": {},\n"
            message += msg.format(column, list(output[column].values))

        for column in self.nominal.columns:
            for pair in zip(self.nominal[column].values,
                            output[column].values):
                self.assertAlmostEqual(*pair, places=5, msg=message)
