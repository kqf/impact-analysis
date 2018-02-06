from test.configurable import Configurable

from impact.datafit import DataFit
from impact.datapoint import DataSet

class TestVisualRepresentation(Configurable):

    def testValues(self):
        df = DataFit('a')
        df.fit(self.dataset)
