from test.configurable import Configurable

from impact.datafit import DataFit
from impact.datapoint import DataSet
from impact.parametrization.numeric import Numeric

class TestVisualRepresentation(Configurable):

    def testValues(self):
        df = DataFit('a', Numeric())
        df.fit(self.dataset)
