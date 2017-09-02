from test.configurable import Configurable

from impact.datafit import DataFit
from impact.datapoint import DataReader

class TestVisualRepresentation(Configurable):

    def testValues(self):
        data = DataReader(self.ENERGY, self.PROCESS).read(self.infile)
        df = DataFit(data, 'a', 'a', self.ENERGY, self.SIGMA, self.RHO)
        df.fit()
