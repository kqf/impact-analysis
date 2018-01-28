from test.configurable import Configurable

from impact.datafit import DataFit
from impact.datapoint import DataPoint

class TestVisualRepresentation(Configurable):

    def testValues(self):
        data = DataPoint.read(self.ENERGY, self.PROCESS, self.infile)
        df = DataFit(data, 'a', 'a', self.ENERGY, self.SIGMA, self.RHO)
        df.fit()
        raw_input('')
