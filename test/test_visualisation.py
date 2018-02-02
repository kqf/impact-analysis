from test.configurable import Configurable

from impact.datafit import DataFit
from impact.datapoint import DataSet

class TestVisualRepresentation(Configurable):

    def testValues(self):
        data = DataSet(self.infile, self.data).data
	    
        df = DataFit('a', 'a', self.ENERGY, self.SIGMA, self.RHO)
        df.fit(data)
        raw_input('')
