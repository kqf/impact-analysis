from test.configurable import Configurable

from impact.datafit import DataFit
from impact.gammacomputor import ComputeGamma

class TestVisualRepresentation(Configurable):

	def testValues(self):
		c = ComputeGamma(self.infile, self.PROCESS, self.ENERGY, self.SIGMA, self.RHO) 
		df = DataFit(c.dataPoints, 'a', 'a', self.ENERGY, self.SIGMA, self.RHO)
		df.fit()
