from impact.ComputeGamma import ComputeGamma
from impact.DataFit import DataFit
from test.configurable import Configurable

class TestVisualRepresentation(Configurable):

	def testValues(self):
		# TODO: Try to add quiet/dead mode
		c = ComputeGamma(self.infile, self.PROCESS, self.ENERGY, self.SIGMA, self.RHO) 
		df = DataFit(c.dataPoints, 'a', 'a', self.ENERGY, self.SIGMA, self.RHO)
		df.fit()
