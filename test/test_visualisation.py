import unittest
from impact.ComputeGamma import ComputeGamma
from impact.DataFit import DataFit

class TestVisualRepresentation(unittest.TestCase):

	def testValues(self):

		ENERGY = 7000
		RHO    = 0.14
		SIGMA  = 98.3
		DSIGMA = 2.23
		DRHO =   0.007
		PROCESS = 'pp'
		
		# TODO: Try to add quiet/dead mode
		c = ComputeGamma("input/alldata_v1_4.dat", PROCESS, ENERGY, SIGMA, RHO) 
		df = DataFit(c.dataPoints, 'a', 'a', ENERGY, SIGMA, RHO)
		df.fit()
