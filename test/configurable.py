import unittest
import json

from impact.impactanalysis import ImpactAnalysis

class Configurable(unittest.TestCase):

	def setUp(self):
		with open('config/test.json') as f:
			self.data = json.load(f)

		p = lambda x: self.data[x]
		self.infile = p('infile')
		self.PROCESS, self.ENERGY, self.SIGMA, self.RHO, self.DSIGMA, self.DRHO = map(p, ["PROCESS", "ENERGY", "SIGMA", "RHO", "DSIGMA", "DRHO"])