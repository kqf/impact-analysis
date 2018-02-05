import random
from impact.model import MODEL 
from test.configurable import Configurable
import numpy as np

class TestAmplitude(Configurable):

	def setUp(self):
		super(TestAmplitude, self).setUp()
		random.seed(1234)
		self.parameters = self.data['initial_parameters'] + [self.data['SIGMA'], self.data['RHO']]
		self.real, self.imag = self.data['real_amplitude'], self.data['imag_amplitude']
		# self.nom_values, self.nom_errors = zip(*self.data[' MODEL.amplitude'])

	def npoints(self):
		return np.linspace(1e-5, 3, 100)


	def testAmplitude(self):
		data = [ MODEL.amplitude(t, self.parameters) for t in self.npoints()]
		real, imag = zip(*[[d.real, d.imag] for d in data])

		for a, b in zip(real, self.real):
				self.assertAlmostEqual(a, b)

		# Test errors
		for a, b in zip(imag, self.imag):
				self.assertAlmostEqual(a, b)
