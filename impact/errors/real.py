from math import sin as Sin
from math import cos as Cos
from math import pow as Power
from math import e as E
from math import sqrt, pi

from impact.constants import k_fm, k_norm
from impact.model import hankel_transform
from partial_derivatives_explicit import PartialExplicit
from partial_derivatives_symbolic import PartialSymbolic

class Error(object):
	def __init__(self, covariance, dsigma, drho):
		super(Error, self).__init__()
		self.covariance = covariance
		self.dsigma = dsigma
		self.drho = drho
		self.partials = PartialExplicit()


		@hankel_transform
		def evaluate(x, p):
			return self.treal_error(x, p)

		self.evaluate = evaluate


	def treal_error(self, t, p):
		a1, a2, a4, b1, b2, b3, b4, a5, b5, b6, a_s, rho = p
		a_s = a_s / (sqrt(pi * k_norm) * 4)

		A = [
				self.partials.d_a1(t, p), 
				self.partials.d_a4(t, p), 
				self.partials.d_b1(t, p), 
				self.partials.d_b2(t, p), 
				0, #b4
				self.partials.d_b5(t, p), 
				# self.partials.d_as(t, p), 
				# self.partials.d_rho(t, p) 
			]

		error_squared = 0

		# REMEMBER that a_s error should be multiplied by 1/(sqrt(pi)*4)
		for i in range(len(A)):
			for j in range(len(A)):
				error_squared += self.covariance[i][j] * A[i] * A[j]

		error_squared += (self.dsigma ** 2) * (self.partials.d_as(t, p)/(sqrt(pi)*4.)) ** 2 + (self.drho ** 2) * (self.partials.d_rho(t, p)) ** 2
		error = sqrt(error_squared)
		return error	