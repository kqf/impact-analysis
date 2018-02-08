from scipy import integrate
from scipy.special import j0, j1
from constants import k_fm, k_norm
from math import sqrt, pi


def hankel_transform(func):
    def impact_version(b, p, limits = (0, float("inf"))):
        f = lambda q : q * j0(b * q / k_fm) *  func(q * q, p) / sqrt(pi * k_norm)
        result = integrate.quad(f, *limits)[0]  # integral from zero to lower bound
        return result

    return impact_version

class Amplitude(object):

	def __init__(self):
		super(Amplitude, self).__init__()

		@hankel_transform
		def imag_gamma(x, p):
			return self.amplitude(x, p).imag
		self.imag_gamma = imag_gamma

		@hankel_transform
		def real_gamma(x, p):
			return -self.amplitude(x, p).real	
		self.real_gamma = real_gamma

	def diff_cs(self, t, p):
		A = self.amplitude(t, p)
		try:
			result = abs(A) ** 2
		except OverflowError:
			result = A.imag
		return result

	def ratio(self, t, p):
	    A = self.amplitude(t,p)
	    return A.real / A.imag

