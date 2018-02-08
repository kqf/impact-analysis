from impact.utils import hankel_transform

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

