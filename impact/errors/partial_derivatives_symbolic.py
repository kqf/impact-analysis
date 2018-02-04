from math import sin as Sin
from math import cos as Cos
from math import pow as Power
from math import e as E
from math import sqrt, pi
import sympy as smp

from impact.constants import k_norm

class PartialSymbolic(object):

	def __init__(self):
		super(PartialSymbolic, self).__init__()
		self.variable_names = 'a1 a2 a4 b1 b2 b3 b4 a5 b5 b6 a_s rho'.split()
		self.variables = smp.symbols(self.variable_names)
		self.t = smp.Symbol('t')
		self.amplitude = self.analytic_formula()

		useful_pars = 'a1 a4 b1 b2 b4 b5 a_s rho'.split()
		self.d_a1, self.d_a4, self.d_b1, self.d_b2, self.d_b4, self.d_b5, self.d_as, self.d_rho = \
			map(self._partial, useful_pars)


	def analytic_formula(self):
		a1, a2, a4, b1, b2, b3, b4, a5, b5, b6, a_s, rho =  self.variables
		t = self.t
		a_s = a_s / (smp.sqrt(smp.pi *  k_norm) * 4)
		alpha = (1 - 1j*rho)*(a_s + a4)
		return 1j*alpha*( a1*smp.exp(-0.5*alpha*b1*t) + (1 - a1)*smp.exp(-0.5*alpha*b2*t) ) - 1j*a4*smp.exp(-0.5*b4*t) - a4*rho/((1 + t/b5)**4)


	def _partial(self, arg):
		fpar = next((i for i in self.amplitude.free_symbols if i.name == arg), None)

		partial_derivative = smp.lambdify(
			(self.t, self.variables), 
			smp.re(self.amplitude.diff(fpar)), 
			'numpy'
		)

		return partial_derivative

 	 