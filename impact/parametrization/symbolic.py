from math import sin as Sin
from math import cos as Cos
from math import pow as Power
from math import e as E
from math import sqrt, pi
import sympy as smp
from cmath import exp as Exp

from impact.constants import k_norm
from impact.amplitude import Amplitude

class Symbolic(Amplitude):

    def __init__(self):
        super(Symbolic, self).__init__()
        self.variable_names = 'a1 a2 b1 b2 b3 b4 a_s rho'.split()
        self.variables = smp.symbols(self.variable_names)
        self.t = smp.Symbol('t')
        self.symbol_amplitude = self.analytic_formula()

        useful_pars = 'a1 a2 b1 b2 b3 b4 a_s rho'.split()
        self.d_a1, self.d_a2, self.d_b1, self.d_b2, self.d_b3, self.d_b4, self.d_as, self.d_rho = \
            map(self._partial, useful_pars)

        # Calculate using sympy
        #
        self.amplitude_symbolic = smp.lambdify(
            (self.t, self.variables),
            self.symbol_amplitude,
            'numpy'
        )


    def analytic_formula(self):
        a1, a2, b1, b2, b3, b4, a_s, rho = self.variables
        t = self.t
        a_s = a_s / (smp.sqrt(smp.pi *  k_norm) * 4)
        alpha = (1 - 1j*rho)*(a_s + a2)
        return 1j*alpha*( a1*smp.exp(-0.5*alpha*b1*t) + (1 - a1)*smp.exp(-0.5*alpha*b2*t) ) - 1j*a2*smp.exp(-0.5*b3*t) - a2*rho/((1 + t/b4)**4)


    def _partial(self, arg):
        fpar = next((i for i in self.symbol_amplitude.free_symbols if i.name == arg), None)


        partial_derivative = smp.lambdify(
            (self.t, self.variables), 
            smp.re(self.symbol_amplitude.diff(fpar)), 
            'numpy'
        )

        return partial_derivative       

        
    def amplitude(self, t, p):
        a1, a2, b1, b2, b3, b4, a_s, rho = p

        a_s = a_s/(sqrt(pi *  k_norm) * 4)
        
        alpha = (1 - 1j*rho)*(a_s + a2)

        try:
            ampl= 1j*alpha*(a1 * Exp(-0.5*alpha*b1*t) + (1 - a1)*Exp(-0.5*alpha*b2*t) ) - 1j*a2*Exp(-0.5*b3*t) - a2*rho/((1 + t/b4)**4)
        except OverflowError:
            ampl = 0
        return ampl 

        
class SymbolicUpdated(Symbolic):

    def __init__(self):
        super(SymbolicUpdated, self).__init__()

    def analytic_formula(self):
        a1, a2, b1, b2, b3, b4, a_s, rho =  self.variables
        a_s = a_s / (smp.sqrt(smp.pi *  k_norm) * 4)

        a3 = -a1 - a2 + 0.5 * a_s / k_norm
        a4 = rho * (a1 + a2 + a3)

        amplitude = (
              a1 * smp.exp(b1 * self.t) * 1j
            + a2 * smp.exp(b2 * self.t) * 1j
            + a3 * smp.exp(b3 * self.t) * 1j
            + a4 * smp.exp(b4 * self.t)
        )
        return amplitude

    def amplitude(self, t, p):
        a1, a2, b1, b2, b3, b4, a_s, rho = p
        # a_s = a_s / (sqrt(pi * k_norm) * 4)
        
        a3 = -a1 - a2 + 0.5 * a_s / k_norm
        a4 = rho * (a1 + a2 + a3)

        from math import exp
        amplitude = lambda tt: (
              a1 * exp(b1 * tt) * 1j
            + a2 * exp(b2 * tt) * 1j
            + a3 * exp(b3 * tt) * 1j
            + a4 * exp(b4 * tt)
        )

        try:
            ampl = amplitude(t) * sqrt(k_norm / 4. / pi)
        except OverflowError:
            ampl = 0

        # print '>>> ', abs(ampl) ** 2
        return ampl

    def diff_cs(self, t, p):
        cs = super(SymbolicUpdated, self).diff_cs(t, p)
        return cs
