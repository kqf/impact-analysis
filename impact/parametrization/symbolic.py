from math import pi, sqrt

import sympy as smp
from cmath import exp as Exp
from impact.amplitude import Amplitude
from impact.constants import k_norm


class Symbolic(Amplitude):
    name = "standard-parametrization"

    def __init__(self):
        super(Symbolic, self).__init__()
        self.name = "standard-parametrization"
        variable_names = 'a1 a2 b1 b2 b3 b4 a_s rho'.split()

        self.variables = smp.symbols(variable_names)
        self.t = smp.Symbol('t')
        self.symbol_amplitude = self.analytic_formula()

        self.d_a1, self.d_a2, self.d_b1, self.d_b2,\
            self.d_b3, self.d_b4, self.d_as, self.d_rho = \
            map(self._partial, variable_names)

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
        # a_s = a_s / (smp.sqrt(smp.pi * k_norm) * 4)
        alpha = (1 - 1j * rho) * (a_s + a2)
        return (
            1j * alpha * (
                smp.exp(-0.5 * alpha * b1 * t) * a1 +
                smp.exp(-0.5 * alpha * b2 * t) * (1 - a1)
            ) -
            1j * smp.exp(-0.5 * b3 * t) * a2 -
            a2 * rho / ((1 + t / b4)**4)
        )

    def _partial(self, arg):
        symbols = self.symbol_amplitude.free_symbols
        fpar = next(
            (i for i in symbols if i.name == arg), None)

        partial_derivative = smp.lambdify(
            (self.t, self.variables),
            smp.re(self.symbol_amplitude.diff(fpar)),
            'numpy')

        return partial_derivative

    def amplitude(self, t, p):
        a1, a2, b1, b2, b3, b4, a_s, rho = p

        a_s = a_s / (sqrt(pi * k_norm) * 4)

        alpha = (1 - 1j * rho) * (a_s + a2)

        try:
            ampl = (
                1j * alpha * (
                    Exp(-0.5 * alpha * b1 * t) * a1 +
                    Exp(-0.5 * alpha * b2 * t) * (1 - a1)
                ) -
                1j * Exp(-0.5 * b3 * t) * a2 - a2 * rho / ((1 + t / b4)**4))
        except OverflowError:
            ampl = 0
        return ampl


class TripleExponent(Symbolic):

    def __init__(self):
        super(TripleExponent, self).__init__()
        self.name = "three-exponents"

    def analytic_formula(self):
        a1, a2, b1, b2, b3, b4, a_s, rho = self.variables
        # a_s = a_s / (smp.sqrt(smp.pi * k_norm) * 4)

        a3 = -a1 - a2 + 0.5 * a_s / k_norm
        a4 = rho * (a1 + a2 + a3)

        amplitude = (
            a1 * smp.exp(b1 * self.t) * 1j +
            a2 * smp.exp(b2 * self.t) * 1j +
            a3 * smp.exp(b3 * self.t) * 1j +
            a4 * smp.exp(b4 * self.t)
        )
        return amplitude

    def amplitude(self, t, p):
        a1, a2, b1, b2, b3, b4, a_s, rho = p
        # a_s = a_s / (sqrt(pi * k_norm) * 4)

        a3 = -a1 - a2 + 0.5 * a_s / k_norm
        a4 = rho * (a1 + a2 + a3)

        from math import exp

        def amplitude(tt):
            return (
                a1 * exp(b1 * tt) * 1j +
                a2 * exp(b2 * tt) * 1j +
                a3 * exp(b3 * tt) * 1j +
                a4 * exp(b4 * tt)
            )

        try:
            ampl = amplitude(t)
        except OverflowError:
            ampl = 0

        # print '>>> ', abs(ampl) ** 2
        return ampl

    def dsigdt_norm(self):
        return k_norm / 4. / pi

    def sigma_norm(self):
        return 2 * k_norm

    def h_norm(self):
        return 4.

    def hdata_norm(self):
        return 4. / 8 / pi


class TripleExponentGeneral(Amplitude):
    name = "three-exponents-extended"

    def __init__(self):
        super(TripleExponentGeneral, self).__init__()
        variable_names = 'a1 a2 a5 b1 b2 b3 b4 b5 a_s rho'.split()
        self.name = "three-exponents-extended"

        self.variables = smp.symbols(variable_names)
        self.t = smp.Symbol('t')
        self.symbol_amplitude = self.analytic_formula()

        self.d_a1, self.d_a2, self.d_a5, self.d_b1, self.d_b2,\
            self.d_b3, self.d_b4, self.d_b5, self.d_as, self.d_rho = \
            map(self._partial, variable_names)

        # Calculate using sympy
        #
        self.amplitude_symbolic = smp.lambdify(
            (self.t, self.variables),
            self.symbol_amplitude,
            'numpy'
        )

    def analytic_formula(self):
        print len(self.variables)
        a1, a2, a5, b1, b2, b3, b4, b5, a_s, rho = self.variables
        # a_s = a_s / (smp.sqrt(smp.pi * k_norm) * 4)

        a3 = -a1 - a2 + 0.5 * a_s / k_norm
        a4 = rho * (a1 + a2 + a3)

        amplitude = (
            a1 * smp.exp(b1 * self.t) * 1j +
            a2 * smp.exp(b2 * self.t) * 1j +
            a3 * smp.exp(b3 * self.t) * 1j +
            a4 * smp.exp(b4 * self.t) +
            a5 * smp.exp(b5 * self.t)
        )
        return amplitude

    def amplitude(self, t, p):
        a1, a2, a5, b1, b2, b3, b4, b5, a_s, rho = p
        # TODO: Fix me
        # a_s = a_s / (sqrt(pi * k_norm) * 4)

        a3 = -a1 - a2 + 0.5 * a_s / k_norm
        a4 = rho * (a1 + a2 + a3)

        from math import exp

        def amplitude(tt):
            return (
                a1 * exp(b1 * tt) * 1j +
                a2 * exp(b2 * tt) * 1j +
                a3 * exp(b3 * tt) * 1j +
                a4 * exp(b4 * tt) +
                a5 * exp(b5 * tt)
            )

        try:
            ampl = amplitude(t)
        except OverflowError:
            ampl = 0

        # print '>>> ', abs(ampl) ** 2
        return ampl

    def _partial(self, arg):
        symbols = self.symbol_amplitude.free_symbols
        fpar = next(
            (i for i in symbols if i.name == arg), None)

        partial_derivative = smp.lambdify(
            (self.t, self.variables),
            smp.re(self.symbol_amplitude.diff(fpar)),
            'numpy')

        return partial_derivative

    def partial_derivatives(self, t, p):
        A = [
            self.d_a1(t, p),
            self.d_a2(t, p),
            self.d_a5(t, p),
            self.d_b1(t, p),
            self.d_b2(t, p),
            self.d_b3(t, p),
            self.d_b4(t, p),
            self.d_b5(t, p),
        ]
        return A

    def dsigdt_norm(self):
        return k_norm / 4. / pi

    def sigma_norm(self):
        return 2 * k_norm

    def h_norm(self):
        return 4.

    def hdata_norm(self):
        return 4. / 8 / pi
