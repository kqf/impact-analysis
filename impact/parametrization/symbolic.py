from math import pi, log

import numpy as np
import sympy as smp
from cmath import exp
from impact.amplitude import Amplitude
from impact.constants import k_norm


class Coulomb(object):
    def coulomb(self, t, p):
        if not t > 0:
            return 0
        amcp = 0.71
        eg = 0.577
        aem = 0.007297
        ffp = 1. / (1. + t / amcp) ** 2
        phase = aem * (log(10.2 * t) + eg)
        coulomb = - 8 * pi * aem * ffp ** 2 * exp(-1j * phase) / t
        return coulomb


class Standard(Amplitude):
    name = "standard-parametrization"
    variable_names = 'a1 a2 b1 b2 b3 b4 a_s rho'.split()

    def __init__(self):
        super(Standard, self).__init__()
        self._partials = None
        self.variables = smp.symbols(self.variable_names)
        self.t = smp.Symbol('t')
        self.symbol_amplitude = self.analytic_formula()
        self.amplitude_symbolic = smp.lambdify(
            (self.t, self.variables),
            self.symbol_amplitude,
            'numpy'
        )

    @property
    def partials(self):
        if not self._partials:
            self._partials = map(self._partial, self.variables)
        return self._partials

    def analytic_formula(self):
        a1, a2, b1, b2, b3, b4, a_s, rho = self.variables
        t = self.t
        a_s = a_s / k_norm
        alpha = (1 - 1j * rho) * (a_s + a2)
        return (
            1j * alpha * (
                smp.exp(-0.5 * alpha * b1 * t) * a1 +
                smp.exp(-0.5 * alpha * b2 * t) * (1 - a1)
            ) -
            1j * smp.exp(-0.5 * b3 * t) * a2 -
            a2 * rho / ((1 + t / b4)**4)
        )

    def _partial(self, fpar):
        partial_derivative = smp.lambdify(
            (self.t, self.variables),
            smp.re(self.symbol_amplitude.diff(fpar)),
            'numpy')
        return partial_derivative

    def partial_derivatives(self, t, p):
        return np.asarray([f(t, p) for f in self.partials])

    def amplitude(self, t, p):
        a1, a2, b1, b2, b3, b4, a_s, rho = p
        a_s = a_s / k_norm

        alpha = (1 - 1j * rho) * (a_s + a2)

        try:
            ampl = (
                1j * alpha * (
                    exp(-0.5 * alpha * b1 * t) * a1 +
                    exp(-0.5 * alpha * b2 * t) * (1 - a1)
                ) -
                1j * exp(-0.5 * b3 * t) * a2 - a2 * rho / ((1 + t / b4)**4))
        except OverflowError:
            ampl = 0
        return ampl


class TripleExponent(Standard):
    name = "three-exponents"

    def analytic_formula(self):
        a1, a2, b1, b2, b3, b4, a_s, rho = self.variables

        a3 = -a1 - a2 + a_s / k_norm
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

        a3 = -a1 - a2 + a_s / k_norm
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


class TripleExponentGeneral(Amplitude):
    name = "three-exponents-extended"
    variable_names = 'a1 a2 a6 b1 b2 b3 b4 b5 a_s rho'.split()

    def analytic_formula(self):
        print len(self.variables)
        a1, a2, a6, b1, b2, b3, b4, b5, a_s, rho = self.variables

        a3 = -a1 - a2 + a_s / k_norm
        a4 = rho * (a1 + a2 + a3)

        amplitude = (
            a1 * smp.exp(b1 * self.t) * 1j +
            a2 * smp.exp(b2 * self.t) * 1j +
            a3 * smp.exp(b3 * self.t) * 1j +
            a4 * smp.exp(b4 * self.t) +
            a6 * smp.exp(b5 * self.t)
        )
        return amplitude

    def amplitude(self, t, p):
        a1, a2, a6, b1, b2, b3, b4, b5, a_s, rho = p
        # TODO: Fix me

        a3 = -a1 - a2 + a_s / k_norm
        a4 = rho * (a1 + a2 + a3) - a6

        from math import exp

        def amplitude(tt):
            return (
                a1 * exp(b1 * tt) * 1j +
                a2 * exp(b2 * tt) * 1j +
                a3 * exp(b3 * tt) * 1j +
                a4 * exp(b4 * tt) +
                a6 * exp(b5 * tt)
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


class FullStandard(Coulomb, Standard):
    pass


class ThreePlusOne(Standard):
    name = "three-plus-one"
    variable_names = 'a1 a2 a6 b1 b2 b3 b4 b6 a_s rho'.split()

    def analytic_formula(self):
        a1, a2, a6, b1, b2, b3, b4, b6, a_s, rho = self.variables
        a3 = -a1 - a2 - a6 + a_s / k_norm
        a4 = rho * (a1 + a2 + a3 + a6)
        t = -self.t
        amplitude = (
            a1 * smp.exp(b1 * t) * 1j +
            a2 * smp.exp(b2 * t) * 1j +
            a3 * smp.exp(b3 * t) * 1j +
            a4 * smp.exp(b4 * t) +
            a6 * 1j / (1 - t / b6) ** 4)
        return amplitude

    def amplitude(self, t, p):
        a1, a2, a6, b1, b2, b3, b4, b6, a_s, rho = p

        a3 = -a1 - a2 - a6 + a_s / k_norm
        a4 = rho * (a1 + a2 + a3 + a6)

        from math import exp

        def amplitude(tt):
            return (
                a1 * exp(b1 * tt) * 1j +
                a2 * exp(b2 * tt) * 1j +
                a3 * exp(b3 * tt) * 1j +
                a4 * exp(b4 * tt) +
                a6 * 1j / (1 - tt / b6) ** 4
            )

        try:
            ampl = amplitude(-t)
        except OverflowError:
            ampl = 0

        return ampl


class FullThreePlusOne(Coulomb, ThreePlusOne):
    name = "full-three-plus-one"
    pass


class ThreePlusTwo(Standard):
    name = "full-three-plus-two"
    variable_names = 'a1 a2 a5 b1 b2 b3 b4 b5 a_s rho'.split()

    def analytic_formula(self):
        a1, a2, a5, b1, b2, b3, b4, b5, a_s, rho = self.variables
        a3 = -a1 - a2 + a_s / k_norm
        a4 = rho * (a1 + a2 + a3)
        t = -self.t
        amplitude = (
            a1 * smp.exp(b1 * t) * 1j +
            a2 * smp.exp(b2 * t) * 1j +
            a3 * smp.exp(b3 * t) * 1j +
            a4 * smp.exp(b4 * t) +
            a5 * smp.exp(b5 * t)
        )
        return amplitude

    def amplitude(self, t, p):
        a1, a2, a5, b1, b2, b3, b4, b5, a_s, rho = p

        a3 = -a1 - a2 + a_s / k_norm
        a4 = rho * (a1 + a2 + a3) - a5

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
            ampl = amplitude(-t)
        except OverflowError:
            ampl = 0

        return ampl


class FullThreePlusTwo(Coulomb, ThreePlusOne):
    name = "full-three-plus-two"
    pass
