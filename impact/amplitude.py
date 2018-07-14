from math import sqrt, pi

from impact.utils import hankel_transform
from impact.constants import k_norm


class Amplitude(object):

    def __init__(self):
        super(Amplitude, self).__init__()

        @hankel_transform
        def imag_gamma(x, p):
            return self.total_amplitude(x, p).imag
        self.imag_gamma = imag_gamma

        @hankel_transform
        def real_gamma(x, p):
            return -self.total_amplitude(x, p).real
        self.real_gamma = real_gamma

    def diff_cs(self, t, p):
        A = self.total_amplitude(t, p)
        try:
            result = self.dsigdt_norm() * abs(A) ** 2
        except OverflowError:
            result = A.imag
        return result

    def ratio(self, t, p):
        A = self.total_amplitude(t, p)
        return A.real / A.imag

    def total_amplitude(self, t, p):
        return self.amplitude(t, p) + self.coulomb(t, p)

    def coulomb(self, t, p):
        return 0

    def partial_derivatives(self, t, p):
        A = [
            self.d_a1(t, p),
            self.d_a2(t, p),
            self.d_b1(t, p),
            self.d_b2(t, p),
            self.d_b3(t, p),
            self.d_b4(t, p)
        ]
        return A

    def treal_error(self, t, dataset):
        p = dataset.parameters

        A = self.partial_derivatives(t, p)
        error_squared = 0

        for i in range(len(A)):
            for j in range(len(A)):
                error_squared += dataset.covariance[i][j] * A[i] * A[j]

        # Sigma and rho aren't parameters of the fit
        # we can't find the covariance matrix between them
        error_squared += (
            (dataset.dsigma ** 2) * (self.d_as(t, p)) ** 2 +
            (dataset.drho ** 2) * (self.d_rho(t, p)) ** 2
        )

        error = sqrt(error_squared)
        return error

    def dsigdt_norm(self):
        return k_norm / 16. / pi

    def sigma_norm(self):
        return k_norm

    def h_norm(self):
        return 2.

    def hdata_norm(self):
        return 2. / 8 / pi
