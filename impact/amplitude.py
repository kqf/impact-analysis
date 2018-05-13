from math import sqrt

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
            result = self.dsisdt_norm() * abs(A) ** 2
        except OverflowError:
            result = A.imag
        return result

    def ratio(self, t, p):
        A = self.amplitude(t, p)
        return A.real / A.imag

    def partial_derivatives(self, t, p):
        A = [
            self.d_a1(t, p),
            self.d_a2(t, p),
            self.d_b1(t, p),
            self.d_b2(t, p),
            self.d_b3(t, p),
            self.d_b4(t, p)
            # self.d_as(t, p),
            # self.d_rho(t, p)
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
