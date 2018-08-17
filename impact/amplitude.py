from math import sqrt, pi

from impact.utils import hankel_transform
from impact.constants import k_norm


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

    def sigma(self, t, p):
        A = self.total_amplitude(t, p)
        return self.sigma_norm() * A.imag

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
        pass

    def treal_error(self, t, dataset):
        partial = self.partial_derivatives(t, dataset.parameters)
        # The covariance is modified in dataset class
        # when sigma and rho are not fitted
        covariance = dataset.covariance
        error_squared = partial.T.dot(covariance.dot(partial))
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

    def find_t_hadron(self, dataset):
        p = dataset.parameters
        for i in dataset.data:
            dsigma = self.dsigdt_norm() * abs(self.amplitude(i.t, p)) ** 2
            ratio = i.ds / dsigma
            if (ratio - 1.) < 0.01:
                return i.t
        return 9999
