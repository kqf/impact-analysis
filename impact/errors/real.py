from math import sin as Sin
from math import cos as Cos
from math import pow as Power
from math import e as E
from math import sqrt, pi

from impact.constants import k_fm, k_norm
from impact.utils import hankel_transform

from impact.parametrization.numeric import Numeric
from impact.parametrization.symbolic import Symbolic

class Error(object):
    def __init__(self, model, outname='real_gamma_error'):
        super(Error, self).__init__()
        self.outname = outname
        self.model = model

        @hankel_transform
        def evaluate_(x, dataset):
            return self.treal_error(x, dataset)

        self.evaluate_ = evaluate_


    def treal_error(self, t, dataset):
        p = dataset.parameters
        
        A = self.model.partial_derivatives(t, p)
        error_squared = 0

        # REMEMBER that a_s error should be multiplied by 1/(sqrt(pi)*4)
        for i in range(len(A)):
            for j in range(len(A)):
                error_squared += dataset.covariance[i][j] * A[i] * A[j]

        error_squared += (dataset.dsigma ** 2) * (self.model.d_as(t, p)/(sqrt(pi)*4.)) ** 2 + (dataset.drho ** 2) * (self.model.d_rho(t, p)) ** 2
        error = sqrt(error_squared)
        return error    

    def evaluate(self, dataset, output):
        output[self.outname] = map(
            lambda x: self.evaluate_(x, dataset),
            output.index
        )
        return output[self.outname]