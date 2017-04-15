#!/usr/bin/python2.7

from cmath import exp as Exp
from math  import exp, pi, fabs, sqrt, log, sin, cos
from scipy import integrate
from scipy.special import j0, j1
from random import gauss

from random import randrange

from math import sin as Sin
from math import cos as Cos
from math import pow as Power
from math import e as E

# Normalization constants
k_fm   = 0.1973269718
k_norm = 0.389379338

import numpy as np

def amplitude(t, p):
    a1, a2, a4, b1, b2, b3, b4, a5, b5, b6, a_s, rho = p
    a_s = a_s/(sqrt(pi *  k_norm) * 4)
    
    alpha = (1 - 1j*rho)*(a_s + a4)

    try:
        ampl= 1j*alpha*( a1*Exp(-0.5*alpha*b1*t) + (1 - a1)*Exp(-0.5*alpha*b2*t) ) - 1j*a4*Exp(-0.5*b4*t) - a4*rho/((1 + t/b5)**4)
    except OverflowError:
        ampl = 0

    return ampl

def getRealError(t, p, covariance, dsigma, drho):
    a1, a2, a4, b1, b2, b3, b4, a5, b5, b6, a_s, rho = p
    a_s = a_s / (sqrt(pi) * 4)
    
 
    d_a1 = (a4 + a_s)*(
                   rho*( Power(E,((-a4 - a_s)*b1*t)/2.)*Cos(((a4 + a_s)*b1*t*rho)/2.)
                       - Power(E,((-a4 - a_s)*b2*t)/2.)*Cos(((a4 + a_s)*b2*t*rho)/2.) )
                       - Power(E,((-a4 - a_s)*b1*t)/2.)*Sin(((a4 + a_s)*b1*t*rho)/2.) 
                       + Power(E,((-a4 - a_s)*b2*t)/2.)*Sin(((a4 + a_s)*b2*t*rho)/2.)
                       )
    d_a4 = ( 
                - (rho/Power(1 + t/b5,4)) 
                + rho*(   a1*Power(E,((-a4 - a_s)*b1*t)/2.)*Cos(((a4 + a_s)*b1*t*rho)/2.) 
                + (1 - a1)*Power(E,((-a4 - a_s)*b2*t)/2.)*Cos(((a4 + a_s)*b2*t*rho)/2.) ) 
                - a1*Power(E,((-a4 - a_s)*b1*t)/2.)*Sin(((a4 + a_s)*b1*t*rho)/2.) 
                - (1 - a1)*Power(E,((-a4 - a_s)*b2*t)/2.)*Sin(((a4 + a_s)*b2*t*rho)/2.) 

                + (a4 + a_s)
                *(
                    -(a1*b1*Power(E,((-a4 - a_s)*b1*t)/2.)*t*rho*Cos(((a4 + a_s)*b1*t*rho)/2.))/2. 
                    - ((1 - a1)*b2*Power(E,((-a4 - a_s)*b2*t)/2.)*t*rho*Cos(((a4 + a_s)*b2*t*rho)/2.))/2. 
                    + (a1*b1*Power(E,((-a4 - a_s)*b1*t)/2.)*t*Sin(((a4 + a_s)*b1*t*rho)/2.))/2. 
                    + ((1 - a1)*b2*Power(E,((-a4 - a_s)*b2*t)/2.)*t*Sin(((a4 + a_s)*b2*t*rho)/2.))/2. 
                    + rho*( - (a1*b1*Power(E,((-a4 - a_s)*b1*t)/2.)*t*Cos(((a4 + a_s)*b1*t*rho)/2.))/2. 
                            - ((1 - a1)*b2*Power(E,((-a4 - a_s)*b2*t)/2.)*t*Cos(((a4 + a_s)*b2*t*rho)/2.))/2. 
                            - (a1*b1*Power(E,((-a4 - a_s)*b1*t)/2.)*t*rho*Sin(((a4 + a_s)*b1*t*rho)/2.))/2. 
                            - ((1 - a1)*b2*Power(E,((-a4 - a_s)*b2*t)/2.)*t*rho*Sin(((a4 + a_s)*b2*t*rho)/2.))/2. 
                          )
                  )
            )

    d_as = ( 
             rho*(a1*Power(E,((-a4 - a_s)*b1*t)/2.)*Cos(((a4 + a_s)*b1*t*rho)/2.)
           + (1 - a1)*Power(E,((-a4 - a_s)*b2*t)/2.)*Cos(((a4 + a_s)*b2*t*rho)/2.))
           - a1*Power(E,((-a4 - a_s)*b1*t)/2.)*Sin(((a4 + a_s)*b1*t*rho)/2.)
           - (1 - a1)*Power(E,((-a4 - a_s)*b2*t)/2.)*Sin(((a4 + a_s)*b2*t*rho)/2.) 
           + (a4 + a_s)*(-(a1*b1*Power(E,((-a4 - a_s)*b1*t)/2.)*t*rho*Cos(((a4 + a_s)*b1*t*rho)/2.))/2.
           - ((1 - a1)*b2*Power(E,((-a4 - a_s)*b2*t)/2.)*t*rho*Cos(((a4 + a_s)*b2*t*rho)/2.))/2.
           + (a1*b1*Power(E,((-a4 - a_s)*b1*t)/2.)*t*Sin(((a4 + a_s)*b1*t*rho)/2.))/2.
           + ((1 - a1)*b2*Power(E,((-a4 - a_s)*b2*t)/2.)*t*Sin(((a4 + a_s)*b2*t*rho)/2.))/2.
           + rho*(-(a1*b1*Power(E,((-a4 - a_s)*b1*t)/2.)*t*Cos(((a4 + a_s)*b1*t*rho)/2.))/2.
           - ((1 - a1)*b2*Power(E,((-a4 - a_s)*b2*t)/2.)*t*Cos(((a4 + a_s)*b2*t*rho)/2.))/2.
           - (a1*b1*Power(E,((-a4 - a_s)*b1*t)/2.)*t*rho*Sin(((a4 + a_s)*b1*t*rho)/2.))/2.
           - ((1 - a1)*b2*Power(E,((-a4 - a_s)*b2*t)/2.)*t*rho*Sin(((a4 + a_s)*b2*t*rho)/2.))/2.))
           )
    d_b1 = (a4 + a_s)*(
                         - (a1*(a4 + a_s)*Power(E,((-a4 - a_s)*b1*t)/2.)*t*rho*Cos(((a4 + a_s)*b1*t*rho)/2.))/2.
                         - (a1*(-a4 - a_s)*Power(E,((-a4 - a_s)*b1*t)/2.)*t*Sin(((a4 + a_s)*b1*t*rho)/2.))/2. 
                         + rho*((a1*(-a4 - a_s)*Power(E,((-a4 - a_s)*b1*t)/2.)*t*Cos(((a4 + a_s)*b1*t*rho)/2.))/2.
                         - (a1*(a4 + a_s)*Power(E,((-a4 - a_s)*b1*t)/2.)*t*rho*Sin(((a4 + a_s)*b1*t*rho)/2.))/2.)
                       )
 

 
    d_b2 = (a4 + a_s)*(
                        - ((1 - a1)*(a4 + a_s)*Power(E,((-a4 - a_s)*b2*t)/2.)*t*rho*Cos(((a4 + a_s)*b2*t*rho)/2.))/2. 
                        - ((1 - a1)*(-a4 - a_s)*Power(E,((-a4 - a_s)*b2*t)/2.)*t*Sin(((a4 + a_s)*b2*t*rho)/2.))/2. 
                        + rho*( ((1 - a1)*(-a4 - a_s)*Power(E,((-a4 - a_s)*b2*t)/2.)*t*Cos(((a4 + a_s)*b2*t*rho)/2.))/2. 
                              - ((1 - a1)*(a4 + a_s)*Power(E,((-a4 - a_s)*b2*t)/2.)*t*rho*Sin(((a4 + a_s)*b2*t*rho)/2.))/2.)
                       )
    d_b5 = (-4*a4*t*rho)/(Power(b5,2)*Power(1 + t/b5,5))

    d_rho = ( -(a4/Power(1 + t/b5,4)) + 
               (a4 + a_s)*(
                            a1*Power(E,((-a4 - a_s)*b1*t)/2.)*Cos(((a4 + a_s)*b1*t*rho)/2.)
                         - (a1*(a4 + a_s)*b1*Power(E,((-a4 - a_s)*b1*t)/2.)*t*Cos(((a4 + a_s)*b1*t*rho)/2.))/2.
                         + (1 - a1)*Power(E,((-a4 - a_s)*b2*t)/2.)*Cos(((a4 + a_s)*b2*t*rho)/2.)
                         - ((1 - a1)*(a4 + a_s)*b2*Power(E,((-a4 - a_s)*b2*t)/2.)*t*Cos(((a4 + a_s)*b2*t*rho)/2.))/2.
                         + rho*(-(a1*(a4 + a_s)*b1*Power(E,((-a4 - a_s)*b1*t)/2.)*t*Sin(((a4 + a_s)*b1*t*rho)/2.))/2.
                         - ((1 - a1)*(a4 + a_s)*b2*Power(E,((-a4 - a_s)*b2*t)/2.)*t*Sin(((a4 + a_s)*b2*t*rho)/2.))/2.)
                         )
            )

    A = [ 
          (d_a1),
          (d_a4),
          (d_b1),
          (d_b2),
          0, #b4
          (d_b5)
           # (d_as),
           # (d_rho)
        ]

    error_squared = 0

    # REMEMBER that a_s error should be multiplied by 1/(sqrt(pi)*4)


    for i in range(len(A)):
        for j in range(len(A)):
            error_squared += covariance[i][j] * A[i] * A[j]

    error_squared += (dsigma ** 2) * (d_as/(sqrt(pi)*4.))**2 + (drho ** 2) * (d_rho)**2
    error = sqrt(error_squared)
    return error


def getRealGammaError(b, p, covariance, dsigma, drho):
    f = lambda q :  q*j0(b*q/k_fm)*getRealError(q*q, p, covariance, dsigma, drho)/sqrt(pi*k_norm)
    result =  integrate.quad(f, 0, np.infty)[0]  # integral from zero to lower bound

    # imGamma = - \int reAmplitude
    return result


def hankel_transform(func):
    def impact_version(b, p, limits = (0, np.infty)):
        f = lambda q : q * j0(b * q / k_fm) *  func(q * q, p) / sqrt(pi * k_norm)
        result = integrate.quad(f, *limits)[0]  # integral from zero to lower bound
        return result

    return impact_version


@hankel_transform
def imag_gamma(x, p):
	return amplitude(x, p).imag


@hankel_transform
def real_gamma(x, p):
	return -amplitude(x, p).real


def diff_cs(t, p):
	"""
		Differential cross-section formula:
			$$\frac{d\sigma}{dt} = |Re A(s, t)^2 + Im A(s, t)^2|$$
	"""
	A = amplitude(t,p)
	try:
		result = abs(A) ** 2
	except OverflowError:
		result = A.imag
	return result

def ratio(t, p):
    A = amplitude(t,p)
    return A.real/A.imag

class GammaApproximation(object):
    
    def __init__(self, data, new_sigma = None):
        self.dataPoints = data
        self.new_sigma = new_sigma 

        # Define these functions in a very strange way to avoid code replication
        @hankel_transform
        def bspace_low_t(b, p):
        	return self.im_amplitude_low_t(b, p)
        self.bspace_low_t = bspace_low_t


    def im_amplitude_low_t(self, t, p):
        """Calculates imaginary part of extrapolated amplitude"""
        parameters = [i for i in p]

        sigma = self.new_sigma if self.new_sigma else parameters[-2] 

        a0 = sigma / (4 * sqrt(pi * k_norm) )
        a1 = sqrt(self.dataPoints[0].ds - amplitude(self.dataPoints[0].t, p).real ** 2)
        b0 = (1./(self.dataPoints[0].t)) * log(a0 / a1)

        result = a0 * exp(-1. * b0 * t)
        return result

    def __call__(self ,b, p):
        """Calculates gamma(b) using data points"""

        # Taking into account low-t extrapolation
        extrapolation1 = self.bspace_low_t(b, p, (0, self.dataPoints[0].lower ** 0.5))

        # High t-contribution
        extrapolation2 = imag_gamma(b, p, (self.dataPoints[-1].upper ** 0.5, np.infty))
        gamma_data  = 0

        # if B <= 0.01 : print  'At b = ', B, 'ext',extrapolation2

        for i in self.dataPoints:
            A_i = sqrt( fabs(i.ds - amplitude(i.t, p).real ** 2) )
            q1 = sqrt(i.lower)
            q2 = sqrt(i.upper)
            gamma_data = gamma_data + ( (1/sqrt(pi*k_norm))
                               *(
                                    q2*j1(b*q2/k_fm)
                                 -  q1*j1(b*q1/k_fm)
                                )
                               *A_i*(k_fm/b))
        result = extrapolation1 + gamma_data + extrapolation2
        return result