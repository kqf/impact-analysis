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

import numpy as np

def amplitude(t_, p):
    k_norm = 0.389379338
    t = t_[0]
    a1 = p[0]
    a2 = p[1]
    a4 = p[2]

    b1 = p[3]
    b2 = p[4]
    b3 = p[5]
    b4 = p[6]

    a5 = p[7]
    b5 = p[8]
    b6 = p[9]
    a_s = p[10]/(sqrt(pi *  k_norm)*4)
    rho = p[11]

    alpha = (1 - 1j*rho)*(a_s + a4)


    try:
        ampl= 1j*alpha*( a1*Exp(-0.5*alpha*b1*t) + (1 - a1)*Exp(-0.5*alpha*b2*t) ) - 1j*a4*Exp(-0.5*b4*t) - a4*rho/((1 + t/b5)**4)
    except OverflowError:
        ampl = 0

    return ampl

def getImage(t, p):
    return amplitude([t], p).imag

def getReal(t, p):
    return amplitude([t], p).real

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
            error_squared += covariance[i][j]*A[i]*A[j]

    error_squared += (dsigma ** 2) * (d_as/(sqrt(pi)*4.))**2 + (drho ** 2) * (d_rho)**2
    error = sqrt(error_squared)
    return error

def getRealGammaError(B, p, covariance, dsigma, drho):
    b = B[0]
    k_fm   = 0.1973269718
    k_norm = 0.389379338

    f = lambda q :  q*j0(b*q/k_fm)*getRealError(q*q, p, covariance, dsigma, drho)/sqrt(pi*k_norm)
    result =  integrate.quad(f, 0, np.infty)[0]  # integral from zero to lower bound

    # imGamma = - \int reAmplitude
    return result

def getRealGamma(B, p):
    b = B[0]
    k_fm   = 0.1973269718
    k_norm = 0.389379338

    f = lambda q :  q*j0(b*q/k_fm)*getReal(q*q, p)/sqrt(pi*k_norm)
    result =  integrate.quad(f, 0, np.infty)[0]  # integral from zero to lower bound
    return -result

def getImagGamma(B, p):
    b = B[0]
    k_fm   = 0.1973269718
    k_norm = 0.389379338

    f = lambda q :  q*j0(b*q/k_fm)*getImage(q*q, p)/sqrt(pi*k_norm)
    result =  integrate.quad(f, 0, np.infty)[0]  # integral from zero to lower bound
    return result

def diff_cs(t, p):
    k_norm = 0.389379338
    A = amplitude(t,p)
    try:
        result =  ( (A.real)**2 + (A.imag)**2 ) # /(k_norm)
    except OverflowError:
        # print 'Real ', A.real, ' Imag ', A.imag
        result = A.imag
        # result =  round(  (A.real)**2 + (A.imag)**2  )
    return result

def ratio(t, p):
    A = amplitude(t,p)
    # r_ = abs(A.real/A.imag )
    r_ = A.imag
    return r_

class GammaApproximation(object):
    k_fm   = 0.1973269718
    k_norm = 0.389379338
    
    def __init__(self, data):
        self.__dataPoints = data
        # Parameters for extrapolation

    def getImAmplitudeExtrNearZero(self, t, p, new_sigma = 0):
        """Calculates imaginary part of extrapolated amplitude"""
        parameters = [ i for i in p]

        # self.__A_0 = sqrt( diff_cs([0], p)  - getReal(0, p)**2 )

        sigma = parameters[-2]

        if int(new_sigma) != 0:
            sigma = new_sigma
            # print 'Using new sigma'

        self.__A_0 = sigma / (4 * sqrt(pi* self.k_norm) )
        self.__A_1 = sqrt( self.__dataPoints[0].ds - getReal(self.__dataPoints[0].t, p)**2 )
        self.__B_0 = ( 1./(self.__dataPoints[0].t) )*log(self.__A_0/self.__A_1)

        result = self.__A_0*exp(-1.*self.__B_0*t)
        return result 

    def getGammaExtrapNearZero(self, b, p, new_sigma = 0): # need to be checked !!!
        """Calculate extrapolation term near zero"""
        B = b[0]

        k_norm = self.k_norm
        k_fm = self.k_fm

        imA = lambda t: self.getImAmplitudeExtrNearZero(t, p, new_sigma)                       #just shortcut
        f = lambda q :  q*j0(B*q/k_fm)*imA(q*q)/sqrt(pi*k_norm)
        result =  integrate.quad(f, 0, self.__dataPoints[0].lower**0.5)[0]  # integral from zero to lower bound
        return result

    def getGammaExtrapNearInf(self, b, p): # need to be checked !!!
        """Calculate extrapolation term near zero"""
        B = b[0]
        k_norm = self.k_norm
        k_fm = self.k_fm

        imA = lambda t: getImage(t, p)
        f = lambda q :  q*j0(B*q/k_fm)*imA(q*q)/sqrt(pi*k_norm)
        # f = lambda q :  q*j0(B*q/k_fm)*getImage(q**2, p)/sqrt(pi*k_norm)
        result =  integrate.quad(f, self.__dataPoints[-1].upper**0.5, np.infty)[0]  # integral from zero to lower bound
        return result


    def gamma(self ,b, p, new_sigma = 0):
        """Gives gamma(b) from points"""
        B = b[0]

        extrapolation1 = self.getGammaExtrapNearZero(b, p, new_sigma) # taking into account extrapolation near zero 
        extrapolation2 = self.getGammaExtrapNearInf(b, p)
        gamma_data  = 0

        # if B <= 0.01 : print  'At b = ', B, 'ext',extrapolation2

        k_norm = self.k_norm
        k_fm = self.k_fm
        for i in self.__dataPoints:
            A_i = sqrt( fabs(i.ds - getReal(i.t, p)**2) )
            q1 = sqrt(i.lower)
            q2 = sqrt(i.upper)
            gamma_data = gamma_data + ( (1/sqrt(pi*k_norm))
                               *(
                                    q2*j1(B*q2/k_fm)
                                 -  q1*j1(B*q1/k_fm)
                                )
                               *A_i*(k_fm/B))
        result = extrapolation1 + gamma_data + extrapolation2
        return result

    def __call__(self, b, p, new_sigma = 0):
      return  self.gamma(b, p, new_sigma)
