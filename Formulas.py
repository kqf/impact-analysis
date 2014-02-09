#!/usr/bin/python2.7

from cmath import exp as Exp
from math  import exp, pi, fabs, sqrt, log, sin, cos
from scipy import integrate
from scipy.special import j0, j1
# from array  import array
from random import randrange

def amplitude(t_, p):
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
    a_s = p[10]/(sqrt(pi)*4)
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

from math import sin as Sin
from math import cos as Cos
from math import pow as Power
from math import e as E

def getRealError(t, p, pe):
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
    a_s = p[10]/(sqrt(pi)*4)
    rho = p[11]
 
    # errors
    a1e = pe[0]
    a2e = pe[1]
    a4e = pe[2]

    b1e = pe[3]
    b2e = pe[4]
    b3e = pe[5]
    b4e = pe[6]

    a5e = pe[7]
    b5e = pe[8]
    b6e = pe[9]
    a_se = pe[10]/(sqrt(pi)*4)
    rhoe = pe[11]

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
    error = sqrt(
                      (d_a1**2)*(a1e**2)
                    + (d_a4**2)*(a4e**2)
                    + (d_as**2)*(a_se**2)
                    + (d_b1**2)*(b1e**2)
                    + (d_b2**2)*(b2e**2)
                    + (d_b5**2)*(b5e**2)
                    + (d_rho**2)*(rhoe**2)
                )

    return error

def getRealGammaError(b, t_max, p, pe):
    k_fm   = 0.1973269718
    k_norm = 0.389379338

    f = lambda q :  q*j0(b*q/k_fm)*getRealError(q*q, p, pe)/sqrt(pi*k_norm)
    result =  integrate.quad(f, 0, t_max**0.5)[0]  # integral from zero to lower bound

    return result

def getRealGamma(B, t_max, p):
    b = B[0]
    k_fm   = 0.1973269718
    k_norm = 0.389379338

    f = lambda q :  q*j0(b*q/k_fm)*getReal(q*q, p)/sqrt(pi*k_norm)
    result =  integrate.quad(f, 0, t_max**0.5)[0]  # integral from zero to lower bound

    return result

def diff_cs(t, p):
    k_norm = 0.389379338
    A = amplitude(t,p)
    try:
        result =  ( (A.real)**2 + (A.imag)**2 )/k_norm
    except OverflowError:
        print 'Real ', A.real, ' Imag ', A.imag
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

    def getImAmplitudeExtrNearZero(self, t, p):
        """Calculates imaginary part of extrapolated amplitude"""
        self.__A_0 = sqrt( diff_cs([0], p)  - getReal(0, p)**2 )
        self.__A_1 = sqrt( self.__dataPoints[0].ds - getReal(self.__dataPoints[0].t, p)**2 )
        self.__B_0 = ( 1./(self.__dataPoints[0].t) )*log(self.__A_0/self.__A_1)

        result = self.__A_0*exp(-1.*self.__B_0*t)
        return result 


    def getGammaExtrapNearZero(self, b, p): # need to be checked !!!
        """Calculate extrapolation term near zero"""
        B = b[0]

        k_norm = self.k_norm
        k_fm = self.k_fm
        imA = lambda t: self.getImAmplitudeExtrNearZero(t, p)                       #just shortcut
        # imA = lambda t: getImage(t, p)
        f = lambda q :  q*j0(B*q/k_fm)*imA(q*q)/sqrt(pi*k_norm)
        # f = lambda q :  q*j0(B*q/k_fm)*getImage(q**2, p)/sqrt(pi*k_norm)
        result =  integrate.quad(f, 0, self.__dataPoints[0].lower**0.5)[0]  # integral from zero to lower bound
        return result

    def gamma(self ,b, p):
        """Gives gamma(b) from points"""
        B = b[0]

        extrapolation1 = self.getGammaExtrapNearZero(b, p) # taking into account extrapolation near zero 
        result = extrapolation1

        k_norm = self.k_norm
        k_fm = self.k_fm
        for i in self.__dataPoints:
            A_i = sqrt( fabs(i.ds - getReal(i.t, p)**2) )
            q1 = sqrt(i.lower)
            q2 = sqrt(i.upper)
            result = result + ( (1/sqrt(pi*k_norm))
                               *(
                                    q2*j1(B*q2/k_fm)
                                 -  q1*j1(B*q1/k_fm)
                                )
                               *A_i*(k_fm/B))
        return result
