#!/usr/bin/python2.7

from cmath import exp as Exp
from math  import exp, pi, fabs, sqrt, log
from scipy import integrate
from scipy.special import j0, j1
# from array  import array
from random import randrange

def amplitude(t, p):
    A1 = p[0]
    A2 = p[1]
    A3 = p[2]
    b1 = p[3]
    b2 = p[4]
    b3 = p[5]
    sigma = p[6]              # from tamble 2
    sigma0 = p[7] 
    rho = p[8]               # from tamble 2
    j = complex(0, 1)

    alpha =  (sigma/sigma0)*(1 - j*rho)
    secondTerm = A1*Exp(-0.5*b1*t[0]*alpha) + A2*Exp(-0.5*b2*t[0]*alpha)
    # print 'second term ', A2*Exp(-0.5*b2*t[0]*alpha)
    thirdTerm = -A3*j*Exp(-0.5*b3*t[0])
    ampl = j*alpha * secondTerm + thirdTerm

    return ampl

def getImage(t, p):
    return amplitude([t], p).imag

def getReal(t, p):
    return amplitude([t], p).real

def diff_cs(t, p):
    A = amplitude(t,p)
    try:
        result =  ( (A.real)**2 + (A.imag)**2 )
    except OverflowError:
        print 'Real ', A.real, ' Imag ', A.imag
        result =  round(  (A.real)**2 + (A.imag)**2  )
    return result

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

        # extrapolation1 = self.getGammaExtrapNearZero(b, p) # taking into account extrapolation near zero 
        # result = extrapolation1

        k_norm = self.k_norm
        k_fm = self.k_fm
        # for i in self.__dataPoints:
            # A_i = sqrt( fabs(i.ds - getReal(i.t, p)**2) )
            # q1 = sqrt(i.lower)
            # q2 = sqrt(i.upper)
            # result = result + ( (1/sqrt(pi*k_norm))
                               # *(
                                    # q2*j1(B*q2/k_fm)
                                 # -  q1*j1(B*q1/k_fm)
                                # )
                               # *A_i*(k_fm/B))

         imA = lambda t: getImage(t, p)
         f = lambda q :  q*j0(B*q/k_fm)*imA(q*q)/sqrt(pi*k_norm)
         result =  integrate.quad(f, 0, self.__dataPoints[-1].upper**0.5)[0]  # integral from zero to lower bound
         return result




