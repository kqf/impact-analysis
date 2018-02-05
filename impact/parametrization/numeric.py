from math import sin as Sin
from math import cos as Cos
from math import pow as Power
from math import e as E
from math import sqrt, pi
from cmath import exp as Exp

from impact.constants import k_norm

class Numeric(object):

    def d_a1(self, t, p):
        a1, a2, b1, b2, b3, b4, a_s, rho = p
        a_s = a_s / (sqrt(pi * k_norm) * 4)

        return  (a2 + a_s) * ( 
                                rho*( Power(E,((-a2 - a_s)*b1*t)/2.)*Cos(((a2 + a_s)*b1*t*rho)/2.) 
                                - Power(E,((-a2 - a_s)*b2*t)/2.)*Cos(((a2 + a_s)*b2*t*rho)/2.) ) 
                                - Power(E,((-a2 - a_s)*b1*t)/2.)*Sin(((a2 + a_s)*b1*t*rho)/2.) 
                                + Power(E,((-a2 - a_s)*b2*t)/2.)*Sin(((a2 + a_s)*b2*t*rho)/2.)
                             )

    def d_a2(self, t, p):
        a1, a2, b1, b2, b3, b4, a_s, rho = p
        a_s = a_s / (sqrt(pi * k_norm) * 4)

        return  ( 
                    - (rho/Power(1 + t/b4,4)) 
                    + rho*(   a1*Power(E,((-a2 - a_s)*b1*t)/2.)*Cos(((a2 + a_s)*b1*t*rho)/2.) 
                    + (1 - a1)*Power(E,((-a2 - a_s)*b2*t)/2.)*Cos(((a2 + a_s)*b2*t*rho)/2.) ) 
                    - a1*Power(E,((-a2 - a_s)*b1*t)/2.)*Sin(((a2 + a_s)*b1*t*rho)/2.) 
                    - (1 - a1)*Power(E,((-a2 - a_s)*b2*t)/2.)*Sin(((a2 + a_s)*b2*t*rho)/2.) 

                    + (a2 + a_s)
                    *(
                        -(a1*b1*Power(E,((-a2 - a_s)*b1*t)/2.)*t*rho*Cos(((a2 + a_s)*b1*t*rho)/2.))/2. 
                        - ((1 - a1)*b2*Power(E,((-a2 - a_s)*b2*t)/2.)*t*rho*Cos(((a2 + a_s)*b2*t*rho)/2.))/2. 
                        + (a1*b1*Power(E,((-a2 - a_s)*b1*t)/2.)*t*Sin(((a2 + a_s)*b1*t*rho)/2.))/2. 
                        + ((1 - a1)*b2*Power(E,((-a2 - a_s)*b2*t)/2.)*t*Sin(((a2 + a_s)*b2*t*rho)/2.))/2. 
                        + rho*( - (a1*b1*Power(E,((-a2 - a_s)*b1*t)/2.)*t*Cos(((a2 + a_s)*b1*t*rho)/2.))/2. 
                                - ((1 - a1)*b2*Power(E,((-a2 - a_s)*b2*t)/2.)*t*Cos(((a2 + a_s)*b2*t*rho)/2.))/2. 
                                - (a1*b1*Power(E,((-a2 - a_s)*b1*t)/2.)*t*rho*Sin(((a2 + a_s)*b1*t*rho)/2.))/2. 
                                - ((1 - a1)*b2*Power(E,((-a2 - a_s)*b2*t)/2.)*t*rho*Sin(((a2 + a_s)*b2*t*rho)/2.))/2. 
                              )
                      )
                )

    def d_as(self, t, p):
        a1, a2, b1, b2, b3, b4, a_s, rho = p
        a_s = a_s / (sqrt(pi * k_norm) * 4)
        alpha = (1 - 1j*rho)*(a_s + a2)
        ralpha = (1 - 1j*rho) / (sqrt(pi * k_norm) * 4)

        first = 1j * (a1 * Exp(-0.5 * alpha * b1 * t) + (1 - a1) * Exp(-0.5 * alpha * b2 * t)) * ralpha 
        second = 1j * alpha * (b1 * a1 * Exp(-0.5 * alpha * b1 * t) + (1 - a1) * b2 * Exp(-0.5 * alpha * b2 * t)) \
            * ralpha * (-0.5 * t)

        result = first + second
        return result.real

    
    def d_b1(self, t, p):
        a1, a2, b1, b2, b3, b4, a_s, rho = p
        a_s = a_s / (sqrt(pi * k_norm) * 4)


        return (a2 + a_s) * (
                                - (a1*(a2 + a_s)*Power(E,((-a2 - a_s)*b1*t)/2.)*t*rho*Cos(((a2 + a_s)*b1*t*rho)/2.))/2.
                                - (a1*(-a2 - a_s)*Power(E,((-a2 - a_s)*b1*t)/2.)*t*Sin(((a2 + a_s)*b1*t*rho)/2.))/2. 
                                + rho*((a1*(-a2 - a_s)*Power(E,((-a2 - a_s)*b1*t)/2.)*t*Cos(((a2 + a_s)*b1*t*rho)/2.))/2.
                                - (a1*(a2 + a_s)*Power(E,((-a2 - a_s)*b1*t)/2.)*t*rho*Sin(((a2 + a_s)*b1*t*rho)/2.))/2.)
                            )
     

    def d_b2(self, t, p):
        a1, a2, b1, b2, b3, b4, a_s, rho = p
        a_s = a_s / (sqrt(pi * k_norm) * 4)
        return  (a2 + a_s) * (
                                 - ((1 - a1)*(a2 + a_s)*Power(E,((-a2 - a_s)*b2*t)/2.)*t*rho*Cos(((a2 + a_s)*b2*t*rho)/2.))/2. 
                                 - ((1 - a1)*(-a2 - a_s)*Power(E,((-a2 - a_s)*b2*t)/2.)*t*Sin(((a2 + a_s)*b2*t*rho)/2.))/2. 
                                 + rho*( ((1 - a1)*(-a2 - a_s)*Power(E,((-a2 - a_s)*b2*t)/2.)*t*Cos(((a2 + a_s)*b2*t*rho)/2.))/2. 
                                 - ((1 - a1)*(a2 + a_s)*Power(E,((-a2 - a_s)*b2*t)/2.)*t*rho*Sin(((a2 + a_s)*b2*t*rho)/2.))/2.)
                             )

    def d_b4(self, t, p):
        a1, a2, b1, b2, b3, b4, a_s, rho = p
        a_s = a_s / (sqrt(pi * k_norm) * 4)

        return (-4*a2*t*rho)/(Power(b4,2)*Power(1 + t/b4,5))


    def d_rho(self, t, p):
        a1, a2, b1, b2, b3, b4, a_s, rho = p
        a_s = a_s / (sqrt(pi * k_norm) * 4)

        return ( -(a2/Power(1 + t/b4,4)) + 
                   (a2 + a_s) * (
                                    a1*Power(E,((-a2 - a_s)*b1*t)/2.)*Cos(((a2 + a_s)*b1*t*rho)/2.)
                                    - (a1*(a2 + a_s)*b1*Power(E,((-a2 - a_s)*b1*t)/2.)*t*Cos(((a2 + a_s)*b1*t*rho)/2.))/2.
                                    + (1 - a1)*Power(E,((-a2 - a_s)*b2*t)/2.)*Cos(((a2 + a_s)*b2*t*rho)/2.)
                                    - ((1 - a1)*(a2 + a_s)*b2*Power(E,((-a2 - a_s)*b2*t)/2.)*t*Cos(((a2 + a_s)*b2*t*rho)/2.))/2.
                                    + rho*(-(a1*(a2 + a_s)*b1*Power(E,((-a2 - a_s)*b1*t)/2.)*t*Sin(((a2 + a_s)*b1*t*rho)/2.))/2.
                                    - ((1 - a1)*(a2 + a_s)*b2*Power(E,((-a2 - a_s)*b2*t)/2.)*t*Sin(((a2 + a_s)*b2*t*rho)/2.))/2.)
                                )
                )
        
    def amplitude(self, t, p):
        a1, a2, b1, b2, b3, b4, a_s, rho = p
        a_s = a_s/(sqrt(pi *  k_norm) * 4)
        
        alpha = (1 - 1j*rho)*(a_s + a2)

        try:
            ampl= 1j*alpha*( a1*Exp(-0.5*alpha*b1*t) + (1 - a1)*Exp(-0.5*alpha*b2*t) ) - 1j*a2*Exp(-0.5*b3*t) - a2*rho/((1 + t/b4)**4)
        except OverflowError:
            ampl = 0
        return ampl 