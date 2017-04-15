from math import sin as Sin
from math import cos as Cos
from math import pow as Power
from math import e as E
from math import sqrt, pi
from impact.constants import k_fm, k_norm
from impact.model import hankel_transform


class RealPartErrorEvaluator(object):
	def __init__(self, covariance, dsigma, drho):
		super(RealPartErrorEvaluator, self).__init__()
		self.covariance = covariance
		self.dsigma = dsigma
		self.drho = drho

		@hankel_transform
		def breal_error(x, p):
			return self.treal_error(x, p)

		self.breal_error = breal_error

	def treal_error(self, t, p):
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
	            error_squared += self.covariance[i][j] * A[i] * A[j]

	    error_squared += (self.dsigma ** 2) * (d_as/(sqrt(pi)*4.)) ** 2 + (self.drho ** 2) * (d_rho) ** 2
	    error = sqrt(error_squared)
	    return error