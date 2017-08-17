from math import sin as Sin
from math import cos as Cos
from math import pow as Power
from math import e as E
from math import sqrt, pi
from impact.constants import k_norm


class PartialExplicit(object):

	def d_a1(self, t, p):
		a1, a2, a4, b1, b2, b3, b4, a5, b5, b6, a_s, rho = p
		a_s = a_s / (sqrt(pi * k_norm) * 4)

		return  (a4 + a_s) * ( 
								rho*( Power(E,((-a4 - a_s)*b1*t)/2.)*Cos(((a4 + a_s)*b1*t*rho)/2.) 
								- Power(E,((-a4 - a_s)*b2*t)/2.)*Cos(((a4 + a_s)*b2*t*rho)/2.) ) 
								- Power(E,((-a4 - a_s)*b1*t)/2.)*Sin(((a4 + a_s)*b1*t*rho)/2.) 
								+ Power(E,((-a4 - a_s)*b2*t)/2.)*Sin(((a4 + a_s)*b2*t*rho)/2.)
						     )

	def d_a4(self, t, p):
		a1, a2, a4, b1, b2, b3, b4, a5, b5, b6, a_s, rho = p
		a_s = a_s / (sqrt(pi * k_norm) * 4)

		return  ( 
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

	# TODO: Fix this contribution
	def d_as(self, t, p):
		a1, a2, a4, b1, b2, b3, b4, a5, b5, b6, a_s, rho = p
		a_s = a_s / (sqrt(pi * k_norm) * 4)

		return ( 
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

	
	def d_b1(self, t, p):
		a1, a2, a4, b1, b2, b3, b4, a5, b5, b6, a_s, rho = p
		a_s = a_s / (sqrt(pi * k_norm) * 4)


		return (a4 + a_s) * (
								- (a1*(a4 + a_s)*Power(E,((-a4 - a_s)*b1*t)/2.)*t*rho*Cos(((a4 + a_s)*b1*t*rho)/2.))/2.
								- (a1*(-a4 - a_s)*Power(E,((-a4 - a_s)*b1*t)/2.)*t*Sin(((a4 + a_s)*b1*t*rho)/2.))/2. 
								+ rho*((a1*(-a4 - a_s)*Power(E,((-a4 - a_s)*b1*t)/2.)*t*Cos(((a4 + a_s)*b1*t*rho)/2.))/2.
								- (a1*(a4 + a_s)*Power(E,((-a4 - a_s)*b1*t)/2.)*t*rho*Sin(((a4 + a_s)*b1*t*rho)/2.))/2.)
						    )
	 

	def d_b2(self, t, p):
		a1, a2, a4, b1, b2, b3, b4, a5, b5, b6, a_s, rho = p
		a_s = a_s / (sqrt(pi * k_norm) * 4)
		return  (a4 + a_s) * (
								 - ((1 - a1)*(a4 + a_s)*Power(E,((-a4 - a_s)*b2*t)/2.)*t*rho*Cos(((a4 + a_s)*b2*t*rho)/2.))/2. 
								 - ((1 - a1)*(-a4 - a_s)*Power(E,((-a4 - a_s)*b2*t)/2.)*t*Sin(((a4 + a_s)*b2*t*rho)/2.))/2. 
								 + rho*( ((1 - a1)*(-a4 - a_s)*Power(E,((-a4 - a_s)*b2*t)/2.)*t*Cos(((a4 + a_s)*b2*t*rho)/2.))/2. 
								 - ((1 - a1)*(a4 + a_s)*Power(E,((-a4 - a_s)*b2*t)/2.)*t*rho*Sin(((a4 + a_s)*b2*t*rho)/2.))/2.)
						     )

	def d_b5(self, t, p):
		a1, a2, a4, b1, b2, b3, b4, a5, b5, b6, a_s, rho = p
		a_s = a_s / (sqrt(pi * k_norm) * 4)

		return (-4*a4*t*rho)/(Power(b5,2)*Power(1 + t/b5,5))


	def d_rho(self, t, p):
		a1, a2, a4, b1, b2, b3, b4, a5, b5, b6, a_s, rho = p
		a_s = a_s / (sqrt(pi * k_norm) * 4)

		return ( -(a4/Power(1 + t/b5,4)) + 
				   (a4 + a_s) * (
									a1*Power(E,((-a4 - a_s)*b1*t)/2.)*Cos(((a4 + a_s)*b1*t*rho)/2.)
									- (a1*(a4 + a_s)*b1*Power(E,((-a4 - a_s)*b1*t)/2.)*t*Cos(((a4 + a_s)*b1*t*rho)/2.))/2.
									+ (1 - a1)*Power(E,((-a4 - a_s)*b2*t)/2.)*Cos(((a4 + a_s)*b2*t*rho)/2.)
									- ((1 - a1)*(a4 + a_s)*b2*Power(E,((-a4 - a_s)*b2*t)/2.)*t*Cos(((a4 + a_s)*b2*t*rho)/2.))/2.
									+ rho*(-(a1*(a4 + a_s)*b1*Power(E,((-a4 - a_s)*b1*t)/2.)*t*Sin(((a4 + a_s)*b1*t*rho)/2.))/2.
									- ((1 - a1)*(a4 + a_s)*b2*Power(E,((-a4 - a_s)*b2*t)/2.)*t*Sin(((a4 + a_s)*b2*t*rho)/2.))/2.)
							 	)
				)

