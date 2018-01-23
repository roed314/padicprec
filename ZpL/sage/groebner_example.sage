#########################
#    Groebner bases     #
#########################

from sage.rings.polynomial.toy_buchberger import *
A = QpLC(2,50,print_mode='terse')
prec = A.precision()
R.<x,y,z> = PolynomialRing(A, order = 'invlex')
F=[2*x + z, x^2 + y^2 - 2*z^2,4*y^2 + y*z + 8*z^2]
I = ideal(F)
g1 = buchberger_improved(I)
g1.sort()

A = Qp(2,50,print_mode='terse')
R.<x,y,z> = PolynomialRing(A, order = 'invlex')
F=[2*x + z, x^2 + y^2 - 2*z^2,4*y^2 + y*z + 8*z^2]
I = ideal(F)
g2 = buchberger_improved(I)
g2.sort()

#sage: g1
#[x^3, x*y + (2251799813685218 + O(2^51))*x^2, y^2 + (1125899906842617 + 
#O(2^50))*x^2, z + (2 + O(2^51))*x]

#sage: g2
#[x^3, x*y + (1125899906842594 + O(2^50))*x^2, y^2 + (281474976710649 + 
#O(2^48))*x^2, z + (2 + O(2^51))*x]
