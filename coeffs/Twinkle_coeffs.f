* The frame is m2 centre frame, where m1 = 1/(1+q) at (-s,0), m2 = q/(1+q) at (0,0)
* s means Binary Separation, q means Mass Ratio. s>0, 0<q<=1
* the input zeta_M is the point source coordinate in the mass centre frame, where m1 is at (-s/(1+q),0), m2 is at (sq/(1+q),0)
* "moving" is the length from the mass centre to m2, which is always positive for m1 > m2

* The solution of "coeffs_out" is image positions in m2 centre frame.
* Adding "moving" to the real part of the solution gives the image positions in mass centre frame.

* coeffs_out[1] + coeffs_out[2]*z + coeffs_out[3]*z^2 + coeffs_out[4]*z^3 + coeffs_out[5]*z^4 + coeffs_out[6]*z^5 = 0



COMPLEX FUNCTION norm2(z)
COMPLEX z
REAL :: re, im
re = REAL(z)
im = AIMAG(z)
norm2 = re**2 + im**2
END FUNCTION norm2


subroutine PolyCoeffTwinkle(s, q, zeta_M, coeffs_out, moving)
implicit none
double precision :: s, q, moving
complex*16 :: zeta_M, coeffs_out(6)
double precision :: m1, m2
complex*16 :: temp(3), y, yc, zeta_c

m1 = 1.0d0 / (1.0d0 + q)
m2 = q * m1

moving = s * m1
y = zeta_M - moving
yc = conjg(y)

zeta_c = conjg(zeta_M)

temp(1) = m2 * s
temp(2) = zeta_c + temp(1)
temp(3) = (1.0d0 - s * s) + s * temp(2)

coeffs_out(6) = - yc * temp(2)
coeffs_out(5) = (norm2(y) - 1.0d0 - 2.0d0 * s * yc) * temp(2) + temp(1)
coeffs_out(4) = (2.0d0 * y - s) * (temp(3) * temp(2) + (0.0d0, 2.0d0 * aimag(y))) + (0.0d0, 4.0d0 * aimag(y)) * (temp(1) - y)
coeffs_out(3) = temp(3) * temp(3) * temp(2) * temp(2) + (0.0d0, (1.0d0 - s * y(1)) * aimag(y)) + s * aimag(y) * aimag(y) * (2.0d0 + s * temp(2)) + temp(1) * (2.0d0 * norm2(y) + s * (0.0d0, -y(1) + 3.0d0 * aimag(y)))
coeffs_out(2) = temp(1) * ((s + 2.0d0 * y) * temp(3) + s * ((0.0d0, 2.0d0 * aimag(y) * s) - m2))
coeffs_out(1) = temp(1) * temp(1) * y

return
end subroutine PolyCoeffTwinkle


