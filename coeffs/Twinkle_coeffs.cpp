#include <complex>

using namespace std;

// The frame is m2 centre frame, where m1 = 1/(1+q) at (-s,0), m2 = q/(1+q) at (0,0)
// s means Binary Separation, q means Mass Ratio. s>0, 0<q<=1
// the input zeta_M is the point source coordinate in the mass centre frame, where m1 is at (-s/(1+q),0), m2 is at (sq/(1+q),0)
// "moving" is the length from the mass centre to m2, which is always positive for m1 > m2

// The solution of "coeffs_out" is image positions in m2 centre frame.
// Adding "moving" to the real part of the solution gives the image positions in mass centre frame.

double PolyCoeffTwinkle(const double s, const double q, const complex<double>& zeta_M, complex<double>* coeffs_out)
{
    double m1,m2;
    complex<double> temp[3],y,yc,zeta_c;
    double moving;

    m1 = 1. / (1.+q);
    m2 = q * m1;

    moving = s*m1;
    y = zeta_M - moving;
    yc = conj(y);

    zeta_c = conj(zeta_M);

    temp[0] = m2*s;
    temp[1] = zeta_c + temp[0];
    temp[2] = (1-s*s) + s*temp[1];

    coeffs_out[5] = - yc * temp[1];
    coeffs_out[4] = (norm(y) -1. - 2.*s*yc ) * temp[1] + temp[0];
    coeffs_out[3] = (2.*y-s) * (  temp[2]*temp[1]  + complex<double>(0.,2.*y.imag()) - temp[0] ) + complex<double>(0.,4.*y.imag()) * (temp[0] - y);
    coeffs_out[2] = temp[2].real() * (temp[2].real() * temp[1].real() + complex<double>(0.,(1.-s*y.real())*y.imag())) + s*y.imag()*y.imag()* (2. + s * temp[1]) + temp[0] * ( 2*norm(y) + s * complex<double>(-y.real(), 3*y.imag()) - 1.);
    coeffs_out[1] = temp[0] * ( (s + 2.*y) * temp[2] + s * (complex<double>(0.,2.*y.imag()*s) - m2)  );
    coeffs_out[0] = temp[0] * temp[0] * y;

    return moving;
}


int main()
{
    complex<double> coeffs[6];
    double moving;
    double s = 0.5;
    double q = 0.1;
    complex<double> zeta_M(0.1, 0.2);
    moving = PolyCoeffTwinkle(s, q, zeta_M, coeffs);


    printf("moving: %.16f\n", moving);
    for(int i=0;i<6;i++)
    {
        printf("coeffs[%d]: %.20f + %.20f i\n", i, coeffs[i].real(), coeffs[i].imag());
    }

    printf("norm test: %.3f\n",norm(complex<double>(1,2))   );

    return 0;
}