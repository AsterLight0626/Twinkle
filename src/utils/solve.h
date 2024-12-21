#pragma once

#include <type_traits>

namespace twinkle::root_finder
{

static const float ROUND_OFF_D=2e-15f;
static const float ROUND_OFF_F=1e-6f;

static const int MAXIT = 80;
static const int FRAC_JUMP_EVERY = 16;      // must be 2^n

template< class c_T >
__device__ bool cmplx_roots_gen(c_T *roots, c_T *poly, int degree, bool polish_roots_after, bool use_roots_as_starting_points);
template< class img_T, class c_T >
__device__ bool cmplx_roots_gen( img_T *imgs, c_T *poly, int degree, bool polish_roots_after, bool use_roots_as_starting_points);
template<class c_T>
__device__ void solve_quadratic_eq(c_T &x0, c_T &x1, c_T *poly);
template< class c_T >
__device__ void cmplx_newton_spec(c_T *poly, int degree, c_T *root, int &iter, bool &success);
template<class c_T>
__device__ void cmplx_laguerre(c_T *poly, int degree, c_T *root, int &iter, bool &success);
template< class c_T >
__device__ void cmplx_laguerre2newton(c_T *poly, int degree, c_T *root, int &iter, bool &success, int starting_mode);

template< class c_T >
__device__ __forceinline__ auto abs_c( const c_T & z )
{
    return z.abs(  );
};
template< class c_T >
__device__ __forceinline__ auto real ( const c_T & z )
{
    return z.re;
};
template< class c_T >
__device__ __forceinline__ auto imag ( const c_T & z )
{
    return z.im;
};
template< class c_T >
__device__ __forceinline__ auto conj ( const c_T & z )
{
    return z.conj(  );
};
template< class c_T >
__device__ __forceinline__ auto sqrt ( const c_T & z )
{
    return z.sqrt(  );
};
template< class c_T >
__device__ __forceinline__ auto norm ( const c_T & z )
{
    return z.norm2(  );
};


template< class c_T >
__device__ bool cmplx_roots_gen(c_T *roots, c_T *poly, int degree, bool polish_roots_after, bool use_roots_as_starting_points)
{
    using f_T = typename c_T::f_t;
    
	//roots - array which will hold all roots that had been found.
	//If the flag 'use_roots_as_starting_points' is set to
	//.true., then instead of point(0, 0) we use value from
	//this array as starting point for cmplx_laguerre

	//poly - is an array of polynomial cooefs, length = degree + 1,
	//poly[0] x ^ 0 + poly[1] x ^ 1 + poly[2] x ^ 2 + ...

	//degree - degree of the polynomial and size of 'roots' array

	//polish_roots_after - after all roots have been found by dividing
	//original polynomial by each root found,
	//you can opt in to polish all roots using full
	//polynomial

	//use_roots_as_starting_points - usually we start Laguerre's 
	//method from point(0, 0), but you can decide to use the
	//values of 'roots' array as starting point for each new
	//root that is searched for.This is useful if you have
	//very rough idea where some of the roots can be.
	//

	// c_T* poly2 = new c_T[degree+1];
	c_T poly2[6];
	// static int i, j, n, iter;
    int iter;
	bool success;
	c_T coef, prev;

	if (!use_roots_as_starting_points) {
		for (int jj = 0; jj < degree; jj++)
            roots[ jj ].set( 0, 0 );
        // {
		// 	roots[jj].re = 0;
		// 	roots[jj].im = 0; //c_T( 0 );
		// }
	}

	for (int j = 0; j <= degree; j++){ poly2[j] = poly[j]; }

	// Don't do Laguerre's for small degree polynomials
	if (degree <= 1) {
		if (degree == 1) roots[0] = -poly[0] / poly[1];
		return false;
	}

	for (int n = degree; n >= 3; n--) {
		cmplx_laguerre2newton(poly2, n, &roots[n - 1], iter, success, 2);        // iter means iteration steps, 2 is "starting mode"
		if (!success) {
			roots[n - 1].set(0, 0);
			cmplx_laguerre(poly2, n, &roots[n - 1], iter, success);
		}

		// Divide by root
		coef = poly2[n];
		for (int i = n - 1; i >= 0; i--) {
			prev = poly2[i];
			poly2[i] = coef;
			coef *= prev + roots[n - 1];
		}
	}


	//Find the to last 2 roots
	solve_quadratic_eq(roots[1], roots[0], poly2);
	//cmplx_laguerre2newton(poly2, 2, &roots[1], iter, success, 2);
	//if (!success) {
	//	solve_quadratic_eq(roots[1], roots[0], poly2);
	//}
	//else {
	//	roots[0] = -(roots[1] + poly2[1] / poly2[2]); // Viete's Formula for the last root
	//}



	if (polish_roots_after) {
		for (int n = 0; n < degree; n++) {
			cmplx_newton_spec(poly, degree, &roots[n], iter, success); // Polish roots with full polynomial
		}
	}

	c_T& sum = coef;
	c_T& multi = prev;
	// use Viete's Formula to test if get repeated solution
	sum.set(0,0);
	multi.set(1,0);
	bool failed( false );
	// sign = ((degree & 1)==1)?-1:1; // degree is odd: negative, even: positive
	
	for(int j=0;j<degree;j++)
	{
		sum   += roots[ j ];
		multi += roots[ j ];
	}


	if constexpr ( std::is_same_v< f_T, double > )
	{
		if((degree & 1)==1)
		{
			failed = (norm(multi + poly[0]/poly[degree])<1e-10f);
		}
		else
		{
			failed = (norm(multi - poly[0]/poly[degree])<1e-10f);
		}
		failed = !(norm(sum + poly[degree-1]/poly[degree])<1e-10  && failed);

	}
	else
	{
		if constexpr (std::is_same_v< f_T, float >)
		{
			if((degree & 1)==1)
			{
				failed = (norm(multi + poly[0]/poly[degree])<1e-5f);
			}
			else
			{
				failed = (norm(multi - poly[0]/poly[degree])<1e-5f);
			}
			failed = !(norm(sum + poly[degree-1]/poly[degree])<1e-5f  && failed);

			// if(! (norm(sum + poly[degree-1]/poly[degree])<1e-5 && norm(multi - sign*poly[0]/poly[degree])<1e-5))
			// {
			// 	failed = true;
			// }
		}
		else{failed = true;}		
	}


    // delete[] poly2;

	return failed;
}

template< class img_T, class c_T >
__device__ bool cmplx_roots_gen( img_T *imgs, c_T *poly, int degree, bool polish_roots_after, bool use_roots_as_starting_points)
{
    using f_T = typename c_T::f_t;
    
	//roots - array which will hold all roots that had been found.
	//If the flag 'use_roots_as_starting_points' is set to
	//.true., then instead of point(0, 0) we use value from
	//this array as starting point for cmplx_laguerre

	//poly - is an array of polynomial cooefs, length = degree + 1,
	//poly[0] x ^ 0 + poly[1] x ^ 1 + poly[2] x ^ 2 + ...

	//degree - degree of the polynomial and size of 'roots' array

	//polish_roots_after - after all roots have been found by dividing
	//original polynomial by each root found,
	//you can opt in to polish all roots using full
	//polynomial

	//use_roots_as_starting_points - usually we start Laguerre's 
	//method from point(0, 0), but you can decide to use the
	//values of 'roots' array as starting point for each new
	//root that is searched for.This is useful if you have
	//very rough idea where some of the roots can be.
	//

	// c_T* poly2 = new c_T[degree+1];
	c_T poly2[6];
	// static int i, j, n, iter;
    int iter;
	bool success;
	c_T coef, prev;

	if (!use_roots_as_starting_points) {
		for (int jj = 0; jj < degree; jj++) {
			imgs[jj].position.set(0, 0);
		}
	}

	for (int j = 0; j <= degree; j++){ poly2[j] = poly[j]; }

	// Don't do Laguerre's for small degree polynomials
	if (degree <= 1) {
		if (degree == 1) imgs[0].position = -poly[0] / poly[1];
		return false;
	}

	for (int n = degree; n >= 3; n--) {
		cmplx_laguerre2newton(poly2, n, &(imgs[n - 1].position), iter, success, 2);        // iter means iteration steps, 2 is "starting mode"
		if (!success) {
			imgs[n - 1].position.set(0, 0);
			cmplx_laguerre(poly2, n, &(imgs[n - 1].position), iter, success);
		}

		// Divide by root
		coef = poly2[n];
		for (int i = n - 1; i >= 0; i--) {
			prev = poly2[i];
			poly2[i] = coef;
			coef = prev + imgs[n - 1].position * coef;
		}
	}


	//Find the to last 2 roots
	solve_quadratic_eq(imgs[1].position, imgs[0].position, poly2);
	//cmplx_laguerre2newton(poly2, 2, &roots[1], iter, success, 2);
	//if (!success) {
	//	solve_quadratic_eq(roots[1], roots[0], poly2);
	//}
	//else {
	//	roots[0] = -(roots[1] + poly2[1] / poly2[2]); // Viete's Formula for the last root
	//}



	if (polish_roots_after) {
		for (int n = 0; n < degree; n++) {
			cmplx_newton_spec(poly, degree, &(imgs[n].position), iter, success); // Polish roots with full polynomial
		}
	}

	c_T& sum = coef;
	c_T& multi = prev;
	// use Viete's Formula to test if get repeated solution
	sum.set(0,0);
	multi.set(1,0);
	bool failed( false );
	// sign = ((degree & 1)==1)?-1:1; // degree is odd: negative, even: positive
	
	for(int j=0;j<degree;j++)
	{
		sum = sum + imgs[j].position;
		multi = multi * imgs[j].position;
	}


	if constexpr ( std::is_same_v< f_T, double > )
	{
		if((degree & 1)==1)
		{
			failed = (norm(multi + poly[0]/poly[degree])<1e-10f);
		}
		else
		{
			failed = (norm(multi - poly[0]/poly[degree])<1e-10f);
		}
		failed = !(norm(sum + poly[degree-1]/poly[degree])<1e-10  && failed);

	}
	else
	{
        if constexpr ( std::is_same_v< f_T, float > )
		{
			if((degree & 1)==1)
			{
				failed = (norm(multi + poly[0]/poly[degree])<1e-5f);
			}
			else
			{
				failed = (norm(multi - poly[0]/poly[degree])<1e-5f);
			}
			failed = !(norm(sum + poly[degree-1]/poly[degree])<1e-5  && failed);

			// if(! (norm(sum + poly[degree-1]/poly[degree])<1e-5 && norm(multi - sign*poly[0]/poly[degree])<1e-5))
			// {
			// 	failed = true;
			// }
		}
		else{failed = true;}		
	}


    // delete[] poly2;

	return failed;
}

template<class c_T>
__device__ void solve_quadratic_eq(c_T &x0, c_T &x1, c_T *poly)
{
    using f_T = typename c_T::f_t;
    
    c_T delta;
	c_T& a = poly[2];
	c_T& b = poly[1];
	c_T& c = poly[0];
	delta = sqrt(b*b - f_T(4) * a*c);
	if (real(conj(b)*delta) >= 0) {
		x0 = f_T(-0.5f)*(b + delta);
	}
	else {
		x0 = f_T(-0.5f)*(b - delta);
	}
	if (x0 == c_T(0, 0)) {
		x1.set(0, 0);
	}
	else { //Viete's formula
		x1 = c / x0;
		x0 = x0 / a;
	}
}

template< class c_T >
__device__ void cmplx_newton_spec(c_T *poly, int degree, c_T *root, int &iter, bool &success)
{
    using f_T = typename c_T::f_t;    
	//Subroutine finds one root of a c_T polynomial
	//Newton's method. It calculates simplified Adams' stopping 
	//criterion for the value of the polynomial once per 10 iterations (!),
	//after initial iteration. This is done to speed up calculations
	//when polishing roots that are known preety well, and stopping
	// criterion does significantly change in their neighborhood.

	//Uses 'root' value as a starting point (!!!!!)
	//Remember to initialize 'root' to some initial guess.
	//Do not initilize 'root' to point (0,0) if the polynomial 
	//coefficients are strictly real, because it will make going 
	//to imaginary roots impossible.

	// poly - is an array of polynomial cooefs
	//	length = degree+1, poly(1) is constant 
	//0					1				2
	//poly[0] x^0 + poly[1] x^1 + poly[2] x^2 + ...
	//degree - a degree of the polynomial
	// root - input: guess for the value of a root
	//		  output: a root of the polynomial
	//iter - number of iterations performed (the number of polynomial evaluations)
	//success - is false if routine reaches maximum number of iterations

	//For a summary of the method go to: 
	//http://en.wikipedia.org/wiki/Newton's_method

    // 关于jump：太占内存，回头再说
	// int FRAC_JUMP_EVERY = 10;
	// const int FRAC_JUMP_LEN = 10;
	// double FRAC_JUMPS[FRAC_JUMP_LEN] = { 0.64109297, 0.91577881, 0.25921289, 0.50487203, 0.08177045, 0.13653241, 0.306162, 0.37794326, 0.04618805, 0.75132137 }; //some random numbers
	c_T faq; //jump length

    f_T FRAC_ERR;
	if constexpr ( std::is_same_v< f_T, float > )    
    {
        FRAC_ERR = ROUND_OFF_F;
    }
    else{FRAC_ERR = ROUND_OFF_D;}
	c_T p; //value of polynomial
	c_T dp; //value of 1st derivative
	int i, k;
	bool good_to_go;
	c_T dx, newroot;
	f_T ek, absroot, abs2p;
	c_T zero(0, 0);
	f_T stopping_crit2;

	iter = 0;
	success = true;

	// //the next if block is an EXTREME failsafe, not usually needed, and thus turned off in this version
	// if (false) { //change false to true if you would like to use caustion about haveing first coefficient == 0
	// 	if (degree < 0) {
	// 		printf("Error: cmplx_newton_spec: degree<0");
	// 		return;
	// 	}
	// 	if (poly[degree] == zero) {
	// 		if (degree == 0) return;
	// 		cmplx_newton_spec(poly, degree, root, iter, success);
	// 		return;
	// 	}
	// 	if (degree <= 1) {
	// 		if (degree == 0) {
	// 			success = false;
	// 			printf("Warning: cmplx_newton_spec: degree=0 and poly[0]!=0, no roots");
	// 			return;
	// 		}
	// 		else {
	// 			*root = -poly[0] / poly[1];
	// 			return;
	// 		}
	// 	}
	// }
	// //end EXTREME Failsafe


	good_to_go = false;

	stopping_crit2 = 0; //value not important, will be initialized anyway on the first loop
	for (i = 1; i <= MAXIT; i++) {
		//prepare stoping criterion
		//calculate value of polynomial and its first two derivatives
		p = poly[degree];
		dp = zero;
		if ((i&7) == 1) { //calculate stopping criterion every tenth iteration
			ek = abs_c(poly[degree]);
			absroot = abs_c(*root);
			for (k = degree - 1; k >= 0; k--) {
				dp = p + dp * (*root);
				p = poly[k] + p * (*root); //b_k
										   //Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
										   //Communications of ACM, Volume 10 Issue 10, Oct. 1967, p. 655
										   //ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
										   //Eq. 8
				ek = absroot * ek + abs_c(p);
			}
			// stopping_crit2 = pow(FRAC_ERR * ek, 2);
            stopping_crit2 = FRAC_ERR * ek * FRAC_ERR * ek;
		}
		else { // calculate just the value and derivative
			for (k = degree - 1; k >= 0; k--) { //Horner Scheme, see for eg. Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
				dp = p + dp * (*root);
				p = poly[k] + p * (*root);
			}
		}

		iter = iter + 1;

		abs2p = real(conj(p) * p);
		if (abs2p == 0) return;
		if (abs2p < stopping_crit2) { //simplified a little Eq. 10 of Adams 1967
			if (dp == zero) return; //if we have problem with zero, but we are close to the root, just accept
									//do additional iteration if we are less than 10x from stopping criterion
			if (abs2p < 0.01f * stopping_crit2) return; //return immediatley because we are at very good place
			else {
				good_to_go = true; //do one iteration more
			}
		}

		else {
			good_to_go = false; //reset if we are outside the zone of the root
		}
		if (dp == zero) {
			//problem with zero
			// dx = (abs(*root) + 1) * expcmplx(c_T(0, FRAC_JUMPS[i% FRAC_JUMP_LEN] * 2 * M_PI));  // just some random numbers is enough
            // dx = f_T(abs_c(*root) + 1) * c_T(0.626f, 0.212f);        // some random () numbers
            dx.set( 0.626f, 0.212f );
            dx *= ( 1 + abs_c( * root ) );
		}
		else {
			dx = p / dp; // Newton method, see http://en.wikipedia.org/wiki/Newton's_method
		}
		newroot = *root - dx;
		if (newroot == *root) return; //nothing changes -> return
		if (good_to_go) {//this was jump already after stopping criterion was met
			*root = newroot;
			return;
		}
		if ((i & (FRAC_JUMP_EVERY - 1 ))== 0) { // decide whether to do a jump of modified length (to break cycles)
			faq.set(0.108f,0.525f);  // some random () numbers
			newroot = *root - faq * dx;
		}
		*root = newroot;
	}
	success = false;
	//too many iterations here
}

template<class c_T>
__device__ void cmplx_laguerre(c_T *poly, int degree, c_T *root, int &iter, bool &success)
{
    using f_T = typename c_T::f_t;    
	//Subroutine finds one root of a c_T polynomial using
	//Laguerre's method. In every loop it calculates simplified 
	//Adams' stopping criterion for the value of the polynomial.
	//
	//Uses 'root' value as a starting point(!!!!!)
	//Remember to initialize 'root' to some initial guess or to
	//point(0, 0) if you have no prior knowledge.
	//
	//poly - is an array of polynomial cooefs
	//
	//length = degree + 1, poly(1) is constant
	//	1              2				3
	//poly(1) x ^ 0 + poly(2) x ^ 1 + poly(3) x ^ 2 + ...
	//
	//degree - a degree of the polynomial
	//
	//root - input: guess for the value of a root
	//output : a root of the polynomial
	//iter - number of iterations performed(the number of polynomial
	//evaluations and stopping criterion evaluation)
	//
	//success - is false if routine reaches maximum number of iterations
	//
	//For a summary of the method go to :
	//http://en.wikipedia.org/wiki/Laguerre's_method
	//
	// static int FRAC_JUMP_EVERY = 10;
	// const int FRAC_JUMP_LEN = 10;
	// double FRAC_JUMPS[FRAC_JUMP_LEN] = { 0.64109297,
	// 	0.91577881, 0.25921289, 0.50487203,
	// 	0.08177045, 0.13653241, 0.306162,
	// 	0.37794326, 0.04618805, 0.75132137 }; // some random numbers

	c_T faq; //jump length

    f_T FRAC_ERR;
	if constexpr ( std::is_same_v< f_T, float > )
    {
        FRAC_ERR = ROUND_OFF_F;
    }
    else{FRAC_ERR = ROUND_OFF_D;} //Fractional Error for double precision

	c_T p, dp, d2p_half; //value of polynomial, 1st derivative, and 2nd derivative
	static int i, k;
	bool good_to_go;
	c_T denom, denom_sqrt, dx, newroot;
	f_T ek, absroot, abs2p;
	c_T fac_newton, fac_extra, F_half, c_one_nth;
	f_T one_nth, n_1_nth, two_n_div_n_1;
	c_T c_one(1, 0);
	c_T zero(0, 0);
	f_T stopping_crit2;

	//--------------------------------------------------------------------------------------------

	// //EXTREME FAILSAFE! not usually needed but kept here just to be on the safe side. Takes care of first coefficient being 0
	// if (false) {
	// 	if (degree < 0) {
	// 		printf("Error: cmplx_laguerre: degree<0");
	// 		return;
	// 	}
	// 	if (poly[degree] == c_T(0, 0)) {
	// 		if (degree == 0) return;
	// 		cmplx_laguerre(poly, degree - 1, root, iter, success);
	// 	}
	// 	if (degree <= 1) {
	// 		if (degree == 0) {
	// 			success = false; // we just checked if poly[0] is zero and it isnt
	// 			printf("Warning: cmplx_laguerre: degree = 0 and poly[0] does not equal zero, no roots");
	// 			return;
	// 		}
	// 		else {
	// 			*root = -poly[0] / poly[1];
	// 			return;
	// 		}
	// 	}
	// } // End of EXTREME failsafe

	good_to_go = false;
	one_nth = 1 / f_T(degree);			// 1/n
	n_1_nth = (degree - 1)*one_nth;		// (n-1)/n
	two_n_div_n_1 = 2 / n_1_nth;			// 2n/(n-1)
	c_one_nth.set(one_nth, 0);		// (1/n,0)
	for (i = 1; i <= MAXIT; i++) {
		ek = abs_c(poly[degree]); // Preparing stopping criterion
		absroot = abs_c(*root);
		// Calculate the values of polynomial and its first and second derivatives
		p = poly[degree];
		dp = zero;
		d2p_half = zero;
		for (k = degree - 1; k >= 0; k--) {
			d2p_half = dp + d2p_half*(*root);
			dp = p + dp * (*root);
			p = poly[k] + p*(*root); // b_k
									 //Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
									 //Communications of the ACM, Volume 10 Issue 10, Oct. 1967, p. 655
									 //ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
									 //Eq 8.
			ek = absroot*ek + abs_c(p);
		}
		iter += 1;

		abs2p = real(conj(p)*p);
		if (abs2p == 0) return;
		stopping_crit2 = FRAC_ERR*ek * FRAC_ERR*ek;
		if (abs2p < stopping_crit2) {
			//(simplified a little Eq. 10 of Adams 1967)
			//do additional iteration if we are less than 10x from stopping criterion
			if (abs2p < 0.01f*stopping_crit2) {
				return; // we are at a good place!
			}
			else {
				good_to_go = true;
			}
		}
		else {
			good_to_go = false;
		}

		// faq = 1.0;
		denom = zero;
		if (dp != zero) {
			fac_newton = p / dp;
			fac_extra = d2p_half / dp;
			F_half = fac_newton*fac_extra;
			denom_sqrt = sqrt(c_one - two_n_div_n_1*F_half);

			//NEXT LINE PROBABLY CAN BE COMMENTED OUT. Check if compiler outputs positive real
			if (real(denom_sqrt) >= 0) {
				denom = c_one_nth + n_1_nth*denom_sqrt;
			}
			else {
				denom = c_one_nth - n_1_nth*denom_sqrt;
			}
		}

		if (denom == zero) {
			// dx = (absroot + 1.0)*expcmplx(c_T(0, FRAC_JUMPS[i % FRAC_JUMP_LEN] * 2 * M_PI));
			// dx = f_T(abs_c(*root) + 1) * c_T(0.626f, 0.212f);        // some random () numbers
            dx.set( 0.626f, 0.212f );
            dx *= ( 1 + abs_c( * root ) );            
		}
		else {
			dx = fac_newton / denom;
		}


		newroot = *root - dx;
		if (newroot == *root) return; //nothing changes so return
		if (good_to_go) {
			*root = newroot;
			return;
		}
		if ((i & (FRAC_JUMP_EVERY-1)) == 0) { //decide whether to do a jump of modified length (to break cycles)
			// faq = FRAC_JUMPS[(i / FRAC_JUMP_EVERY - 1) % FRAC_JUMP_LEN];
			faq.set(0.108f,0.525f);  // some random () numbers
			newroot = *root - faq*dx; // do jump of semi-random length
		}
		*root = newroot;
	}
	success = false; // too many iterations here
	return;
}

template< class c_T >
__device__ void cmplx_laguerre2newton(c_T *poly, int degree, c_T *root, int &iter, bool &success, int starting_mode)
{
    using f_T = typename c_T::f_t;
    
	//Subroutine finds one root of a c_T polynomial using
	//Laguerre's method, Second-order General method and Newton's
	//method - depending on the value of function F, which is a 
	//combination of second derivative, first derivative and
	//value of polynomial [F=-(p"*p)/(p'p')].

	//Subroutine has 3 modes of operation. It starts with mode=2
	//which is the Laguerre's method, and continues until F
	//becames F<0.50, at which point, it switches to mode=1,
	//i.e., SG method (see paper). While in the first two
	//modes, routine calculates stopping criterion once per every
	//iteration. Switch to the last mode, Newton's method, (mode=0)
	//happens when becomes F<0.05. In this mode, routine calculates
	//stopping criterion only once, at the beginning, under an
	//assumption that we are already very close to the root.
	//If there are more than 10 iterations in Newton's mode,
	//it means that in fact we were far from the root, and
	//routine goes back to Laguerre's method (mode=2).

	//Uses 'root' value as a starting point (!!!!!)
	//Remember to initialize 'root' to some initial guess or to 
	//point (0,0) if you have no prior knowledge.

	//poly - is an array of polynomial cooefs
	//	0					1				2
	//	poly[0] x^0 + poly[1] x^1 + poly[2] x^2
	//degree - a degree of the polynomial
	//root - input: guess for the value of a root
	//		output: a root of the polynomial
	//iter - number of iterations performed (the number of polynomial
	//		 evaluations and stopping criterion evaluation)
	//success - is false if routine reaches maximum number of iterations
	//starting_mode - this should be by default = 2. However if you  
	//				  choose to start with SG method put 1 instead.
	//				  Zero will cause the routine to
	//				  start with Newton for first 10 iterations, and
	//				  then go back to mode 2.

	//For a summary of the method see the paper: Skowron & Gould (2012)

	// int FRAC_JUMP_EVERY = 10;
	// const int FRAC_JUMP_LEN = 10;
	// double FRAC_JUMPS[FRAC_JUMP_LEN] = { 0.64109297, 0.91577881, 0.25921289, 0.50487203, 0.08177045, 0.13653241, 0.306162, 0.37794326, 0.04618805, 0.75132137 }; //some random numbers

	c_T faq; //jump length
	f_T FRAC_ERR;
	if constexpr ( std::is_same_v< f_T, float > )
    {
        FRAC_ERR = ROUND_OFF_F;
    }
    else{FRAC_ERR = ROUND_OFF_D;} //Fractional Error for double precision

	c_T p; //value of polynomial
	c_T dp; //value of 1st derivative
	c_T d2p_half; //value of 2nd derivative
	int i, j, k;
	bool good_to_go;
	//c_T G, H, G2;
	c_T denom, denom_sqrt, dx, newroot;
	f_T ek, absroot, abs2p, abs2_F_half;
	c_T fac_netwon, fac_extra, F_half, c_one_nth;
	f_T one_nth, n_1_nth, two_n_div_n_1;
	int mode;
	c_T c_one(1, 0);
	c_T zero(0, 0);
	f_T stopping_crit2;

	iter = 0;
	success = true;
	stopping_crit2 = 0; //value not important, will be initialized anyway on the first loop

						//next if block is an EXTREME failsafe, not usually needed, and thus turned off in this version.

	// if (false) {//change false to true if you would like to use caution about having first coefficent == 0
	// 	if (degree < 0) {
	// 		printf("Error: cmplx_laguerre2newton: degree < 0");
	// 		return;
	// 	}
	// 	if (poly[degree] == zero) {
	// 		if (degree == 0) return;
	// 		cmplx_laguerre2newton(poly, degree, root, iter, success, starting_mode);
	// 		return;
	// 	}
	// 	if (degree <= 1) {
	// 		if (degree == 0) {//// we know from previous check that poly[0] not equal zero
	// 			success = false;
	// 			printf("Warning: cmplx_laguerre2newton: degree = 0 and poly[0] = 0, no roots");
	// 			return;
	// 		}
	// 		else {
	// 			*root = -poly[0] / poly[1];
	// 			return;
	// 		}
	// 	}
	// }
	// //end EXTREME failsafe

	j = 1;
	good_to_go = false;

	mode = starting_mode; // mode = 2 full laguerre, mode = 1 SG, mode = 0 newton

	while(true) { //infinite loop, just to be able to come back from newton, if more than 10 iteration there

			   ////////////
			   ///mode 2///
			   ////////////

		if (mode >= 2) {//Laguerre's method
			one_nth = 1 / f_T(degree); ///
			n_1_nth = (degree - 1) * one_nth; ////
			two_n_div_n_1 = 2 / n_1_nth;
			c_one_nth.set(one_nth, 0);

			for (i = 1; i <= MAXIT; i++) {
				// faq = 1.0;

				//prepare stoping criterion
				ek = abs_c(poly[degree]);			// error of polynomial calculation
				absroot = abs_c(*root);
				//calculate value of polynomial and its first two derivative
				p = poly[degree];
				dp = zero;
				d2p_half = zero;		// p''(z)/2
				for (k = degree; k >= 1; k--) {//Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
					d2p_half = dp + d2p_half * (*root);
					dp = p + dp * (*root);
					p = poly[k - 1] + p * (*root); // b_k
												   //Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
												   //Communications of the ACM, Volume 10 Issue 10, Oct. 1967, p. 655
												   //ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
												   //Eq 8.
					ek = absroot * ek + abs_c(p);
				}
				abs2p = real(conj(p) * p); // abs(p)
				iter = iter + 1;
				if (abs2p == 0) return;

				stopping_crit2 = FRAC_ERR * ek * FRAC_ERR * ek;		// 终止条件
				if (abs2p < stopping_crit2) {//(simplified a little Eq. 10 of Adams 1967)
											 //do additional iteration if we are less than 10x from stopping criterion
					if (abs2p < 0.01f*stopping_crit2) return; // ten times better than stopping criterion
															 //return immediately, because we are at very good place
					else {
						good_to_go = true; //do one iteration more
					}
				}
				else {
					good_to_go = false; //reset if we are outside the zone of the root
				}

				denom = zero;
				if (dp != zero) {
					fac_netwon = p / dp;
					fac_extra = d2p_half / dp;
					F_half = fac_netwon * fac_extra;

					abs2_F_half = real(conj(F_half) * F_half);
					if (abs2_F_half <= 0.0625f) {//F<0.50, F/2<0.25
												//go to SG method
						if (abs2_F_half <= 0.000625f) {//F<0.05, F/2<0.025
							mode = 0; //go to Newton's
						}
						else {
							mode = 1; //go to SG
						}
					}

					denom_sqrt = sqrt(c_one - two_n_div_n_1*F_half);	// sqrt part in delta_Lag

					//NEXT LINE PROBABLY CAN BE COMMENTED OUT 
					if (real(denom_sqrt) > 0) {
						//real part of a square root is positive for probably all compilers. You can \F9
						//test this on your compiler and if so, you can omit this check
						denom = c_one_nth + n_1_nth * denom_sqrt;
					}
					else {
						denom = c_one_nth - n_1_nth * denom_sqrt;
					}
				}
				if (denom == zero) {//test if demoninators are > 0.0 not to divide by zero
					// dx = (abs(*root) + 1.0) + expcmplx(c_T(0.0, FRAC_JUMPS[i% FRAC_JUMP_LEN] * 2 * M_PI)); //make some random jump
					// dx = f_T(abs_c(*root) + 1.0) * c_T(0.626f, 0.212f);        // some random () numbers
                    dx.set( 0.626f, 0.212f );
                    dx *= ( 1 + abs_c( * root ) );                    
				}
				else {
					dx = fac_netwon / denom;
				}
				newroot = *root - dx;
				if (newroot == *root) return; // nothing changes -> return
				if (good_to_go) {//this was jump already after stopping criterion was met
					*root = newroot;
					return;
				}
				if (mode != 2) {
					*root = newroot;
					j = i + 1; //remember iteration index
					break; //go to Newton's or SG
				}
				if ((i & (FRAC_JUMP_EVERY-1)) == 0) {//decide whether to do a jump of modified length (to break cycles)
					// faq = FRAC_JUMPS[((i / FRAC_JUMP_EVERY - 1) % FRAC_JUMP_LEN)];
					faq.set(0.108f,0.525f);  // some random () numbers
					newroot = *root - faq * dx; // do jump of some semi-random length (0 < faq < 1)
				}
				*root = newroot;
			} //do mode 2

			if (i >= MAXIT) {
				success = false;
				return;
			}
		}

		////////////
		///mode 1///
		////////////

		if (mode == 1) {//SECOND-ORDER GENERAL METHOD (SG)

			for (i = j; i <= MAXIT; i++) {
				// faq = 1.0;
				//calculate value of polynomial and its first two derivatives
				p = poly[degree];
				dp = zero;
				d2p_half = zero;
				if (((i - j) & 7) == 0) {
					//prepare stopping criterion
					ek = abs_c(poly[degree]);
					absroot = abs_c(*root);
					for (k = degree; k >= 1; k--) {//Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
						d2p_half = dp + d2p_half * (*root);
						dp = p + dp * (*root);
						p = poly[k - 1] + p * (*root); //b_k
													   //Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
													   //Communications of the ACM, Volume 10 Issue 10, Oct. 1967, p. 655
													   //ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
													   //Eq 8.
						ek = absroot * ek + abs_c(p);
					}
					stopping_crit2 = FRAC_ERR*ek * FRAC_ERR*ek;
				}
				else {
					for (k = degree; k >= 1; k--) {//Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
						d2p_half = dp + d2p_half * (*root);
						dp = p + dp * (*root);
						p = poly[k - 1] + p * (*root); //b_k
					}
				}
				abs2p = real(conj(p) * p); //abs(p)**2
				iter = iter + 1;
				if (abs2p == 0) return;

				if (abs2p < stopping_crit2) {//(simplified a little Eq. 10 of Adams 1967)
					if (dp == zero) return;
					//do additional iteration if we are less than 10x from stopping criterion
					if (abs2p < 0.01f*stopping_crit2) return; //ten times better than stopping criterion
															 //ten times better than stopping criterion
					else {
						good_to_go = true; //do one iteration more
					}
				}
				else {
					good_to_go = false; //reset if we are outside the zone of the root
				}
				if (dp == zero) {//test if denominators are > 0.0 not to divide by zero
					// dx = (abs(*root) + 1.0) * expcmplx(c_T(0.0, FRAC_JUMPS[i% FRAC_JUMP_LEN] * 2 * M_PI)); //make some random jump
					// dx = f_T(abs_c(*root) + 1) * c_T(0.626f, 0.212f);        // some random () numbers
                    dx.set( 0.626f, 0.212f );
                    dx *= ( 1 + abs_c( * root ) );                    
				}
				else {
					fac_netwon = p / dp;
					fac_extra = d2p_half / dp;
					F_half = fac_netwon * fac_extra;

					abs2_F_half = real(conj(F_half) * F_half);
					if (abs2_F_half <= 0.000625f) {//F<0.05, F/2<0.025
						mode = 0; //set Newton's, go there after jump
					}
					dx = fac_netwon * (c_one + F_half); //SG
				}
				newroot = *root - dx;
				if (newroot == *root) return; //nothing changes -> return
				if (good_to_go) {
					*root = newroot; //this was jump already after stopping criterion was met
					return;
				}
				if (mode != 1) {
					*root = newroot;
					j = i + 1; //remember iteration number
					break; //go to Newton's
				}
				if ((i &  (FRAC_JUMP_EVERY-1)) == 0) {// decide whether to do a jump of modified length (to break cycles)
					// faq = FRAC_JUMPS[(i / FRAC_JUMP_EVERY - 1) % FRAC_JUMP_LEN];
					faq.set(0.108f,0.525f);  // some random () numbers
					newroot = *root - faq * dx; //do jump of some semi random lenth (0 < faq < 1)		
				}
				*root = newroot;
			}
			if (i >= MAXIT) {
				success = false;
				return;
			}

		}
		//------------------------------------------------------------------------------- mode 0
		if (mode == 0) { // Newton's Method

			for (i = j; i <= j + 10; i++) { // Do only 10 iterations the most then go back to Laguerre
				// faq = 1.0;

				//calc polynomial and first two derivatives
				p = poly[degree];
				dp = zero;
				if (i == j) { // Calculating stopping criterion only at the beginning
					ek = abs_c(poly[degree]);
					absroot = abs_c(*root);
					for (k = degree; k >= 1; k--) {
						dp = p + dp*(*root);
						p = poly[k - 1] + p*(*root);
						ek = absroot*ek + abs_c(p);
					}
					stopping_crit2 = FRAC_ERR*ek * FRAC_ERR*ek;
				}
				else {
					for (k = degree; k >= 1; k--) {
						dp = p + dp*(*root);
						p = poly[k - 1] + p*(*root);
					}
				}
				abs2p = real(conj(p)*p);
				iter = iter + 1;
				if (abs2p == 0) return;

				if (abs2p < stopping_crit2) {
					if (dp == zero) return;
					// do additional iteration if we are less than 10x from stopping criterion
					if (abs2p < 0.01f*stopping_crit2) {
						return; // return immediately since we are at a good place
					}
					else {
						good_to_go = true; // do one more iteration
					}
				}
				else {
					good_to_go = false;
				}

				if (dp == zero) {
					// dx = (abs(*root) + 1.0)*expcmplx(c_T(0.0, 2 * M_PI*FRAC_JUMPS[i % FRAC_JUMP_LEN])); // make a random jump
					// dx = f_T(abs_c(*root) + 1) * c_T(0.626f, 0.212f);        // some random () numbers
                    dx.set( 0.626f, 0.212f );
                    dx *= ( 1 + abs_c( * root ) );                    
				}
				else {
					dx = p / dp;
				}

				newroot = *root - dx;
				if (newroot == *root) return;
				if (good_to_go) {
					*root = newroot;
					return;
				}
				*root = newroot;
			}
			if (iter >= MAXIT) {
				//too many iterations
				success = false;
				return;
			}
			mode = 2; //go back to Laguerre's. Happens when could not converge with 10 steps of Newton
		}

	}/// end of infinite loop
}


}; // namespace twinkle::root_finder
