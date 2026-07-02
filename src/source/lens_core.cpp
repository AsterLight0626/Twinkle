#include "source_base.h"
#include "../utils/solve.h"

namespace twinkle
{

using f_t = source_base_t::f_t;
using c_t = source_base_t::c_t;

template< class f_t, class g_T > __device__ __forceinline__
auto max_f( const f_t & a, const g_T & b )
{
    return ( a > b ? a : b );
}

__device__ f_t source_base_t::yield_lens_coef
( c_t * coef, const c_t & zeta, const f_t lens_s ) const
{
    // const auto m1 = 1 / ( q + 1 );
    // const auto m2 = q / ( q + 1 );

    c_t  temp[ 3 ] , zeta_c;
    
	auto displacement   ( lens_s * m1 );
	auto y  =  zeta - displacement;
	// auto yc = conj(      y );
	// zeta_c  = conj(   zeta );
	auto yc =    y.conj(  ) ;
	zeta_c  = zeta.conj(  ) ;

	temp[ 0 ] = m2 * lens_s;
	temp[ 1 ] = zeta_c + temp[ 0 ];
	temp[ 2 ] = ( 1 - lens_s * lens_s ) + lens_s * temp[ 1 ];

	coef[ 5 ] = - yc * temp[1];
	coef[ 4 ] = ( y.norm2(   ) - 1 - 2 * lens_s * yc )
              * temp[ 1 ] + temp[ 0 ];

    c_t ctmp( 0, 2 * y.im );
	coef[ 3 ]  = ( 2 * y - lens_s )
        * ( temp [ 2 ] * temp[ 1 ] + ctmp - temp[ 0 ] );
    ctmp.set( 0, 4 * y.im );
    coef[ 3 ] += ( temp[ 0 ] - y ) * ctmp ;

    ctmp.set( 0, 1 - lens_s * y.re );
	coef[ 2 ]  = temp[ 2 ].re
        * ( temp[ 2 ].re * temp[ 1 ].re  + y.im * ctmp )
        + lens_s * y.im * y.im * ( 2 + lens_s * temp[ 1 ] );
    ctmp.set( - y.re, 3 * y.im );
    coef[ 2 ] += temp [ 0 ]
        * ( 2 * y.norm2(   ) + lens_s * ctmp - 1 );

    ctmp.set( 0, 2 * lens_s * y.im );
	coef[ 1 ] = temp[ 0 ]
        * ( ( lens_s + 2 * y ) * temp[ 2 ] + lens_s * ( ctmp - m2 ) );
    
	coef[ 0 ] = temp[ 0 ] * temp[ 0 ] * y ;

	return displacement;    
}

__device__ bool source_base_t::solve_lens_eq
( img_t * img, c_t * coef, const bool use_init_roots) const
{
    return root_finder::cmplx_roots_gen
         ( img, coef, order - 1, true, use_init_roots );
}

__device__ bool source_base_t::solve_imgs
( img_t * img, const c_t & zeta, const bool use_init_roots, const f_t lens_s, c_t * roots_point) const
{
    c_t coeffs[order];
    bool fail;

    const auto & displacement = yield_lens_coef( coeffs , zeta, lens_s );

    fail = solve_lens_eq( img , coeffs , use_init_roots);

    if(roots_point != nullptr)
    {
        for(int j=0;j<order-1;j++)
        {
            roots_point[j] = img[j].position;
        }
    }

    for(int j=0;j<5;j++)
    {
        img[j].position.re += displacement;
    }
    return fail;
}





__host__ __device__ c_t source_base_t::f_zbar
( const c_t & z, const f_t lens_s ) const
{
    return -(m1 / (z.conj(  ) + lens_s*m2) + m2 / (z.conj(  ) - lens_s*m1));
}

// <<<<<, 1<2<3<4<5, smaller to bigger
__device__ void source_base_t::bubble_arg_sort
( int* index, f_t* arr, const int len ) const
{
	f_t temp;
	int temp_idx;
	int i, j;
	for (i = 0; i < len - 1; i++)
	{
		for (j = 0; j < len - 1 - i; j++)
		{
			if (arr[j] > arr[j + 1])
			{
				// swap(arr[j], arr[j + 1]);
				temp = arr[j];
				arr[j] = arr[j+1];
				arr[j+1] = temp;
				temp_idx = index[j];
				index[j] = index[j+1];
				index[j+1] = temp_idx;
			}		
		}
	} 
    return ;   
}

__device__ void source_base_t::is_physical_all
( img_t * imgs, const c_t & zeta, int& Nphys, const f_t lens_s ) const
{
    f_t errors[order-1];
    c_t temp;
    int rank[order-1];
    int count;
    bool jump_test;
    bool residual_test;

    // if constexpr (std::is_same_v< f_t, float >)
    // {
    //     dlmin = 1e-3;
    // }

	for(int i=0;i<order-1;i++)
	{
		temp = f_zbar(imgs[i].position, lens_s ) + imgs[i].position - zeta;	// temp = z + f(zbar) - zeta
		errors[i] = temp.abs(  );
		rank[i] = i;
	}
    
    bubble_arg_sort(rank,errors,order-1);

    if(Nphys==-1)		// 前面没做判断
	{
		for(int j=0;j<5;j++)
		{
			imgs[j].physical = true;		// 先全都是 true，查到false再改
		}

        if constexpr (std::is_same_v< f_t, float >)
        {
            jump_test = errors[3]*1e-3 > (errors[2] + 1e-9);
        }
        else
        {
            jump_test = errors[3]*1e-4 > (errors[2] + 1e-12);
        }
        
		// if(errors[3]*1e-4 > (errors[2] + 1e-12))		// same criterion as VBB, "1e-12" is for double
        if(jump_test)
		// 在很远离焦散线的地方，由于数值不稳定，有时残差会相差比较大（甚至达到1e-7），在这种安全区域可以用 jump 来判断，效果也比较好
		{
			imgs[rank[4]].physical = false;
			imgs[rank[3]].physical = false;
			Nphys = 3;
		}
		else
		{
			count = 0;
			for(int j=4;j>0;j--)
			{
                if constexpr (std::is_same_v< f_t, float >)
                {
                    residual_test = errors[j]>2e-5 && errors[j-1]>2e-5 && fabs(errors[j] / errors[j-1] - 1)<1e-6;
                }
                else
                {
                    residual_test = errors[j]>2e-12 && errors[j-1]>2e-12 && fabs(errors[j] / errors[j-1] - 1)<2e-15;
                }                
				if( residual_test )
				{
					imgs[rank[j]].physical = false;
					imgs[rank[j-1]].physical = false;
					break;
				}
				count++;
			}
			Nphys = ( ( count == 4 ) ? 5 : 3 );
		}
    }
    // else{printf("phy_tst not -1 at first!!!, = %d. ", Nphys);}


    return ;
}


__device__ f_t source_base_t::jacobian
( const c_t & z, const f_t lens_s ) const
{
    c_t d_f = m1/( ( z.conj(  ) + lens_s * m2 ) * ( z.conj(  ) + lens_s * m2 ) ) + m2 / ( (z.conj(  ) - lens_s * m1 ) * (z.conj(  ) - lens_s * m1) );
    return 1 - d_f.norm2(   );
}

__device__ f_t source_base_t::mu_qc
( const c_t & z, const f_t lens_s ) const
{
	c_t items[2];
	c_t f1,f2,f3;
	f_t J;
    // f_t temp0;
    // c_t& temp1 = items[0];
    // c_t& temp2 = items[1];
    // f_t& temp3 = items[1].im;


	items[0] = m1 / (( z + m2 * lens_s ) * ( z + m2 * lens_s ));
	items[1] = m2 / (( z - m1 * lens_s ) * ( z - m1 * lens_s ));
	f1 = items[0] + items[1];

	J = 1 - f1.norm2(  );

	items[0] *= -2 / (z + m2 * lens_s );
	items[1] *= -2 / (z - m1 * lens_s );
	f2 = items[0] + items[1];

	items[0] *= -3 / (z + m2 * lens_s );
	items[1] *= -3 / (z - m1 * lens_s );	
	f3 = items[0] + items[1];

	f1 = f1.conj(  );

    // c_t f1_2_f3(f1);
    // f1_2_f3 *= f1;
    // c_t f1_3_f2_2(f1);
    // f1_3_f2_2 *= f1_2_f3;
    // f1_3_f2_2 *= f2;
    // f1_3_f2_2 *= f2;

    // f1_2_f3 *= f3;


    // return 1/(J*J*J*J*J) * (fabs(-2*(3* f1_3_f2_2.re  -  (3-3*J + J*J / 2)* f2.norm2(  )  +  J * f1_2_f3.re)) + fabs(18* f1_3_f2_2.im));

	return 1/fabs(J*J*J*J*J) * (fabs(-2*(3* f1*f1*f1 * f2*f2  -  (3-3*J + J*J / 2)* f2.norm2(  )  +  J*f1*f1 * f3).re) + fabs(6*(3* f1*f1*f1 * f2*f2).im));
    // return 1/(J*J*J*J*J) * (fabs(-2*(3* (f1*f1*f1 * f2*f2).re  -  (3-3*J + J*J / 2)* f2.norm2(  )  +  J*(f1*f1 * f3).re)) + fabs(18* (f1*f1*f1 * f2*f2).im));

}

__device__ f_t source_base_t::ghost_tst
( const c_t & z_G, const c_t & z_hat, const f_t lens_s ) const
{
	c_t fz1,fz2,fz_hat1,J_hat;
	f_t J;		// J(z_G)
	c_t items[2];


	items[0] = m1 / (( z_G + m2 * lens_s ) * ( z_G + m2 * lens_s ));
	items[1] = m2 / (( z_G - m1 * lens_s ) * ( z_G - m1 * lens_s ));
	fz1 = items[0] + items[1];

	J = 1 - fz1.norm2(  );

	items[0] = -2 * items[0] / (z_G + m2 * lens_s );
	items[1] = -2 * items[1] / (z_G - m1 * lens_s );
	fz2 = items[0] + items[1];

	items[0] = m1 / (( z_hat + m2 * lens_s ) * ( z_hat + m2 * lens_s ));
	items[1] = m2 / (( z_hat - m1 * lens_s ) * ( z_hat - m1 * lens_s ));
	fz_hat1 = items[0] + items[1];


	J_hat = 1 - fz1*fz_hat1;

    c_t temp;	// the inverse of LHS of Eq(47) in VBB2

	temp =   J_hat * fz2.conj(  ) * fz1;
	temp = ( temp - temp.conj(  ) * fz_hat1 ) / (J*J_hat*J_hat);

	return temp.abs(  );
}

__device__ bool source_base_t::ghost_test_all
( const img_t * const imgs, const int & Nphys, const c_t & zeta, const f_t & Rho, const f_t lens_s ) const
{
    bool pass( true );
    if(Nphys==3)
    {
        f_t GT( 0 );
        for(int j=0;j<order-1;j++)
        {
            if(imgs[j].physical==false)
            {
                const c_t & z_hat = zeta.conj(  ) - f_zbar(imgs[j].position, lens_s );     // f_zbar: input z, return f(z_conj)
                GT = max_f(GT, f_t(ghost_tst(imgs[j].position,z_hat, lens_s )));
            }
        }
        pass = (GT* (Rho+1e-3)*2*C_G < 1);
    }
    return pass;
}

__device__ bool source_base_t::safe_distance_test
( const c_t & zeta, const f_t & Rho, const f_t lens_s ) const
{
    f_t safedist = 10;
    if(lens_q<0.01)     // q<0.01
    {
        const c_t zeta_pc( -1/lens_s + m1*lens_s, 0 );                          // slightly different from VBB2 Eq49, because the frame centre is different. params[1] is z_1
        const f_t & delta_pc2 = 9*lens_q /lens_s / lens_s;
        safedist = (zeta-zeta_pc).norm2(  ) - C_P*delta_pc2;
    }
    return (safedist > C_P * Rho*Rho);          
}

}; // namespace twinkle
