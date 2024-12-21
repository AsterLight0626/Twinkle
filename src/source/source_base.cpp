#include "source_base.h"
#include "../utils/solve.h"

namespace twinkle
{
////////////////////////////////////////////////////////////
//

__host__ source_base_t:: source_base_t(  )
{
    // n_th  = 32;
    // n_src = 64*2*2;
    // n_src = 1;
    n_src = 256;
    // n_src = 900;

    // Temporarily hard-coded
    n_point_max = 4096;
    n_cross_max =  n_th;
    return;
}

__host__ source_base_t::~source_base_t(  )
{
    return;
}

__host__ void source_base_t::init( device_t & f_dev )
{
    pool_margin.setup_mem( f_dev, n_src * n_point_max,
                             n_point_max  );
    pool_cross .setup_mem( f_dev, n_src * n_cross_max,
                             n_cross_max  );
    pool_center.setup_mem( f_dev, n_src );
    pool_mag   .setup_mem( f_dev, n_src );
    pool_extended   .setup_mem( f_dev, n_src );

    stream = f_dev.yield_stream(  );

    // s = ss;
    // q = qq;


    // RelTol    =   1e-4;
    // m1 = 1 / ( q + 1 );
    // m2 = q / ( q + 1 );
    // int Nx = 30;
    // int Ny = 30;

    // for(int y_idx=0;y_idx<Ny;y_idx++)
    // {
    //     for(int x_idx=0;x_idx<Nx;x_idx++)
    //     {
    //         auto & shape  = pool_center.dat_h[ x_idx + y_idx * Nx ];
    //         shape.rho = rho;
    //         shape.loc_centre.re = ((f_t(x_idx)+0.5) / Nx) * (xmax - xmin) + xmin;
    //         shape.loc_centre.im = ((f_t(y_idx)+0.5) / Ny) * (ymax - ymin) + ymin;            
    //     }
    // }


    // pool_center.cp_h2d( f_dev );
    // // std::cerr << "Finishing setting up initial condition"
    // //           << " for tests\n";
    return;
}

__host__ void source_base_t::set_params( device_t & f_dev, f_t ss, f_t qq, f_t rho, f_t xmax, f_t xmin, f_t ymax, f_t ymin )
{
    s = ss;
    q = qq;

    RelTol    =   1e-4;
    m1 = 1 / ( q + 1 );
    m2 = q / ( q + 1 );
    int Nx = 30;
    int Ny = 30;

    for(int y_idx=0;y_idx<Ny;y_idx++)
    {
        for(int x_idx=0;x_idx<Nx;x_idx++)
        {
            auto & shape  = pool_center.dat_h[ x_idx + y_idx * Nx ];
            shape.rho = rho;
            shape.loc_centre.re = ((f_t(x_idx)+0.5) / Nx) * (xmax - xmin) + xmin;
            shape.loc_centre.im = ((f_t(y_idx)+0.5) / Ny) * (ymax - ymin) + ymin;            
        }
    }

    pool_center.cp_h2d( f_dev );
    return;
}

__host__ void source_base_t::set_same( device_t & f_dev, f_t ss, f_t qq, f_t rho, f_t zeta_x, f_t zeta_y )
{
    s = ss;
    q = qq;

    RelTol    =   1e-4;
    m1 = 1 / ( q + 1 );
    m2 = q / ( q + 1 );

    for(int idx=0;idx<n_src;idx++)
    {
        auto & shape  = pool_center.dat_h[ idx ];
        shape.rho = rho;
        shape.loc_centre.re = zeta_x;
        shape.loc_centre.im = zeta_y;            
    }


    pool_center.cp_h2d( f_dev );
    return;
}


__host__ void source_base_t::free( device_t & f_dev )
{
    pool_mag   .free( f_dev );    
    pool_cross .free( f_dev );
    pool_center.free( f_dev );
    pool_margin.free( f_dev );
    pool_extended.free( f_dev );
    return;
}

////////////////////////////////////////////////////////////
//

using f_t = source_base_t::f_t;
using c_t = source_base_t::c_t;


template< class f_t, class g_T > __device__ __forceinline__
auto max_f( const f_t & a, const g_T & b )
{
    return ( a > b ? a : b );
}

__device__ f_t source_base_t::yield_lens_coef
( c_t * coef, const c_t & zeta ) const
{
    // const auto m1 = 1 / ( q + 1 );
    // const auto m2 = q / ( q + 1 );

    c_t  temp[ 3 ] , zeta_c;
    
	auto displacement   ( s * m1 );
	auto y  =  zeta - displacement;
	// auto yc = conj(      y );
	// zeta_c  = conj(   zeta );
	auto yc =    y.conj(  ) ;
	zeta_c  = zeta.conj(  ) ;

	temp[ 0 ] = m2 * s;
	temp[ 1 ] = zeta_c + temp[ 0 ];
	temp[ 2 ] = ( 1 - s * s ) + s * temp[ 1 ];

	coef[ 5 ] = - yc * temp[1];
	coef[ 4 ] = ( y.norm2(   ) - 1 - 2 * s * yc )
              * temp[ 1 ] + temp[ 0 ];

    c_t ctmp( 0, 2 * y.im );
	coef[ 3 ]  = ( 2 * y - s )
        * ( temp [ 2 ] * temp[ 1 ] + ctmp - temp[ 0 ] );
    ctmp.set( 0, 4 * y.im );
    coef[ 3 ] += ( temp[ 0 ] - y ) * ctmp ;

    ctmp.set( 0, 1 - s * y.re );
	coef[ 2 ]  = temp[ 2 ].re
        * ( temp[ 2 ].re * temp[ 1 ].re  + y.im * ctmp )
        + s * y.im * y.im * ( 2 + s * temp[ 1 ] );
    ctmp.set( - y.re, 3 * y.im );
    coef[ 2 ] += temp [ 0 ]
        * ( 2 * y.norm2(   ) + s * ctmp - 1 );

    ctmp.set( 0, 2 * s * y.im );
	coef[ 1 ] = temp[ 0 ]
        * ( ( s + 2 * y ) * temp[ 2 ] + s * ( ctmp - m2 ) );
    
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
( img_t * img, const c_t & zeta, const bool use_init_roots, c_t * roots_point) const
{
    c_t coeffs[order];
    bool fail;

    const auto & displacement = yield_lens_coef( coeffs , zeta );

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




__device__ c_t source_base_t::f_zbar
( const c_t & z ) const
{
    return -(m1 / (z.conj(  ) + s*m2) + m2 / (z.conj(  ) - s*m1));
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
( img_t * imgs, const c_t & zeta, int& Nphys ) const
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
		temp = f_zbar(imgs[i].position) + imgs[i].position - zeta;	// temp = z + f(zbar) - zeta
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
( const c_t & z ) const
{
    c_t d_f = m1/( ( z.conj(  ) + s*m2 ) * ( z.conj(  ) + s*m2 ) ) + m2 / ( (z.conj(  ) - s*m1) * (z.conj(  ) - s*m1) );
    return 1 - d_f.norm2(   );
}

__device__ f_t source_base_t::mu_qc
( const c_t & z ) const
{
	c_t items[2];
	c_t f1,f2,f3;
	f_t J;
    // f_t temp0;
    // c_t& temp1 = items[0];
    // c_t& temp2 = items[1];
    // f_t& temp3 = items[1].im;


	items[0] = m1 / (( z + m2 * s) * ( z + m2 * s));
	items[1] = m2 / (( z - m1 * s) * ( z - m1 * s));
	f1 = items[0] + items[1];

	J = 1 - f1.norm2(  );

	items[0] *= -2 / (z + m2 * s);
	items[1] *= -2 / (z - m1 * s);
	f2 = items[0] + items[1];

	items[0] *= -3 / (z + m2 * s);
	items[1] *= -3 / (z - m1 * s);	
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

	return 1/(J*J*J*J*J) * (fabs(-2*(3* f1*f1*f1 * f2*f2  -  (3-3*J + J*J / 2)* f2.norm2(  )  +  J*f1*f1 * f3).re) + fabs(6*(3* f1*f1*f1 * f2*f2).im));
    // return 1/(J*J*J*J*J) * (fabs(-2*(3* (f1*f1*f1 * f2*f2).re  -  (3-3*J + J*J / 2)* f2.norm2(  )  +  J*(f1*f1 * f3).re)) + fabs(18* (f1*f1*f1 * f2*f2).im));

}

__device__ f_t source_base_t::ghost_tst
( const c_t & z_G, const c_t & z_hat ) const
{
	c_t fz1,fz2,fz_hat1,J_hat;
	f_t J;		// J(z_G)
	c_t items[2];


	items[0] = m1 / (( z_G + m2 * s) * ( z_G + m2 * s));
	items[1] = m2 / (( z_G - m1 * s) * ( z_G - m1 * s));
	fz1 = items[0] + items[1];

	J = 1 - fz1.norm2(  );

	items[0] = -2 * items[0] / (z_G + m2 * s);
	items[1] = -2 * items[1] / (z_G - m1 * s);
	fz2 = items[0] + items[1];

	items[0] = m1 / (( z_hat + m2 * s) * ( z_hat + m2 * s));
	items[1] = m2 / (( z_hat - m1 * s) * ( z_hat - m1 * s));
	fz_hat1 = items[0] + items[1];


	J_hat = 1 - fz1*fz_hat1;

    c_t temp;	// the inverse of LHS of Eq(47) in VBB2

	temp =   J_hat * fz2.conj(  ) * fz1;
	temp = ( temp - temp.conj(  ) * fz_hat1 ) / (J*J_hat*J_hat);

	return temp.abs(  );
}

__device__ bool source_base_t::ghost_test_all
( const img_t * const imgs, const int & Nphys, const c_t & zeta, const f_t & Rho) const
{
    bool pass( true );
    if(Nphys==3)
    {
        f_t GT( 0 );
        for(int j=0;j<order-1;j++)
        {
            if(imgs[j].physical==false)
            {
                const c_t & z_hat = zeta.conj(  ) - f_zbar(imgs[j].position);     // f_zbar: input z, return f(z_conj)
                GT = max_f(GT, f_t(ghost_tst(imgs[j].position,z_hat)));
            }
        }
        pass = (GT* (Rho+1e-3)*2*C_G < 1);
    }
    return pass;
}

__device__ bool source_base_t::safe_distance_test
( const c_t & zeta, const f_t & Rho) const
{
    f_t safedist = 10;
    if(q<0.01)     // q<0.01
    {
        const c_t zeta_pc( -1/s + m1*s, 0 );                          // slightly different from VBB2 Eq49, because the frame centre is different. params[1] is z_1
        const f_t & delta_pc2 = 9*q /s / s;
        safedist = (zeta-zeta_pc).norm2(  ) - C_P*delta_pc2;
    }
    return (safedist > C_P * Rho*Rho);          
}



__device__ void source_base_t::margin_Q
( c_t & marginloc, const src_shape_t<f_t>& shape, const float& q ) const
{
    f_t theta = q*2*PI;
    marginloc.set( shape.loc_centre.re, shape.loc_centre.im );

    marginloc.re += cos(theta) * shape.rho;
    marginloc.im += sin(theta) * shape.rho;
    return;
}

__device__ void source_base_t::set_skip_info
( src_pt_t < f_t > & src ) const
{
    for(int j=0;j<order-1;j++)
    {
        src.deltaS[j] = 0;
        src.deltaS_Err[j] = 0;
        // src.deltaS[j] = 0;
        // src.deltaS_Err[j] = 0; 
        src.images[j].position = 0;
        src.images[j].physical = false;
        src.next_idx[j] = -10;
        src.next_j[j] = -10;
        src.wedge[j] = 0;
        src.dz[j] = 0;
    }     
    src.error_interval = 0;
    src.area_interval = 0;
    return;
}


__device__ void source_base_t::dz_wedge
( c_t & dz, f_t & wedge, f_t & jacobian, const img_t & img, const c_t & zeta, const c_t & loc_centre, const f_t rho) const
{
	c_t pzeta_pconjz,ppzeta_ppconjz;
	c_t items[2];
    c_t zc;
    zc.set(img.position.re, -img.position.im);

	items[0] = m1 / (( zc + m2 * s) * ( zc + m2 * s));
	items[1] = m2 / (( zc - m1 * s) * ( zc - m1 * s));
	pzeta_pconjz = items[0] + items[1];

	jacobian = 1 - pzeta_pconjz.norm2(  );

	items[0] = -2 * items[0] / (zc + m2 * s);
	items[1] = -2 * items[1] / (zc - m1 * s);
	ppzeta_ppconjz = items[0] + items[1];

    c_t dzeta;
    dzeta.set(0,1);
    dzeta *= (zeta - loc_centre);
    dz = (dzeta - pzeta_pconjz * ( dzeta.conj(  )) ) / jacobian;		// equation20 in VBB

    wedge = (rho*rho + ( dz * dz * dzeta * (ppzeta_ppconjz.conj(  )) ).im) / jacobian;

    return;
}

__device__ bool source_base_t::polish_pp
( img_t * imgs, const f_t * const jacobians) const
{
	int neat_parity=0;
	bool fail=false;
	
	for(int j=0;j<order-1;j++)
	{
		if(imgs[j].physical)
		{
			neat_parity += (imgs[j].parity) ? -1 : 1;
		}
	}
	if(neat_parity != -1)
	{
		int order_p[order-1];
		f_t temp[order-1];
		for(int j=0;j<order-1;j++)
		{
			order_p[j] = j;
			temp[j] = 1/jacobians[j];
		}
		bubble_arg_sort(order_p,temp,order-1);
		// order_p[i]是第i小的Mag（包括正负号）对应的指标


		if(abs(-1-neat_parity)>3)		// 0，2，-2是允许的，更多不行，例如4，-4
		{
			fail=true;
		}
		else
		{
			// 如果正手征多了，说明有的负手征跨过散焦线变成了正手征，此时放大率最大的嫌疑最大
			if(neat_parity > -1)
			{
				for(int jj=order-1-1;jj>=0;jj--)
				{
					if(imgs[order_p[jj]].physical)
					{
						imgs[order_p[jj]].parity = true;
						break;
					}
				}
			}
			// 如果负手征多了，说明有的正手征跨过散焦线变成了负手征，此时放大率最小的嫌疑最大
			else{
				for(int jj=0;jj<=order-1-1;jj++)
				{
					if(imgs[order_p[jj]].physical)
					{
						imgs[order_p[jj]].parity = false;
						break;
					}					
				}
			}			
		}
	}
	return fail;    
}

__device__ bool source_base_t::polish_pp2
( img_t * imgs, const f_t * const jacobians, const c_t & zeta) const
{
	int LCR[3];
	f_t LCR_re[3];

	int fifth;
	int same_sign = -1;
	f_t Jmax;
	Jmax=0;
	f_t tempf;

	int count=0;

	for(int j=0;j<5;j++)
	{
		if((imgs[j].position.im>0) == (zeta.im>0))
		{
			same_sign = j;			// only one candidate
			continue;
		}
		tempf = fabs(jacobians[j]);
		if(tempf > Jmax)
		{
			Jmax = tempf;
			fifth = j;
		}
	}
	for(int j=0;j<5;j++)
	{
		if(j!=same_sign && j!=fifth)
		{
			LCR[count] = j;
			LCR_re[count] = imgs[j].position.re;
			count++;
		}
	}

	bubble_arg_sort(LCR,LCR_re,3);

	if(imgs[LCR[0]].parity != 1)
	{
		imgs[LCR[0]].parity = 1;
	}
	if(imgs[LCR[1]].parity != 0)
	{
		imgs[LCR[1]].parity = 0;
	}
	if(imgs[LCR[2]].parity != 1)
	{
		imgs[LCR[2]].parity = 1;
	}
    return false;
} 


__device__ bool source_base_t::check_parity
( img_t * imgs, const f_t * const jacobians, const c_t & zeta, const int Nphys) const
{
    bool ambiguous = false;
    for(int j=0;j<5;j++)
    {
        if(fabsf(jacobians[j])<1e-5f){ambiguous = true;}
    }
    if(Nphys==5 && ambiguous){return polish_pp2( imgs, jacobians, zeta);}
    else{return polish_pp( imgs, jacobians );}
}

__device__ bool source_base_t::margin_solve_cal                // return skip
( img_t * imgs, c_t * dz, f_t * wedge, int & phys_tst, int * temp_j_from_out_j, const c_t & zeta, const c_t & loc_centre, const f_t Rho) const
{
    bool skip = solve_imgs(imgs, zeta, true);
    if(skip){ return true; }

    f_t temp_jacobi[ 5 ];
    for(int j=0;j<order-1;j++)
    {   
        dz_wedge( dz[j], wedge[j], temp_jacobi[j], imgs[j], zeta, loc_centre, Rho);
        imgs[j].parity = ( (temp_jacobi[j]>0) ? 0 : 1 );
    }   
    is_physical_all( imgs, zeta, phys_tst); 
    skip = check_parity( imgs, temp_jacobi, zeta, phys_tst);

    int temp_j_from_order_j[5];       // temp_j = temp_j_from_order_j[order_j]
    for(int order_i=0;order_i<5;order_i++)
    {
        temp_j_from_order_j[order_i] = order_i;
    }
    bubble_arg_sort( temp_j_from_order_j, temp_jacobi, 5);

    int out_j=0;
    int order_j_from_out_j[5];        // out_j = order_j_from_out_j[order_j]
    for(int order_j=0;order_j<5;order_j++)
    {
        int temp_j = temp_j_from_order_j[order_j];
        if(imgs[temp_j].physical && (imgs[temp_j].parity==1))
        {
            order_j_from_out_j[out_j] = order_j;
            out_j += 1;
        }
        if((out_j + int(phys_tst==3))==3){ break; }    // N=5, 3 break; N=3, 2 break
    }
    out_j = 4;
    for(int order_j=4;order_j>=0;order_j--)
    {
        int temp_j = temp_j_from_order_j[order_j];
        if(imgs[temp_j].physical && (imgs[temp_j].parity==0))
        {
            order_j_from_out_j[out_j] = order_j;
            out_j -= 1;
        }
        if((out_j - int(phys_tst==3))==2){ break; }    // N=5, 2 break, N=3, 3 break
    }
    if(phys_tst==3)
    {
        out_j = 2;
        for(int order_j=0;order_j<5;order_j++)
        {
            int temp_j = temp_j_from_order_j[order_j];
            if(!imgs[temp_j].physical)
            {
                order_j_from_out_j[out_j] = order_j;
                out_j += 1;
            }
            if((out_j)==4){ break; }
        }            
    }

    for(int out_j=0;out_j<5;out_j++)
    {
        temp_j_from_out_j[out_j] = temp_j_from_order_j[order_j_from_out_j[out_j]];
    }

    return skip;
}



__device__ int source_base_t::prev_src_idx_g
( const int idx_here, const int srcidx ) const
{
	bool skip;
	int idx_prev;
	idx_prev = pool_margin[ srcidx * n_point_max + idx_here].prev_src_idx;
	skip = pool_margin[ srcidx * n_point_max + idx_prev ].skip;
	while(skip)
	{
		idx_prev = pool_margin[ srcidx * n_point_max + idx_prev ].prev_src_idx;
		skip = pool_margin[ srcidx * n_point_max + idx_prev ].skip;
	}
	return idx_prev;
} 

__device__ int source_base_t::next_src_idx_g
( const int idx_here, const int srcidx ) const
{
	bool skip;
	int idx_next;
	idx_next = pool_margin[ srcidx * n_point_max + idx_here ].next_src_idx;
	skip = pool_margin[ srcidx * n_point_max + idx_next ].skip;
	while(skip)
	{
		idx_next = pool_margin[ srcidx * n_point_max + idx_next ].next_src_idx;
		skip = pool_margin[ srcidx * n_point_max + idx_next ].skip;
	}
	return idx_next;
}


__device__ int source_base_t::prev_src_idx_local
( const int idx_here, local_info_t<f_t>& local_info, bool * changed ) const
{
	bool skip;
	int idx_prev;
    int batchidx = local_info.batchidx;

    idx_prev = local_info.new_pts[idx_here - batchidx*n_th].prev_src_idx;
    if(idx_prev >= batchidx*n_th)       // 如果是新点，验证是否skip
    {
        skip = local_info.new_pts[idx_prev - batchidx*n_th].skip;
        if(skip){*changed = true;}
        while(skip)
        {
            idx_prev = local_info.new_pts[idx_prev - batchidx*n_th].prev_src_idx;
            if(idx_prev < batchidx*n_th)
            {
                break;
            }
            skip = local_info.new_pts[idx_prev - batchidx*n_th].skip;
        }
    }
    // else{};  // 如果是旧点，一定不会skip，因为上一个循环改好了
    return idx_prev;
}

__device__ int source_base_t::next_src_idx_local
( const int idx_here, local_info_t<f_t>& local_info, bool * changed ) const
{
	bool skip;
	int idx_next;
    int batchidx = local_info.batchidx;

    idx_next = local_info.new_pts[idx_here - batchidx*n_th].next_src_idx;
    if(idx_next >= batchidx*n_th)
    {
        skip = local_info.new_pts[idx_next - batchidx*n_th].skip;
        if(skip){*changed = true;}
        while(skip)
        {
            idx_next = local_info.new_pts[idx_next - batchidx*n_th].next_src_idx;
            if(idx_next < batchidx*n_th)
            {
                break;          // 如果找到旧点了，它一定没有skip，就是它
            }
            skip = local_info.new_pts[idx_next - batchidx*n_th].skip;
        }
    } // 逻辑见 prev 版本
    return idx_next;
}

__device__ int source_base_t::prev_src_idx_local_g
( const int idx_here, local_info_t<f_t>& local_info, const int srcidx ) const
{
	bool skip;
	int idx_prev;
    int batchidx = local_info.batchidx;
    bool global = false;

    idx_prev = local_info.new_pts[idx_here - batchidx*n_th].prev_src_idx;
    if(idx_prev < batchidx*n_th)
    {    
        skip = local_info.new_pts[idx_here - batchidx*n_th].skip;
        while(skip)
        {
            idx_prev = local_info.new_pts[idx_prev - batchidx*n_th].prev_src_idx;
            if(idx_prev < batchidx*n_th)
            {
                global = true;
                break;
            }
            skip = local_info.new_pts[idx_prev - batchidx*n_th].skip;
        }
    }
    if(global)
    {
        skip = pool_margin[ srcidx * n_point_max + idx_prev ].skip;
        while(skip)
        {
            idx_prev = pool_margin[ srcidx * n_point_max + idx_prev ].prev_src_idx;
            skip = pool_margin[ srcidx * n_point_max + idx_prev ].skip;
        }
    }

    return idx_prev;
}

__device__ int source_base_t::next_src_idx_local_g
( const int idx_here, local_info_t<f_t>& local_info, const int srcidx ) const
{
	bool skip;
	int idx_next;
    int batchidx = local_info.batchidx;
    bool global = false;

    idx_next = local_info.new_pts[idx_here - batchidx*n_th].next_src_idx;
    if(idx_next < batchidx*n_th)
    {
        skip = local_info.new_pts[idx_here - batchidx*n_th].skip;
        while(skip)
        {
            idx_next = local_info.new_pts[idx_next - batchidx*n_th].next_src_idx;
            if(idx_next < batchidx*n_th)
            {
                global = true;
                break;
            }
            skip = local_info.new_pts[idx_next - batchidx*n_th].skip;
        }        
    }
    if(global)
    {
        skip = pool_margin[ srcidx * n_point_max + idx_next ].skip;
        while(skip)
        {
            idx_next = pool_margin[ srcidx * n_point_max + idx_next ].next_src_idx;
            skip = pool_margin[ srcidx * n_point_max + idx_next ].skip;
        }        
    }
    
    return idx_next;
}



__device__ void source_base_t::connect_order
( const c_t* posi_here, const c_t* posi_other, int len_here, int len_other, int* order) const
{
	f_t Distance[6];
	int jbest=0;
	f_t temp;
	c_t temp_c;
	if((len_here==2))
	{
		if(len_other==2)
		{
			// f_t Distance[2];	// 01,10
			Distance[0] = (posi_here[0]-posi_other[0]).abs(  ) + (posi_here[1]-posi_other[1]).abs(  );
			Distance[1] = (posi_here[0]-posi_other[1]).abs(  ) + (posi_here[1]-posi_other[0]).abs(  );

			if(Distance[0] < Distance[1])
			{
				order[0] = 0;
				order[1] = 1;
			}
			else
			{
				order[0] = 1;
				order[1] = 0;
			}			
		}
		else	// len_other==3
		{
			// f_t Distance[6];	// 01,02,12,10,20,21
			Distance[0] = (posi_here[0]-posi_other[0]).abs(  ) + (posi_here[1]-posi_other[1]).abs(  );
			Distance[1] = (posi_here[0]-posi_other[0]).abs(  ) + (posi_here[1]-posi_other[2]).abs(  );
			Distance[2] = (posi_here[0]-posi_other[1]).abs(  ) + (posi_here[1]-posi_other[2]).abs(  );
			Distance[3] = (posi_here[0]-posi_other[1]).abs(  ) + (posi_here[1]-posi_other[0]).abs(  );
			Distance[4] = (posi_here[0]-posi_other[2]).abs(  ) + (posi_here[1]-posi_other[0]).abs(  );
			Distance[5] = (posi_here[0]-posi_other[2]).abs(  ) + (posi_here[1]-posi_other[1]).abs(  );

			temp = Distance[0];
			for(int j=1;j<6;j++)
			{
				if(Distance[j] < temp)
				{
					temp = Distance[j];
					jbest = j;
				}
			}
			order[0] = jbest/2;		// 0/2=0,1/2=0, 2/2=1.3/2=1, 4/2=2,5/2=2
			order[1] = ((jbest+3)/2) % 3;	// 1,2,2,0,0,1
		}
	}
	else		// len_here==3
	{
		if(len_other==2)
		{
			// f_t Distance[6];	// 01X，0X1, 1X0, 10X, X01, X10
			Distance[0] = (posi_here[0]-posi_other[0]).abs(  ) + (posi_here[1]-posi_other[1]).abs(  );
			Distance[1] = (posi_here[0]-posi_other[0]).abs(  ) + (posi_here[2]-posi_other[1]).abs(  );
			Distance[2] = (posi_here[0]-posi_other[1]).abs(  ) + (posi_here[2]-posi_other[0]).abs(  );
			Distance[3] = (posi_here[0]-posi_other[1]).abs(  ) + (posi_here[1]-posi_other[0]).abs(  );
			Distance[4] = (posi_here[1]-posi_other[0]).abs(  ) + (posi_here[2]-posi_other[1]).abs(  );
			Distance[5] = (posi_here[1]-posi_other[1]).abs(  ) + (posi_here[2]-posi_other[0]).abs(  );

			temp = Distance[0];
			for(int j=1;j<6;j++)
			{
				if(Distance[j] < temp)
				{
					temp = Distance[j];
					jbest = j;
				}
			}
			order[0] = jbest/2;		// 0/2=0,1/2=0, 2/2=1.3/2=1, 4/2=2,5/2=2
			order[1] = ((jbest+3)/2) % 3;	// 1,2,2,0,0,1
			order[2] = 2 - (jbest%3);	// 2,1,0,2,1,0

			for(int j=0;j<6;j++)
			{
				if(order[j]==2){order[j]=-1;}	// 2为虚解，标记为未连接
			}
		}
		else	// len_other==3
		{
			// f_t Distance[6];	// 012, 021, 120, 102, 201, 210
			Distance[0] = (posi_here[0]-posi_other[0]).abs(  ) + (posi_here[1]-posi_other[1]).abs(  ) + (posi_here[2]-posi_other[2]).abs(  );
			Distance[1] = (posi_here[0]-posi_other[0]).abs(  ) + (posi_here[1]-posi_other[2]).abs(  ) + (posi_here[2]-posi_other[1]).abs(  );
			Distance[2] = (posi_here[0]-posi_other[1]).abs(  ) + (posi_here[1]-posi_other[2]).abs(  ) + (posi_here[2]-posi_other[0]).abs(  );
			Distance[3] = (posi_here[0]-posi_other[1]).abs(  ) + (posi_here[1]-posi_other[0]).abs(  ) + (posi_here[2]-posi_other[2]).abs(  );
			Distance[4] = (posi_here[0]-posi_other[2]).abs(  ) + (posi_here[1]-posi_other[0]).abs(  ) + (posi_here[2]-posi_other[1]).abs(  );
			Distance[5] = (posi_here[0]-posi_other[2]).abs(  ) + (posi_here[1]-posi_other[1]).abs(  ) + (posi_here[2]-posi_other[0]).abs(  );

			temp = Distance[0];
			for(int j=1;j<6;j++)
			{
				if(Distance[j] < temp)
				{
					temp = Distance[j];
					jbest = j;
				}
			}
			order[0] = jbest/2;		// 0/2=0,1/2=0, 2/2=1,3/2=1, 4/2=2,5/2=2
			order[1] = ((jbest+3)/2) % 3;	// 1,2,2,0,0,1
			order[2] = 2 - (jbest%3);	// 2,1,0,2,1,0
		}
	}
}

__device__ int source_base_t::positive_connect_j34     // return j_positive_candidate
( int* next_idx_out, int* next_j_out, const int Nphys_here, const int Nphys_next,
    const c_t * posi_here_positive, const c_t * posi_next, const int here_idx, const int next_idx) const
{
    int j_positive_candidate = -1;
    if(Nphys_here==3)
    {
        next_idx_out[3] = -13;
        next_j_out[3] = -13;  
    }
    if(Nphys_here==Nphys_next)
    {
        if(Nphys_here==3)
        {
            // output
            next_idx_out[4] = next_idx;
            next_j_out[4] = 4;
        }
        else    // Nphys_here == 5, 两个比较大小
        {
            int order_e[2];
            connect_order(posi_here_positive,posi_next,2,2,order_e);
            // output
            // j_here_positive = j_next = {4,3}
            next_idx_out[4] = next_idx;
            next_idx_out[3] = next_idx;
            next_j_out[4] = (order_e[0]==0 ? 4 : 3);
            next_j_out[3] = (order_e[1]==1 ? 3 : 4);
        }
    }
    else    // Nphys_here!=Nphys_next
    {
        if(Nphys_here > Nphys_next)     // here: 2, next: 1
        {
            if((posi_here_positive[0]-posi_next[0]).norm2(  ) < (posi_here_positive[1]-posi_next[0]).norm2(  )) // 连近的，现在00是近的
            {
                // output
                next_idx_out[4] = next_idx;
                next_j_out[4] = 4;
                next_idx_out[3] = here_idx;        // connect to here
                j_positive_candidate = 3;       // j3 还没连
            }
            else    // 10是近的
            {
                // output
                next_idx_out[3] = next_idx;
                next_j_out[3] = 4;
                next_idx_out[4] = here_idx;        // connect to here                            
                j_positive_candidate = 4;       // j4 还没连
            }
        }
        else            // here: 1, next: 2
        {
            if((posi_here_positive[0]-posi_next[0]).norm2(  ) < (posi_here_positive[0]-posi_next[1]).norm2(  )) // 连近的，现在00是近的
            {
                // output
                // pt_here.next_idx[4] = next_idx;
                // pt_here.next_j[4] = 4;
                next_idx_out[4] = next_idx;
                next_j_out[4] = 4;
            }            
            else    // 01是近的
            {
                // output
                // pt_here.next_idx[4] = next_idx;
                // pt_here.next_j[4] = 3;   
                next_idx_out[4] = next_idx;
                next_j_out[4] = 3;   
            }            
        }
    }
    return j_positive_candidate;
}

__device__ int source_base_t::negative_connect_j012    // return j_negative_candidate
( int* next_idx_out, int* next_j_out, const int Nphys_here, const int Nphys_prev,
    const c_t * posi_here_negative, const c_t * posi_prev, const int here_idx, const int prev_idx) const
{
    int order_o[3];
    int j_negative_candidate = -1;

    if(Nphys_here==3)
    {
        next_idx_out[2] = -14;
        next_j_out[2] = -14;  
    }
    if(Nphys_here==Nphys_prev)
    {
        if(Nphys_here==3)
        {
            connect_order(posi_here_negative, posi_prev,2,2,order_o);
            // output
            for(int ii=0;ii<2;ii++)
            {
                next_idx_out[ii] = prev_idx;
                next_j_out[ii] = order_o[ii];
            }           
        }
        else
        {
            connect_order(posi_here_negative, posi_prev,3,3,order_o);
            for(int ii=0;ii<3;ii++)
            {
                next_idx_out[ii] = prev_idx;
                next_j_out[ii] = order_o[ii];
            }                        
        }
    }
    else
    {
        if(Nphys_here > Nphys_prev) // here: 3, prev: 2 
        {
            connect_order(posi_here_negative,posi_prev,3,2,order_o);      // order is 0,1 or -1 (cross candidate)
            // output
            for(int ii=0;ii<3;ii++)
            {
                if(order_o[ii]>=0)  // real
                {
                    next_idx_out[ii] = prev_idx;
                    next_j_out[ii] = order_o[ii];                                
                }
                else    // cross candidate
                {
                    next_idx_out[ii] = here_idx;    // connect to here
                    j_negative_candidate = ii;
                }
            }                          
        }
        else        // here: 2, prev: 3
        {

            connect_order(posi_here_negative,posi_prev,2,3,order_o);

            // output
            for(int ii=0;ii<2;ii++)
            {
                next_idx_out[ii] = prev_idx;
                next_j_out[ii] = order_o[ii];
            }                        
        }
    }
    return j_negative_candidate;
}

__device__ bool source_base_t::cross_j34              // 可能修改 012
( int* next_j_out, const int j_positive_candidate,
const c_t * posi_here_negative, const c_t * posi_next_negative, const int here_idx) const
{
    // int temp1;
    int order_c[3];     // c means cross caustic
    bool fail = false;

    connect_order(posi_here_negative,posi_next_negative,3,2,order_c);
    for(int ii=0;ii<3;ii++)
    {
        if(order_c[ii]==-1)     // 连接会有一个落单的，即离每个nextnegative都很远，这个点被positive cross连接
        {
            next_j_out[j_positive_candidate] = ii;
        }
    }
    return fail;
}

__device__ bool source_base_t::cross_j012               // 可能修改34
( int* next_j_out, const int j_negative_candidate,
const c_t * posi_here_positive, const c_t * posi_prev_positive, const int here_idx) const
{
    // int temp0;
    bool fail = false;

    // posi_here_positive[2] connect to posi_prev_positive[1]
    if((posi_here_positive[0]-posi_prev_positive[0]).norm2(  ) < (posi_here_positive[1]-posi_prev_positive[0]).norm2(  ))
    {
        next_j_out[j_negative_candidate] = 3;
    }
    else    // posi_here_positive[0]落单
    {
        next_j_out[j_negative_candidate] = 4;
    }
    return fail;
}


// j_here_negative = j_prev = {0,1,2} or {0,1,-1}
__device__ bool source_base_t::connect_prev_j012
( int* next_idx_out, int* next_j_out,
const img_t * imgs_here, const img_t * imgs_prev, const int here_idx, const int prev_idx, 
const int Nphys_here, const int Nphys_prev) const
{
    c_t posi_here_negative[3];
    c_t posi_prev[3];
    int j_negative_candidate;
    bool fail = false;

    posi_here_negative[0] = imgs_here[0].position;
    posi_here_negative[1] = imgs_here[1].position;
    if(Nphys_here == 5)
    {
        posi_here_negative[2] = imgs_here[2].position;
    }
    else
    {
        posi_here_negative[2] = 0;             
    }
    posi_prev[0] = imgs_prev[0].position;
    posi_prev[1] = imgs_prev[1].position;
    if(Nphys_prev == 5)
    {
        posi_prev[2] = imgs_prev[2].position;
    }
    else
    {
        posi_prev[2] =0;
    }

    j_negative_candidate = negative_connect_j012
    ( next_idx_out, next_j_out, Nphys_here, Nphys_prev,
    posi_here_negative, posi_prev, here_idx, prev_idx);      // only change next_idx_out[0,1,2], next_j_out[0,1,2]

    if(j_negative_candidate>=0 && j_negative_candidate<order-1)
    {
        c_t posi_prev_positive[1];   // Nphys==3
        c_t posi_here_positive[2];
        posi_prev_positive[0] = imgs_prev[4].position;
        posi_here_positive[0] = imgs_here[4].position;
        posi_here_positive[1] = imgs_here[3].position;
        fail = cross_j012
        ( next_j_out, j_negative_candidate,
        posi_here_positive, posi_prev_positive, here_idx);
    }
    return fail;
}

// j_here_positive = j_next = {4,3} or {4, -1}
__device__ bool source_base_t::connect_next_j34
( int* next_idx_out, int* next_j_out,
const img_t * imgs_here, const img_t * imgs_next, const int here_idx, const int next_idx,
const int Nphys_here, const int Nphys_next) const
{
    c_t posi_here_positive[2];
    c_t posi_next[2];
    int j_positive_candidate;
    bool fail = false;

    posi_here_positive[0] = imgs_here[4].position;
    if(Nphys_here == 5)
    {
        posi_here_positive[1] = imgs_here[3].position;
    }
    else
    {
        posi_here_positive[1] = 0;      // 保险起见的初始化
    }
    posi_next[0] = imgs_next[4].position;
    if(Nphys_next == 5)
    {
        posi_next[1] = imgs_next[3].position;
    }
    else
    {
        posi_next[1] = 0;
    }            

    // positive parity, 1c1 or 2c2
    j_positive_candidate = positive_connect_j34
    ( next_idx_out, next_j_out, Nphys_here, Nphys_next, 
    posi_here_positive, posi_next, here_idx, next_idx);      // only change next_idx_out[3,4], next_j_out[3,4]
    // 检查跨焦散
    // positive parity connect to here
    if(j_positive_candidate>=0 && j_positive_candidate<order-1 )
    {
        c_t posi_next_negative[2];    // Nphys==3
        c_t posi_here_negative[3];
        posi_next_negative[0] = imgs_next[0].position;
        posi_next_negative[1] = imgs_next[1].position;
        posi_here_negative[0] = imgs_here[0].position;
        posi_here_negative[1] = imgs_here[1].position;
        posi_here_negative[2] = imgs_here[2].position;
        fail = cross_j34
        ( next_j_out, j_positive_candidate,
        posi_here_negative, posi_next_negative, here_idx );
    }
    return fail;
}

__device__ f_t source_base_t::delta_s_1
( const c_t & z_here, const c_t & z_next) const
{
    return -((z_next.im + z_here.im) * (z_next.re - z_here.re) - (z_next.re + z_here.re) * (z_next.im - z_here.im)) / 4;
}


__device__ f_t source_base_t::wedge_product
( const c_t &z1, const c_t &z2 ) const
{
    return z1.re*z2.im - z1.im*z2.re;
}


__device__ void source_base_t::deltaS_error_image
( f_t & deltaS, f_t & Error1234, const int j,
    const src_pt_t < f_t > * pt_here, const src_pt_t < f_t > * pt_other,
    const f_t & theta1, const f_t & theta2, const f_t & theta3) const
{
    f_t deltaS_t, deltaS_p1, deltaS_p2, deltaS_p;
    f_t E_1, E_2, E_3, E_4;
    f_t deter;

    c_t zs0,zs1,dzs0,dzs1;
    f_t wedge0,wedge1;
    int next_j;

    next_j = pt_here->next_j[j];
    zs0 = pt_here->images[j].position;
    zs1 = pt_other->images[next_j].position;
    wedge0 = pt_here->wedge[j];
    wedge1 = pt_other->wedge[next_j];
    dzs0 = pt_here->dz[j];
    dzs1 = pt_other->dz[next_j];

    deltaS_t = delta_s_1(zs0,zs1);
    deltaS_p1 = (wedge0 + wedge1) * theta3 /24;
    deltaS_p2 = wedge_product((zs1-zs0) , (dzs1-dzs0)) * theta1 / 12;
    deltaS_p = (deltaS_p1 + deltaS_p2) / 2;
    E_4 = fabs((deltaS_p1 - deltaS_p2));
    deter = theta2 * (dzs0*dzs1).abs(  );
    if(fabs(deter)>2e-12){E_2 = 1.5 * fabs( deltaS_p *  ((zs0-zs1).norm2(  ) / deter -1));}
    else
    {
        E_2=0;
    }        // usually caused by numericial error
    E_1 = fabs((wedge0 - wedge1) * theta3) /48;
    // mode = 0;
    E_3 = 0.1 * fabs(deltaS_p) * theta2;

    deltaS = deltaS_t + deltaS_p;
    // deltaS = deltaS_t;
    Error1234 = E_1 + E_2 + E_3 + E_4;

    return;
}

__device__ void source_base_t::deltaS_error_parity
( f_t * deltaS_out, f_t * deltaS_Err_out, const bool parity,
const int here_idx, const int other_idx,
const src_pt_t < f_t > * pt_here, const src_pt_t < f_t > * pt_other)    const
{
    f_t Qother_src, Qhere;
    f_t theta1, theta2, theta3;
    int next_i;
    f_t deltaS_tp, E_1234;

    Qhere = (pt_here->Q); 
    Qother_src = (pt_other->Q); 

    // theta1 = ?;
    theta1 = Qother_src - Qhere;
    if( parity == 1 && here_idx==0)
        theta1 = Qother_src - 1;
    if( parity == 0 && other_idx==0)
        theta1 = 1 - Qhere;


    theta1 *= 2 * PI;
    theta2 = theta1 * theta1;
    theta3 = theta1 * theta2;

    int jstart = ( parity==1 ? 0 : 3 );
    int jend   = ( parity==1 ? 3 : 5 );
    int jghost = ( parity==1 ? 2 : 3 );

    for(int j=jstart;j<jend;j++)
    {
        if(j==jghost && pt_here->Nphys==3)
        {
            deltaS_out[jghost] = 0;
            deltaS_Err_out[jghost] = 0;
            continue;
        }
        next_i = pt_here->next_idx[j];
        if(here_idx==next_i)
        {
            deltaS_error_cross(deltaS_tp, E_1234, j, pt_here->next_j[j], (j>=3), pt_here);      // ghost_direction = bool (j>=3)
            deltaS_out[j] = deltaS_tp;
            deltaS_Err_out[j] = E_1234;
            continue;
        }
        deltaS_error_image( deltaS_tp, E_1234, j, pt_here, pt_other, theta1, theta2, theta3);

        deltaS_out[j] = deltaS_tp;
        deltaS_Err_out[j] = E_1234;        
    }

    // // parity == 1
    // theta1 = ( here_idx==0 ? Qother_src - 1. : Qother_src - Qhere ) * 2 * PI;


    // // parity==0
    // theta1 = ( other_idx==0 ? 1. - Qhere : Qother_src - Qhere ) * 2 * PI;


}

__device__ void source_base_t::deltaS_error_cross
( f_t & deltaS, f_t & Error1234, const int j0, const int j1, const bool ghost_direction,
const src_pt_t < f_t > * pt_here ) const
{
    // // gd==1, next is ghost, + connect to -, + is j0, - is j1
    // // gd==0, prev is ghost, - connect to +, - is j0, + is j1
    c_t zs[2];
    c_t dzs[2];
    f_t wedge[2];
    // zs[0] connect to zs[1]
    zs[0] = pt_here->images[j0].position;
    zs[1] = pt_here->images[j1].position;
    wedge[0] = pt_here->wedge[j0];
    wedge[1] = pt_here->wedge[j1];
    dzs[0] = pt_here->dz[j0];
    dzs[1] = pt_here->dz[j1];                  

    f_t deltaS_t = delta_s_1(zs[0],zs[1]);
    f_t theta1 = (zs[1]-zs[0]).abs(  ) / sqrt((dzs[0]*dzs[1]).abs(  ));
    f_t theta2 = theta1*theta1;
    f_t theta3 = theta2*theta1;
    int plus = int(!ghost_direction);
    int minus = int(ghost_direction);
    f_t deltaS_p = (wedge[plus] - wedge[minus]) * theta3 /24;
    f_t E_1 = fabs(wedge[plus] + wedge[minus]) * theta3 / 48;
    c_t temp0 = (zs[plus]-zs[minus])*(dzs[plus]-dzs[minus]);
    f_t temp1 = 2* (zs[plus]-zs[minus]).abs(  ) * sqrt((dzs[plus]*dzs[minus]).abs(  ));
    temp0 = ( (ghost_direction==1) ? temp0 + temp1 : temp0 - temp1);
    f_t E_2 = 1.5* (temp0).abs(  ) * theta1;
    f_t E_3 = 0.1 * fabs(deltaS_p) * theta2;
    deltaS = deltaS_t+deltaS_p;
    Error1234 = E_1+E_2+E_3;
    return;
}


__device__ int source_base_t::bisearch_left
(const f_t* values_in_order, const f_t to_insert, const int len) const
{
	// values是第一个点（不是第零个点！）到最后一个点的取值
	// 如果插在第一个点左侧，即0到values[0]之间
	if(to_insert < values_in_order[0])
	{
		return 0;
	}
	else
	{
		int left = 0;;
		int right = len-1;			// values[len-1] = 1; to_insert < 1, 一定不会出界
		int compare_now = len/2;	// 中间开始找，反正也用不了几步

		while( (right-left)>1)
		{
			if(to_insert > values_in_order[compare_now])
			{
				left = compare_now;
			}
			else
			{
				right = compare_now;
			}
			compare_now = (left + right)/2;
		}
		// 上面的left和right是value的指标，但是value[0]对应pt[1]的分度值，需要偏移一点儿
		return left + 1;
	}
}

// return bool fail;
__device__ bool source_base_t::slope_test
( src_pt_t < f_t > * pt_here, src_pt_t < f_t > * pt_other,
const int jj, const int idx_here, bool& Break ) const
{
    Break = false;
    if(pt_here->next_idx[jj] != idx_here)
    {
        return false;
    }

    c_t posi_phys[2];
    c_t posi_ghost[2];

    posi_phys[0] = pt_here->images[jj].position;
    posi_phys[1] = pt_here->images[pt_here->next_j[jj]].position;
    posi_ghost[0] = pt_other->images[2].position;
    posi_ghost[1] = pt_other->images[3].position;
    posi_phys[0] -= posi_phys[1];
    posi_ghost[0] -= posi_ghost[1];

    f_t norm_p = (posi_phys[0]).norm2();
    f_t norm_g = (posi_ghost[0]).norm2();
    if((norm_g>1e-1) ||(norm_p>1e-1))
    {
        if(jj>=3)
            {pt_here->special_Err = (norm_p + norm_g) * 100;}
        else
            {pt_other->special_Err = (norm_p + norm_g) * 100;}
        // printf("trigger!");
        return false;
    }

    f_t temp = posi_phys[0].re * posi_ghost[0].re + posi_phys[0].im * posi_ghost[0].im;
    f_t cos2 = (temp*temp) / (norm_p*norm_g);

    if(cos2>0.1)
    {
        // printf("cos2 too big! loop1, cos2: %.16f, srcidx: %d, idx: %d, j0: %d, j1: %d\n",
        // cos2,blockIdx.x,idx_here,jj,pt_here->next_j[jj]);
        // if(cos2 < 0.9)
        // {
        //     // Break = true;
        // }

        return true;
    }

    return false;
}




////////////////////////////////////////////////////////////
//

__device__ void source_base_t::solve_point_approx(  ) const
{
    const int i_src = threadIdx.x + n_th * blockIdx.x;
    if( i_src >= n_src )
        return;

    auto & src  = pool_center [ i_src ];
    const c_t zeta = src.loc_centre;
    auto  & mag_err = pool_mag [i_src];
    auto  & ext_info = pool_extended [i_src];

    const f_t Rho = src.rho;

    img_t temp_images[order-1];
    int phys_tst = -1;

    if(solve_imgs(temp_images, zeta, false, ext_info.roots_point)){return;}      // solve_imgs return "fail"

    is_physical_all( temp_images , zeta , phys_tst);

    f_t muPS = 0, muQC = 0;
    for(int j=0;j<order-1;j++)
    {
        if(temp_images[j].physical)
        {
            muPS += 1/fabs(jacobian( temp_images[j].position ));
            muQC += fabs(mu_qc(temp_images[j].position) * (Rho+1e-3) * (Rho+1e-3));
        }
    }

    bool is_point_src
     = (muQC * C_Q < RelTol * muPS)
    && ghost_test_all( temp_images, phys_tst, zeta, Rho)
    && safe_distance_test(zeta, Rho);

    if( is_point_src )
    {
        mag_err.mag = muPS;
        mag_err.err = muQC;
        ext_info.SolveSucceed = true;       // default: false
    }

    return;
}


__device__ void source_base_t::init_break_succeed(  ) const
{
    const int i_src = threadIdx.x + n_th * blockIdx.x;
    if( i_src >= n_src )
        return;
    auto  & ext_info = pool_extended [ i_src ];

    ext_info.SolveSucceed = false;
    ext_info.Break = false;
    pool_center[ i_src ].src_area = PI * pool_center[ i_src ].rho * pool_center[ i_src ].rho;
    // ext_info.Ncross = 0;
    ext_info.Area = 0;
    ext_info.Err_A = 0;
    
    auto  * roots_point = pool_extended [ i_src ].roots_point;
    for(int j=0;j<order-1;j++)
    {
        roots_point[j] = 0;
    }
    return;
}

__device__ void source_base_t::margin_set_local  ( local_info_t<f_t>& local_info ) const
{
    const int idx = threadIdx.x;
    auto & margin_pts = local_info.new_pts[ idx ];        // 3n_th 个位置中，中间的 n_th 个
    const auto & src  = (local_info.src_shape);
    f_t Q = f_t(idx) / f_t(n_th);
    c_t marginloc;


    margin_Q(marginloc, src, Q);

    margin_pts.position = marginloc;
    margin_pts.next_src_idx = ((idx+1) & (n_th-1));
    margin_pts.Q = Q;
    margin_pts.special_Err = 0;
    margin_pts.Nphys = -1;
    for(int j=0;j<order-1;j++)
    {
        margin_pts.deltaS[j]=0;
        margin_pts.deltaS_Err[j]=0;
    }
    local_info.new_pts[ ((idx+1) & (n_th-1)) ].prev_src_idx = idx;
    return;
}


__device__ void source_base_t::margin_solve_local( local_info_t<f_t>& local_info ) const
{
    if(local_info.shared_info->Break)
        return;
    // int batchidx = local_info.batchidx;
    // const int idx = batchidx * n_th + threadIdx.x;
    auto & src  = local_info.new_pts [ threadIdx.x ];      // 3n_th 个位置中，中间的 n_th 个
    auto & center = (local_info.src_shape);
    auto & ex_info = (local_info.src_ext);

    img_t temp_images[order-1];
    int phys_tst = -1;

    c_t dz[5];
    f_t wedge[5];
    int temp_j_from_order_j[5];

    // input
    const c_t zeta = src.position;
    c_t loc_centre = center.loc_centre;;
    f_t Rho = center.rho;



    for(int j=0;j<order-1;j++)
    {
        temp_images[j].position = ex_info.roots_point[j];
    }

    //// calculation
    bool skip = margin_solve_cal (temp_images, dz, wedge, phys_tst, temp_j_from_order_j, zeta, loc_centre, Rho);

    //// output
    src.skip = skip;
    if( skip )
    {   
        set_skip_info(src);    
    }
    else
    {
        src.Nphys = phys_tst;
        for(int out_j=0;out_j<order-1;out_j++)
        {
            int temp_j = temp_j_from_order_j[out_j];
            src.dz[out_j] = dz[temp_j];
            src.wedge[out_j] = wedge[temp_j];
            src.images[out_j] = temp_images[temp_j];
        }
    }
    return;
}


__device__ void source_base_t::neighbor3_info_g   ( local_info_t<f_t>& local_info ) const
{
    if(local_info.shared_info->Break)
        return;
    int neighbor3[3];
    const int srcidx = blockIdx.x;
    int batchidx = local_info.batchidx;
    bool changed;

    
    neighbor3[1] = threadIdx.x + batchidx*n_th;
    src_pt_t <f_t> * pt_here = &local_info.new_pts[(neighbor3[1] - batchidx*n_th)];
    if(pt_here->skip)
    {
        local_info.neighbor3[1] = neighbor3[1];
        local_info.neighbor3[0] = -1;
        local_info.neighbor3[2] = -1;
        return;        
    }



    changed = false;                                                                    // changed 机制是当由于 skip 导致的连接顺序修改时，将新的 prev_src_idx 写入（是否有同步问题？）(应该让所有skip的点跳过计算！这样“凡被读入的prev_src_idx，都不会被修改”，因为只有 skip 的点会看它的 prev_src_idx)
    neighbor3[0] = prev_src_idx_local( neighbor3[1] , local_info, &changed );
    if(changed)
        pt_here->prev_src_idx = neighbor3[0];
    changed = false;
    neighbor3[2] = next_src_idx_local( neighbor3[1] , local_info, &changed );
    if(changed)
        pt_here->next_src_idx = neighbor3[2];

    

    local_info.pt_prev = ( (neighbor3[0] < batchidx*n_th) ? pool_margin[ srcidx * n_point_max + neighbor3[0] ] : local_info.new_pts[(neighbor3[0] - batchidx*n_th)]);
    local_info.pt_next = ( (neighbor3[2] < batchidx*n_th) ? pool_margin[ srcidx * n_point_max + neighbor3[2] ] : local_info.new_pts[(neighbor3[2] - batchidx*n_th)]);

    for(int ii=0;ii<3;ii++)
    {
        local_info.neighbor3[ii] = neighbor3[ii];
    }
    return;
}


__device__ void source_base_t::connect_next_local( local_info_t<f_t>& local_info ) const
{
    if(local_info.shared_info->Break)
        return;
    int batchidx = local_info.batchidx;
    int prev_idx = local_info.neighbor3[0]; 
    int idx = local_info.neighbor3[1];
    int next_idx = local_info.neighbor3[2];

    if( local_info.new_pts[ threadIdx.x ].skip ){return;}

    int next_idx_out[5];
    int next_j_out[5];
    int Nphys_prev, Nphys_here, Nphys_next;
    src_pt_t < f_t > * pt_here;
    src_pt_t < f_t > * pt_prev;
    src_pt_t < f_t > * pt_next;

    pt_here = &local_info.new_pts[ threadIdx.x ];
    pt_prev = &local_info.pt_prev;
    pt_next = &local_info.pt_next;


    Nphys_prev = pt_prev->Nphys;
    Nphys_here = pt_here->Nphys;
    Nphys_next = pt_next->Nphys;


    bool fail = false;
    // loop3 = 0    ///////////////////
    if (prev_idx<batchidx*n_th)
    {
        fail = fail || connect_next_j34
        (next_idx_out, next_j_out, 
        pt_prev->images, pt_here->images, prev_idx, idx,
        Nphys_prev, Nphys_here);        
        
        // 输出
        for(int output_j=3;output_j<5;output_j++)
        {
            pt_prev->next_idx[output_j] = next_idx_out[output_j];
            pt_prev->next_j[output_j] = next_j_out[output_j];            
        }
    }
    // loop1 ////////////////////

    fail = fail || connect_prev_j012
    (next_idx_out, next_j_out, 
        pt_here->images, pt_prev->images, idx, prev_idx,
        Nphys_here, Nphys_prev);

    fail = fail || connect_next_j34
    (next_idx_out, next_j_out, 
        pt_here->images, pt_next->images, idx, next_idx,
        Nphys_here, Nphys_next); 

    // 输出
    for(int output_j=0;output_j<5;output_j++)
    {
        pt_here->next_idx[output_j] = next_idx_out[output_j];
        pt_here->next_j[output_j] = next_j_out[output_j];            
    }
    // loop2            ///////////////////

    if(next_idx<batchidx*n_th)
    {
        fail = fail || connect_prev_j012
        (next_idx_out, next_j_out, 
        pt_next->images, pt_here->images, next_idx, idx,
        Nphys_next, Nphys_here);

        // 输出
        for(int output_j=0;output_j<3;output_j++)
        {
            pt_next->next_idx[output_j] = next_idx_out[output_j];
            pt_next->next_j[output_j] = next_j_out[output_j];            
        }
    }

    if(fail)
    {
        local_info.shared_info->Break = true;
    }
    // local_info.src_ext.Break = fail;
}

// __device__ void source_base_t::connect_next_local_test( local_info_t<f_t>& local_info ) const
// {
//     // if(local_info.src_ext.Break)
//     //     return;
//     int batchidx = local_info.batchidx;
//     int prev_idx = local_info.neighbor3[0]; 
//     int idx = local_info.neighbor3[1];
//     int next_idx = local_info.neighbor3[2];

//     if( local_info.new_pts[ threadIdx.x ].skip ){return;}

//     int next_idx_out[5];
//     int next_j_out[5];
//     int Nphys_prev, Nphys_here, Nphys_next;
//     src_pt_t < f_t > * pt_here;
//     src_pt_t < f_t > * pt_prev;
//     src_pt_t < f_t > * pt_next;

//     pt_here = &local_info.new_pts[ threadIdx.x ];
//     pt_prev = &local_info.pt_prev;
//     pt_next = &local_info.pt_next;


//     Nphys_prev = pt_prev->Nphys;
//     Nphys_here = pt_here->Nphys;
//     Nphys_next = pt_next->Nphys;


//     bool fail = false;


//     int* Ncross_ptr = &local_info.shared_info->Ncross;

//     int temp0;
//     fail = false;
//     // 记录跨焦散位置
//     printf("previous temp0: %d, blockIdx: %d, threadIdx: %d\n",*Ncross_ptr,blockIdx.x, threadIdx.x);

//     temp0 = atomicAdd(Ncross_ptr,1);

//     printf("temp0: %d, blockIdx: %d, threadIdx: %d\n",temp0,blockIdx.x, threadIdx.x);

//     // if(fail)
//     //     printf("srcidx: %d, idx: %d\n",blockIdx.x, threadIdx.x);

// }


__device__ void source_base_t::slope_test_local  (local_info_t<f_t>& local_info ) const
{
    if( local_info.new_pts[ threadIdx.x ].skip ){return;}
    // 基本架构接近 connect_next_local
    int batchidx = local_info.batchidx;
    int prev_idx = local_info.neighbor3[0]; 
    int idx = local_info.neighbor3[1];
    int next_idx = local_info.neighbor3[2];
    bool Break;
    src_pt_t < f_t > * pt_here;
    src_pt_t < f_t > * pt_prev;
    src_pt_t < f_t > * pt_next;
    pt_here = &local_info.new_pts[ threadIdx.x ];
    pt_prev = &local_info.pt_prev;
    pt_next = &local_info.pt_next;

    // loop3 = 0, j34   ///////////////////
    if (prev_idx<batchidx*n_th)
    {
        for(int jj=3;jj<5;jj++)
        {
            if(slope_test(pt_prev, pt_here, jj, prev_idx, Break))
            {
                local_info.shared_info->Break = true;
                if(!Break)
                {
                    int temp0 = atomicAdd(&local_info.shared_info->Ncross, 1);
                    local_info.cross_info[ temp0 ].idx_cross = prev_idx;
                    local_info.cross_info[ temp0 ].j_cross = jj;
                    local_info.cross_info[ temp0 ].j1 = pt_prev->next_j[jj];
                    local_info.cross_info[ temp0 ].additional = false;                    
                }
                else
                {
                    local_info.shared_info->Break = true;
                }
            }
        }
    }

    // loop3 = 1, j012
    for(int jj=0;jj<3;jj++)
    {
        if(slope_test(pt_here, pt_prev, jj, idx, Break))
        {
            local_info.shared_info->Break = true;
            if(!Break)
            {
                int temp0 = atomicAdd(&local_info.shared_info->Ncross, 1);
                local_info.cross_info[ temp0 ].idx_cross = idx;
                local_info.cross_info[ temp0 ].j_cross = jj;
                local_info.cross_info[ temp0 ].j1 = pt_here->next_j[jj];
                local_info.cross_info[ temp0 ].additional = false;
            }
            else
            {
                local_info.shared_info->Break = true;
            }
        }

    } 

    // loop3 = 1, j34
    for(int jj=3;jj<5;jj++)
    {
        if(slope_test(pt_here, pt_next, jj, idx, Break))
        {
            local_info.shared_info->Break = true;
            if(!Break)
            {            
                int temp0 = atomicAdd(&local_info.shared_info->Ncross, 1);
                local_info.cross_info[ temp0 ].idx_cross = idx;
                local_info.cross_info[ temp0 ].j_cross = jj;
                local_info.cross_info[ temp0 ].j1 = pt_here->next_j[jj];
                local_info.cross_info[ temp0 ].additional = false;
            }
            else
            {
                local_info.shared_info->Break = true;
            }
        }
    }   

    // loop3 = 2, j012   ///////////////////
    if(next_idx<batchidx*n_th)
    {
        for(int jj=0;jj<3;jj++)
        {
            if(slope_test(pt_next, pt_here, jj, next_idx, Break))
            {
                local_info.shared_info->Break = true;
                if(!Break)
                {
                    int temp0 = atomicAdd(&local_info.shared_info->Ncross, 1);
                    local_info.cross_info[ temp0 ].idx_cross = next_idx;
                    local_info.cross_info[ temp0 ].j_cross = jj;
                    local_info.cross_info[ temp0 ].j1 = pt_next->next_j[jj];
                    local_info.cross_info[ temp0 ].additional = false;
                }
                else
                {
                    local_info.shared_info->Break = true;
                }
            }
        }
    }    
}



__device__ void source_base_t::slope_detector_g    ( local_info_t<f_t>& local_info ) const
{
    if(local_info.shared_info->Break)
        return;
    if(threadIdx.x != 0)            // 0 号单线程执行
        return;

    const int srcidx = blockIdx.x;

    int Ncross = local_info.shared_info->Ncross;
    if(Ncross<=0)
        return;


    int idx_cross_p;        // physical
    bool ghost_direction;     // 1 means next, 0 means previous
    int idx_cross_g;       // ghost
    complex_t<f_t> posi_p[2]; 
    complex_t<f_t> posi_g[2];
    int n_ghost;
    f_t cos2;
    f_t temp;
    int j_test_p[2];
    int temp_j;
    int iter_i;
    bool unclear;       // if cos2 > 0.1 but < 0.9
    bool test_pass;

    int idx_change[n_th+1];       // ghost images candidates' idx
    int j0[n_th+1];
    int j1[n_th+1];
    f_t norms[n_th+1];
    f_t norm_g;
    f_t& norm_p = norm_g;

    int* j_test_g = j_test_p;

    bool additional_piece;
    bool tempbool;

    // src_pt_t<f_t> * pt_phys;
    // src_pt_t<f_t> * pt_ghost;


    for(int i=0;i<Ncross;i++)
    {      
        additional_piece = false;
        test_pass=0;
        unclear = 0;        // set to zero for each cross caustic

        idx_cross_p = local_info.cross_info[ i ].idx_cross;      // 第一次的一定在 local
        j_test_p[0] = local_info.cross_info[ i ].j_cross;
        j_test_p[1] = local_info.cross_info[ i ].j1;

        ghost_direction = ( (j_test_p[0]<3) ? 0 : 1 );
        idx_cross_g = ( (j_test_p[0]<3) ? prev_src_idx_local_g(idx_cross_p, local_info, blockIdx.x) : next_src_idx_local_g(idx_cross_p, local_info, blockIdx.x) );

    
        if(pool_margin[srcidx * n_point_max + idx_cross_p].Nphys==5)       // 这个点判断之初肯定是N5，但是有一种情况（虚像碎片连成环）会在之前的i循环里修改掉
        // if(false)
        {
            
            posi_p[0] = pool_margin[srcidx * n_point_max + idx_cross_p].images[j_test_p[0]].position;
            posi_p[1] = pool_margin[srcidx * n_point_max + idx_cross_p].images[j_test_p[1]].position;

            if(pool_margin[srcidx * n_point_max + idx_cross_g].Nphys==5)
            {
                local_info.cross_info[i].additional = true;
                continue;                
            }
            posi_g[0] = pool_margin[srcidx * n_point_max + idx_cross_g].images[2].position;
            posi_g[1] = pool_margin[srcidx * n_point_max + idx_cross_g].images[3].position;

            iter_i = 0;
            posi_p[0] = posi_p[1] - posi_p[0];      // vector between two points
            posi_g[0] = posi_g[1] - posi_g[0];
            norm_g = (posi_g[0]).norm2(  );
            norms[iter_i] = (posi_p[0]).norm2(  );

            temp = posi_p[0].re * posi_g[0].re + posi_p[0].im * posi_g[0].im;
            cos2 = (temp*temp) / (norms[iter_i]*norm_g);

            if(cos2<=0.1){test_pass = true;}
            else{
                if(cos2<0.9)
                {
                    unclear = true;
                }
            }

            idx_change[iter_i] = idx_cross_p;    // 需要在最后修改的虚点。如果通过检测，iter_i = 0，只会修改指标<0的点，即一个都不改
            j0[iter_i] = j_test_p[0];
            j1[iter_i] = j_test_p[1];
            // cos2s[iter_i] = cos2;                        

            while((!test_pass) && (iter_i<n_th+1))        // 夹角不够垂直，没有jump，说明当前的physical也是ghost（假设只有N5判断会出错，N3判断严格正确）
            {
                iter_i++;
                if(iter_i == n_th+1){break;}
                // 是上面注释掉的if-else的整合版
                j_test_p[1] = pool_margin[srcidx * n_point_max + idx_cross_p].next_j[j_test_p[1]];
                if(ghost_direction)
                {
                    idx_cross_p = prev_src_idx_local_g(idx_cross_p, local_info, blockIdx.x);
                }
                else
                {       
                    idx_cross_p = next_src_idx_local_g(idx_cross_p, local_info, blockIdx.x);                  
                }
                if(pool_margin[srcidx * n_point_max + idx_cross_p].Nphys!=5)
                {
                    // 该片段已走到尽头，也没通过
                    test_pass = true;
                    additional_piece = true;
                    break;
                }
                for(int j=0;j<5;j++)        // 连接：因为nextj是单向的，有一个方向需要遍历查找
                {
                    if(pool_margin[srcidx * n_point_max + idx_cross_p].images[j].parity != ghost_direction)
                    {
                        temp_j = pool_margin[srcidx * n_point_max + idx_cross_p].next_j[j];
                        if((temp_j == j_test_p[0]))
                        {
                            j_test_p[0] = j;
                            break;
                        }                                    
                    }

                }
                posi_p[0] = pool_margin[srcidx * n_point_max + idx_cross_p].images[j_test_p[0]].position;
                posi_p[1] = pool_margin[srcidx * n_point_max + idx_cross_p].images[j_test_p[1]].position;                        

                posi_p[0] = posi_p[1] - posi_p[0];
                norms[iter_i] = (posi_p[0]).norm2(  );

                if(norms[iter_i]>1e-1)        // too far away, pass
                {
                    test_pass = true;
                    // pool_margin[srcidx * n_point_max + idx_cross_p].special = 100.*norms[iter_i];   // 强行要求在附近加点

                    // 目的是在 ghost 和 physical 之间加点
                    if(ghost_direction)     // next is ghost
                    {
                        pool_margin[srcidx * n_point_max + idx_cross_p].special_Err = (norms[iter_i])*100;      // 权重指的是 here 和 next 的区间，对应 physical 和 ghost
                    }
                    else                    // previous is ghost
                    {
                        pool_margin[srcidx * n_point_max + idx_change[iter_i-1]].special_Err = (norms[iter_i])*100;      // 权重指的是 here 和 next 的区间，对应 ghost 和 physical
                    }

                    break;
                }

                temp = posi_p[0].re * posi_g[0].re + posi_p[0].im * posi_g[0].im;
                cos2 = (temp*temp) / (norms[iter_i]*norm_g);

                idx_change[iter_i] = idx_cross_p;    // 需要在最后修改的虚点
                j0[iter_i] = j_test_p[0];
                j1[iter_i] = j_test_p[1];
                // cos2s[iter_i] = cos2;  


                if(unclear)
                {
                    // unclear 区判断标准：距离回升
                    if(norms[iter_i] > norms[iter_i-1])
                    {
                        test_pass=true;
                        // printf("unclear trigger\n");
                    }
                }
                else
                {
                    if(cos2<=0.1){test_pass=true;}
                    else{
                        if(cos2<0.9)
                        {
                            unclear=true;
                            // unclear_start = iter_i;       // +1 因为iter_i++放在后面了
                        }
                    }
                    // modify_len++;
                }
            }


            if((iter_i != n_th+1) && (norms[iter_i]<=1e-1))    // 成功了才改，不成功反向检查
            {
                if(!additional_piece)
                {
                    // connect modified points
                    if(!unclear)
                    {
                        pool_margin[srcidx * n_point_max + idx_cross_p].next_idx[j_test_p[0]] = idx_cross_p;
                        pool_margin[srcidx * n_point_max + idx_cross_p].next_j[j_test_p[0]] = j_test_p[1];

                        // if(ghost_direction)
                        // {
                        //     srcs[srcidx].idx_cross[i] = idx_cross_p;
                        // }
                        // else{
                        //     srcs[srcidx].idx_cross[i] = -idx_cross_p-1;
                        // }
                        local_info.cross_info[i].idx_cross = idx_cross_p;
                        local_info.cross_info[i].additional = false;
                        local_info.cross_info[i].j_cross = j_test_p[0];
                    }
                    else
                    {
                        if(iter_i>1)
                        {                                 
                            pool_margin[srcidx * n_point_max + idx_change[iter_i-1]].next_idx[j0[iter_i-1]] = idx_change[iter_i-1];
                            pool_margin[srcidx * n_point_max + idx_change[iter_i-1]].next_j[j0[iter_i-1]] = j_test_p[1];
                            // if(ghost_direction)
                            // {
                            //     srcs[srcidx].idx_cross[i] = idx_change[iter_i-1];
                            // }
                            // else{
                            //     srcs[srcidx].idx_cross[i] = -idx_change[iter_i-1]-1;
                            // }
                            local_info.cross_info[i].idx_cross = idx_change[iter_i-1];
                            local_info.cross_info[i].j_cross = j0[iter_i-1]; 
                            local_info.cross_info[i].additional = false;                                   
                        }
                        // else{}       // else, 说明第一个点就unclear，查询下一个点发现没问题，所以谁也不修改
                        // iter_i = 0: 有unclear至少为1
                        // iter_i = 1: 一开始进了unclear，检测一次发现下一个更实，所以0号是实的
                        
                    }                            
                }
                else        // additional piece = true
                {
                    // 需要帮人家再补起来
                    // 就是，进入broken之前，连接出broken之后
                    // 

                    // 找到连接开头的那个点
                    int prev_pt_idx = -1;
                    int prev_pt_j = -1;
                    bool found_prev = false;

                    // 试试是不是自己连自己
                    for(int prevj=0;prevj<5;prevj++)
                    {
                        prev_pt_idx = idx_change[0];
                        if((pool_margin[srcidx * n_point_max + prev_pt_idx].next_idx[prevj] == prev_pt_idx) && (pool_margin[srcidx * n_point_max + prev_pt_idx].next_j[prevj] == j0[0]))
                        {
                            prev_pt_j = prevj;
                            found_prev = true;
                        }
                    }
                    // 没找到，说明是前后连过来的，取决于 ghost_direction (这玩意好像就是手征啊，我存这么个莫名其妙的东西干啥？)
                    if(!found_prev)
                    {
                        for(int prevj=0;prevj<5;prevj++)
                        {
                            // prev_pt_idx = PrevSrcIdx(srcs[srcidx].margin_pts, idx_change[0]);
                            prev_pt_idx = prev_src_idx_local_g(idx_change[0], local_info, blockIdx.x);
                            if((pool_margin[srcidx * n_point_max + prev_pt_idx].next_idx[prevj] == prev_pt_idx) && (pool_margin[srcidx * n_point_max + prev_pt_idx].next_j[prevj] == j0[0]))
                            {
                                prev_pt_j = prevj;
                                found_prev = true;
                            }
                        }
                    }
                    if(!found_prev)
                    {
                        for(int prevj=0;prevj<5;prevj++)
                        {
                            // prev_pt_idx = NextSrcIdx(srcs[srcidx].margin_pts, idx_change[0]);
                            prev_pt_idx = next_src_idx_local_g(idx_change[0], local_info, blockIdx.x);   
                            if((pool_margin[srcidx * n_point_max + prev_pt_idx].next_idx[prevj] == prev_pt_idx) && (pool_margin[srcidx * n_point_max + prev_pt_idx].next_j[prevj] == j0[0]))
                            {
                                prev_pt_j = prevj;
                                found_prev = true;
                            }
                        }
                    }

                    pool_margin[srcidx * n_point_max + prev_pt_idx].next_idx[prev_pt_j] = pool_margin[srcidx * n_point_max + idx_change[0]].next_idx[j1[0]];
                    pool_margin[srcidx * n_point_max + prev_pt_idx].next_j[prev_pt_j] = pool_margin[srcidx * n_point_max + idx_change[0]].next_j[j1[0]];


                    local_info.cross_info[i].additional = true;
                    pool_margin[srcidx * n_point_max + idx_change[0]].Nphys = 3;       // 实际上下面也有，只是强调一下，增加鲁棒性
                }



                // printf("%d,%d,%d,iters: %d, bool: %d, cos2: %.16f, j0: %d, j1: %d\n",idx_cross_p,idx_cross_g,batchidx,iter_i,ghost_direction,cos2,j_test[0],j_test[1]);

                // 修改错误的信息
                // 如果除非是unclear段的成功位置，那里旧点也是实像
                // 其他情况旧点都是虚像
                for(int iii=0;iii<iter_i - unclear;iii++)       // unclear==0: iter_i 是实像，故修改0～iter_i-1; unclear: 再少一个
                {
                    // if(srcidx==16 && batchidx==8)
                    // {
                    //     printf("i: %d, idx = %d, j0 = %d, j1 = %d\n",i,idx_change[iii],j0[iii],j1[iii]);
                    // }
                    // change the connection of wrong points
                    // j_test[0] is always the chain to connect another
                    pool_margin[srcidx * n_point_max + idx_change[iii]].images[j0[iii]].physical = false;
                    pool_margin[srcidx * n_point_max + idx_change[iii]].images[j1[iii]].physical = false;
                    pool_margin[srcidx * n_point_max + idx_change[iii]].Nphys = 3;
                    pool_margin[srcidx * n_point_max + idx_change[iii]].next_idx[j0[iii]] = -2;
                    pool_margin[srcidx * n_point_max + idx_change[iii]].next_idx[j1[iii]] = -2;

                    pool_margin[srcidx * n_point_max + idx_change[iii]].deltaS[j0[iii]] = 0;
                    pool_margin[srcidx * n_point_max + idx_change[iii]].deltaS[j1[iii]] = 0;
                    pool_margin[srcidx * n_point_max + idx_change[iii]].deltaS_Err[j0[iii]] = 0;
                    pool_margin[srcidx * n_point_max + idx_change[iii]].deltaS_Err[j1[iii]] = 0;

                    // pool_margin[srcidx * n_point_max + idx_change[iii]].deltaS_new[j0[iii]] = 0;
                    // pool_margin[srcidx * n_point_max + idx_change[iii]].deltaS_new[j1[iii]] = 0;
                    // pool_margin[srcidx * n_point_max + idx_change[iii]].Err_new[j0[iii]] = 0;
                    // pool_margin[srcidx * n_point_max + idx_change[iii]].Err_new[j1[iii]] = 0;
                    // if constexpr (DetailRecord)
                    // #if DetailedRecord
                    // // {
                    //     pool_margin[srcidx * n_point_max + idx_change[iii]].E1[j0[iii]] = 0;
                    //     pool_margin[srcidx * n_point_max + idx_change[iii]].E1[j1[iii]] = 0;
                    //     pool_margin[srcidx * n_point_max + idx_change[iii]].E2[j0[iii]] = 0;
                    //     pool_margin[srcidx * n_point_max + idx_change[iii]].E2[j1[iii]] = 0;
                    //     pool_margin[srcidx * n_point_max + idx_change[iii]].E3[j0[iii]] = 0;
                    //     pool_margin[srcidx * n_point_max + idx_change[iii]].E3[j1[iii]] = 0;                                
                    // // }
                    // #endif


                    // if(srcidx==13)
                    // {
                    //     printf("idx to N3: %d\n",idx_change[iii]);
                    // }                                                                                  
                }
            }
            else      // 认为只有虚像N3错判为实N5失败了，实际上有实像N5错判为虚像N3
            {

                // if(srcidx==280)
                // {
                //     printf("here! idx_cross_p: %d\n",idx_change[0]);
                // }
                // 标记方法为：idx_change 系列数组，0号还是srcs[srcidx].idx_cross[i]，从1号起是反向的N3点，注意j0和j1分别连成链
                test_pass = false;
                iter_i = 1;
                idx_cross_p = idx_change[0];        // 前面拿过了
                posi_p[0] = pool_margin[srcidx * n_point_max + idx_cross_p].images[j0[0]].position;
                posi_p[1] = pool_margin[srcidx * n_point_max + idx_cross_p].images[j1[0]].position;
                                        

                idx_cross_g = (ghost_direction) ? next_src_idx_local_g(idx_cross_p, local_info, blockIdx.x) \
                                                : prev_src_idx_local_g(idx_cross_p, local_info, blockIdx.x) ;
                idx_change[iter_i] = idx_cross_g;

                // n_ghost=0;
                // for(int j=0;j<5;j++)
                // {
                //     if((!pool_margin[srcidx * n_point_max + idx_cross_g].images[j].physical) && (n_ghost<2))
                //     {
                //         posi_g[n_ghost] = pool_margin[srcidx * n_point_max + idx_cross_g].images[j].position;
                //         j_test_g[n_ghost] = j;
                //         n_ghost++;                                
                //     }
                // }           // 这个不用检查了，实际上前半截已经算过了，只是没存，该报的错都报了
                posi_g[0] = pool_margin[srcidx * n_point_max + idx_cross_g].images[2].position;
                posi_g[1] = pool_margin[srcidx * n_point_max + idx_cross_g].images[3].position;


                posi_p[0] = posi_p[1] - posi_p[0];      // vector between two points
                posi_g[0] = posi_g[1] - posi_g[0];
                norm_p = (posi_p[0]).norm2(  );
                norms[iter_i] = (posi_g[0]).norm2(  );
                temp = posi_p[0].re * posi_g[0].re + posi_p[0].im * posi_g[0].im;
                // if(temp>0)
                // {
                //     j0[iter_i] = j_test_g[0]; j1[iter_i] = j_test_g[1];
                // }
                // else
                // {
                //     j0[iter_i] = j_test_g[1]; j1[iter_i] = j_test_g[0];
                // }
                tempbool = (temp>0);
                j0[iter_i] = j_test_g[int(!tempbool)];
                j1[iter_i] = j_test_g[int( tempbool)];

                cos2 = (temp*temp) / (norms[iter_i-1]*norms[iter_i]);
                if(cos2<=0.1){test_pass = true;}
                else{
                    if(cos2<0.9)
                    {
                        printf("unclear! srcidx: %d\n",blockIdx.x);
                        unclear = true;
                        // unclear_start = 0;
                    }
                }     
                while((!test_pass) && (iter_i<n_th+1-1))
                {
                    // 只改posi_g
                    iter_i++;
                    if(iter_i == n_th+1){break;}
                    idx_cross_p = idx_cross_g;
                    idx_cross_g = (ghost_direction) ? next_src_idx_local_g(idx_cross_p, local_info, blockIdx.x) \
                                                    : prev_src_idx_local_g(idx_cross_p, local_info, blockIdx.x) ;
                    idx_change[iter_i] = idx_cross_g;


                    if(pool_margin[srcidx * n_point_max + idx_cross_g].Nphys!=3)
                    {
                        // 该片段已走到尽头，也没通过
                        test_pass = true;
                        additional_piece = true;

                        break;
                    }



                    n_ghost=0;
                    for(int j=0;j<5;j++)
                    {
                        if((!pool_margin[srcidx * n_point_max + idx_cross_g].images[j].physical) && (n_ghost<2))
                        {
                            posi_g[n_ghost] = pool_margin[srcidx * n_point_max + idx_cross_g].images[j].position;
                            j_test_g[n_ghost] = j;
                            n_ghost++;                                
                        }
                    }
                    posi_g[0] = posi_g[1] - posi_g[0];
                    norms[iter_i] = (posi_g[0]).norm2(  );
                    // if(n_ghost!=2)
                    // {
                    //     // if(n_ghost==0)      // 到头了！
                    //     // {
                    //     //     // 实际上应该在前面  if(pool_margin[srcidx * n_point_max + idx_cross_g].Nphys!=3) 就触发了，不应该
                    //     // }
                    //     // else
                    //     // {
                    //     if(!muted)
                    //     {
                    //         printf("fatal error: ghost images number wrong, srcidx: %d, batchidx: %d, idx_cross_g: %d, Nghost: %d\n",srcidx,batchidx,idx_cross_g,n_ghost);
                    //     }
                    //     srcs[srcidx].Break = true;
                    //     srcs[srcidx].SolveSucceed = 1;          // 用于跳过后续计算，会在最后一个SumArea设置为fail
                    //     // #if DetailedRecord
                    //     //     srcs[srcidx].points_used = (batchidx+1) * n_th+1;
                    //     // #endif                                    
                    //     // }
                    // }


                    if(norms[iter_i]>1e-1)        // too far away, pass
                    {
                        test_pass = true;
                        // pool_margin[srcidx * n_point_max + idx_cross_p].special = 100.*norms[iter_i];   // 强行要求在附近加点

                        // 目的是在 ghost 和 physical 之间加点
                        if(ghost_direction)     // next is ghost
                        {
                            pool_margin[srcidx * n_point_max + idx_change[iter_i-1]].special_Err = (norms[iter_i])*100;      // 权重指的是 here 和 next 的区间，对应 physical 和 ghost
                        }
                        else                    // previous is ghost
                        {
                            pool_margin[srcidx * n_point_max + idx_change[ iter_i ]].special_Err = (norms[iter_i])*100;      // 权重指的是 here 和 next 的区间，对应 ghost 和 physical
                        }
                        break;
                    }

                    temp = posi_p[0].re * posi_g[0].re + posi_p[0].im * posi_g[0].im;
                    cos2 = (temp*temp) / (norms[iter_i]*norm_p);
                    // printf("iter_i: %d, cos2: %.12f\n",iter_i,cos2);

                    tempbool = (temp>0);
                    j0[iter_i] = j_test_g[int(!tempbool)];
                    j1[iter_i] = j_test_g[int( tempbool)];    

                    if(unclear)
                    {
                    // unclear 区判断标准：距离回升
                        if(norms[iter_i] > norms[iter_i-1])
                        {
                            test_pass=true;
                            // printf("unclear trigger\n");
                        }
                    }
                    else
                    {
                        if(cos2<=0.1){test_pass=true;}
                        else{
                            if(cos2<0.9)
                            {
                                printf("unclear! srcidx: %d\n",blockIdx.x);
                                unclear=true;
                            }
                        }
                    }

                }       // end while

                // printf("iter_i: %d\n",iter_i);
                if((iter_i != n_th+1) && (norms[iter_i]<=1e-1))    // 成功了才改，还不成功就报错
                {
                    // if(srcidx==654 && batchidx==12)
                    // {
                    //     printf("i: %d, idx_cross_marked: %d\n",i,srcs[srcidx].idx_cross[i]);
                    // }

                    // 修改错误的信息
                    // 如果除非是unclear段的成功位置，那里旧点也是实像
                    // 其他情况旧点都是虚像
                    tempbool = pool_margin[srcidx * n_point_max + idx_change[0]].images[j0[0]].parity;
                    pool_margin[srcidx * n_point_max + idx_change[0]].next_idx[j0[0]] = idx_change[0+1];
                    pool_margin[srcidx * n_point_max + idx_change[0]].next_j[j0[0]]   = j0[0+1];

                    for(int iii=1;iii<iter_i;iii++)       // unclear==0: iter_i 是实像，故修改0～iter_i-1; unclear: 再少一个
                    {
                        // if(srcidx==16 && batchidx==8)
                        // {
                        //     printf("i: %d, idx = %d, j0 = %d, j1 = %d\n",i,idx_change[iii],j0[iii],j1[iii]);
                        // }
                        // change the connection of wrong points
                        // j_test[0] is always the chain to connect another
                        pool_margin[srcidx * n_point_max + idx_change[iii]].images[j0[iii]].physical = true;
                        pool_margin[srcidx * n_point_max + idx_change[iii]].images[j1[iii]].physical = true;
                        pool_margin[srcidx * n_point_max + idx_change[iii]].Nphys = 5;
                        pool_margin[srcidx * n_point_max + idx_change[iii]].images[j0[iii]].parity =  tempbool;
                        pool_margin[srcidx * n_point_max + idx_change[iii]].images[j1[iii]].parity = !tempbool;

                        // ghost = 1: idx_change along next; parity = 0: 0 chain connect to next
                        pool_margin[srcidx * n_point_max + idx_change[iii]].next_idx[j0[iii]] = idx_change[iii+1];
                        pool_margin[srcidx * n_point_max + idx_change[iii]].next_j[j0[iii]]   = j0[iii+1];
                        // ghost = 1: idx_change along next; parity = 1: 1 chain connect to next
                        pool_margin[srcidx * n_point_max + idx_change[iii]].next_idx[j1[iii]] = idx_change[iii-1];
                        pool_margin[srcidx * n_point_max + idx_change[iii]].next_j[j1[iii]]   = j1[iii-1];

                                                                        
                    }

                    if(!additional_piece)
                    {
                        // connect modified points

                        
                        // pool_margin[srcidx * n_point_max + idx_cross_p].next_idx[j_test_p[0]] = idx_cross_p;
                        // pool_margin[srcidx * n_point_max + idx_cross_p].next_j[j_test_p[0]] = j_test_p[1];

                        // if(ghost_direction)
                        // {
                        //     srcs[srcidx].idx_cross[i] = idx_cross_p;
                        // }
                        // else{
                        //     srcs[srcidx].idx_cross[i] = -idx_cross_p-1;
                        // }                                

                        // srcs[srcidx].additional[i] = false;
                        
                        // srcs[srcidx].j_cross[i] = j_test_p[0];

                        // 自己封口

                        pool_margin[srcidx * n_point_max + idx_change[iter_i-1]].next_idx[j0[iter_i-1]] = idx_change[iter_i-1];
                        pool_margin[srcidx * n_point_max + idx_change[iter_i-1]].next_j[j0[iter_i-1]] = j1[iter_i-1];

                        // srcs[srcidx].idx_cross[i] = (ghost_direction)   ? idx_change[iter_i-1]    \
                                                                        // : -idx_change[iter_i-1]-1 ;
                        local_info.cross_info[i].idx_cross = idx_change[iter_i-1];
                        local_info.cross_info[i].additional = false;
                        local_info.cross_info[i].j_cross = j0[iter_i-1];
                
                    }
                    else        // additional piece = true
                    {
                        // 把两端的封口重新打开，想办法连上
                        // 重新进行一次简化的findnext？但findnext不是以srcpt为单位进行的...

                        // if(srcidx==654&&batchidx==12)
                        // {
                        //     printf("idx_change[iter_i-1]: %d, idx_change[iter_i]: %d\n",idx_change[iter_i-1], idx_change[iter_i]);
                        // }

                        // 示例：
                        // idx:             758, 759, 760, 761, 762
                        // real N:            5,   5,   5,   5,   3
                        // judged N:          5,   3,   5,   5,   3
                        // idx_change[iter_i-1] == 759, idx_change[iter_i] == 760, additional piece

                        for(int jjjj=0;jjjj<5;jjjj++)
                        {
                            if(pool_margin[srcidx * n_point_max + idx_change[iter_i]].next_idx[jjjj]==idx_change[iter_i])      // 如果 760 的某个 j 连接自己
                            {
                                if(pool_margin[srcidx * n_point_max + idx_change[iter_i]].images[jjjj].parity == pool_margin[srcidx * n_point_max + idx_change[iter_i-1]].images[j1[iter_i-1]].parity)       // 如果 760 发出自连的这个 j 手征和 759 被连接的还是一样的
                                {
                                    // 就是它了（仅限2体）
                                    // 被 jjjj 连的像现在被 759 的 j0 连
                                    pool_margin[srcidx * n_point_max + idx_change[iter_i-1]].next_idx[j0[iter_i-1]] = idx_change[iter_i];
                                    pool_margin[srcidx * n_point_max + idx_change[iter_i-1]].next_j[j0[iter_i-1]] = pool_margin[srcidx * n_point_max + idx_change[iter_i]].next_j[jjjj];
                                    // 760 发起自连的像现在连 759 的 j1
                                    pool_margin[srcidx * n_point_max + idx_change[iter_i]].next_idx[jjjj] = idx_change[iter_i-1];
                                    pool_margin[srcidx * n_point_max + idx_change[iter_i]].next_j[jjjj] = j1[iter_i];
                                    break;

                                }
                            }
                            
                        }

                        local_info.cross_info[i].additional = true;


                    }

                }

            }

        // printf("%d,%d\n",idx_cross_p,idx_cross_g);
        }
        else{       // 这个片段在前面就被改掉了，现在是N3，整个小片段都是虚的
            // srcs[srcidx].additional[i] = true;
            local_info.cross_info[ i ].additional = true;
        }
    }


}


__device__ void source_base_t::area_err_local    ( local_info_t<f_t>& local_info ) const
{
    if(local_info.shared_info->Break)
        return;
    int batchidx = local_info.batchidx;
    bool parity;
    f_t deltaS_out[5], deltaS_Err_out[5];
    int idx_here = local_info.neighbor3[1];
    int idx_prev   = local_info.neighbor3[0]; 
    int idx_next   = local_info.neighbor3[2];
    src_pt_t < f_t > * pt_here;
    src_pt_t < f_t > * pt_prev;
    src_pt_t < f_t > * pt_next;

    pt_here = &local_info.new_pts[ threadIdx.x ];
    pt_prev = &local_info.pt_prev;
    pt_next = &local_info.pt_next;


    if(local_info.new_pts[ threadIdx.x ].skip){ return; }

    // loop3 == 0   ///////////////////
    if(idx_prev<batchidx*n_th)
    {
        // positive
        parity = 0;
        deltaS_error_parity(deltaS_out, deltaS_Err_out, parity, idx_prev, idx_here, pt_prev, pt_here);

        for(int output_j=3; output_j<5;output_j++)
        {
            pt_prev->deltaS[output_j] = deltaS_out[output_j];
            pt_prev->deltaS_Err[output_j] = deltaS_Err_out[output_j];
        }
    }

    // loop3 == 1       ////////////////////////

    // negative
    parity = 1;                      
    deltaS_error_parity(deltaS_out, deltaS_Err_out, parity, idx_here, idx_prev, pt_here, pt_prev);                    
    // positive
    parity = 0;
    deltaS_error_parity(deltaS_out, deltaS_Err_out, parity, idx_here, idx_next, pt_here, pt_next);

    for(int output_j=0; output_j<5;output_j++)
    {
        pt_here->deltaS[output_j] = deltaS_out[output_j];
        pt_here->deltaS_Err[output_j] = deltaS_Err_out[output_j];
    }

    // loop3 == 2 /////////////////////
    if(idx_next<batchidx*n_th)
    {
        // negative
        parity = 1;                        
        deltaS_error_parity(deltaS_out, deltaS_Err_out, parity, idx_next, idx_here, pt_next, pt_here);

        for(int output_j=0; output_j<3;output_j++)
        {
            pt_next->deltaS[output_j] = deltaS_out[output_j];
            pt_next->deltaS_Err[output_j] = deltaS_Err_out[output_j];
        }
    }
}

__device__ void source_base_t::area_err_local_test ( local_info_t<f_t>& local_info ) const
{
    if(local_info.shared_info->Break)
        return;
    int batchidx = local_info.batchidx;
    bool parity;
    f_t deltaS_out[5], deltaS_Err_out[5];
    int idx_here = local_info.neighbor3[1];
    int idx_prev   = local_info.neighbor3[0]; 
    int idx_next   = local_info.neighbor3[2];
    src_pt_t < f_t > * pt_here;
    src_pt_t < f_t > * pt_prev;
    src_pt_t < f_t > * pt_next;

    pt_here = &local_info.new_pts[ threadIdx.x ];
    pt_prev = &local_info.pt_prev;
    pt_next = &local_info.pt_next;


    if(local_info.new_pts[ threadIdx.x ].skip){ return; }

    // loop3 == 0   ///////////////////
    if(idx_prev<batchidx*n_th)
    {
        // positive
        parity = 0;
        deltaS_error_parity(deltaS_out, deltaS_Err_out, parity, idx_prev, idx_here, pt_prev, pt_here);

        for(int output_j=3; output_j<5;output_j++)
        {
            pt_prev->deltaS[output_j] = deltaS_out[output_j];
            pt_prev->deltaS_Err[output_j] = deltaS_Err_out[output_j];
        }
    }

    // loop3 == 1       ////////////////////////

    // negative
    parity = 1;                      
    deltaS_error_parity(deltaS_out, deltaS_Err_out, parity, idx_here, idx_prev, pt_here, pt_prev);                    
    // positive
    parity = 0;
    deltaS_error_parity(deltaS_out, deltaS_Err_out, parity, idx_here, idx_next, pt_here, pt_next);

    for(int output_j=0; output_j<5;output_j++)
    {
        pt_here->deltaS[output_j] = deltaS_out[output_j];
        pt_here->deltaS_Err[output_j] = deltaS_Err_out[output_j];
    }

    // loop3 == 2 /////////////////////
    if(idx_next<batchidx*n_th)
    {
        // negative
        parity = 1;                        
        deltaS_error_parity(deltaS_out, deltaS_Err_out, parity, idx_next, idx_here, pt_next, pt_here);

        for(int output_j=0; output_j<3;output_j++)
        {
            pt_next->deltaS[output_j] = deltaS_out[output_j];
            pt_next->deltaS_Err[output_j] = deltaS_Err_out[output_j];
        }
    }
}

// __device__ void source_base_t::area_err_c_local( local_info_t<f_t>& local_info ) const
// {
//     if(local_info.shared_info->Break)
//         return;
//     const int Ncross = local_info.shared_info->Ncross;
//     if(threadIdx.x >= Ncross){return;}

//     cross_info_t * cross_info = local_info.cross_info;
//     if(cross_info[ threadIdx.x ].additional){return;}

//     src_pt_t <f_t> * pt_here;
//     int batchidx = local_info.batchidx;
//     // const int srcidx = blockIdx.x;
//     int idx;        // physical
//     bool ghost_direction;     // 1 means next, 0 means previous
//     idx = cross_info[ threadIdx.x ].idx_cross;
//     // ghost_direction = ( (idx<0) ? 0 : 1);


//     if((idx < batchidx*n_th))
//     {
//         local_info.shared_info->in_global = true;
//         return;
//     }

//     pt_here = &local_info.new_pts[ idx - batchidx*n_th ];

    
//     // pt_here = ( (idx >= batchidx*n_th) ? &local_info.new_pts[ idx - batchidx*n_th ] : &pool_margin[ srcidx * n_point_max + idx ] );


//     int j_test_p[2];
//     j_test_p[0] = cross_info[ threadIdx.x ].j_cross;
//     // j_test_p[1] = pt_here->next_j[j_test_p[0]];
//     j_test_p[1] = cross_info[ threadIdx.x ].j1;
//     ghost_direction = (j_test_p[0]>=3);
//     if(pt_here->images[j_test_p[0]].physical)
//     {
//         f_t deltaS, E123;
//         deltaS_error_cross(deltaS, E123, j_test_p[0], j_test_p[1], ghost_direction, pt_here);
//         pt_here->deltaS[j_test_p[0]]= deltaS;
//         pt_here->deltaS_Err[j_test_p[0]]= E123;
//     }
// }

__device__ void source_base_t::sum_area_0_g  ( local_info_t<f_t>& local_info ) const
{
    f_t* deltaS = (f_t*) local_info.deltaS_sum;
    f_t* Err = (f_t*) local_info.Err_sum;

    const int batchidx = local_info.batchidx;
    const int idx = threadIdx.x;
    const int i_src = blockIdx.x;
    int nextsrc_idx;
    int sum_length = 2;
    const int sum_size = (batchidx+1) * n_th;


    f_t temp_deltaS;
    f_t temp_Err;
    int pointidx;
    auto & ex_info = local_info.src_ext;
    auto & ret_info = local_info.src_ret;
    if(local_info.shared_info->Break)
    {

        ex_info.SolveSucceed = false;
        // if(!muted)
        // {
        //     printf("srcidx: %d: fail, break\n",srcidx);
        // }
        
        // ret_info.mag = -1;
        // ret_info.err = -1;
        // ex_info.Area = -0;
        // ex_info.Err_A = -0;

        return;
    }

    if(ex_info.SolveSucceed == false){

        deltaS[idx] = 0;
        Err[idx] = 0;

        // new update

        pointidx = local_info.neighbor3[ 0 ]; 
        if(pointidx < batchidx*n_th)
        {
            for(int j=3;j<5;j++)
            {
                pool_margin[i_src * n_point_max + pointidx].deltaS_Err[j] = local_info.pt_prev.deltaS_Err[j];
                pool_margin[i_src * n_point_max + pointidx].deltaS[j] = local_info.pt_prev.deltaS[j];
            }            
        }
        pointidx = local_info.neighbor3[ 1 ];
        for(int j=0;j<5;j++)
        {
            pool_margin[i_src * n_point_max + pointidx].deltaS_Err[j] = local_info.new_pts[ pointidx - batchidx*n_th ].deltaS_Err[j];
            pool_margin[i_src * n_point_max + pointidx].deltaS[j] = local_info.new_pts[ pointidx - batchidx*n_th ].deltaS[j];
        }
        pointidx = local_info.neighbor3[ 2 ];
        if(pointidx<batchidx*n_th)
        {
            for(int j=0;j<3;j++)
            {
                pool_margin[i_src * n_point_max + pointidx].deltaS_Err[j] = local_info.pt_next.deltaS_Err[j];
                pool_margin[i_src * n_point_max + pointidx].deltaS[j] = local_info.pt_next.deltaS[j];
            }
        }
   


        __syncthreads();


        for(int i=0;i<(batchidx+1);i++)        // ceil(NPOINTS / n_th)
        {
            pointidx = i*n_th + idx;
            if(pointidx >= sum_size){continue;}

            if(!pool_margin[i_src * n_point_max + pointidx].skip)
            {
                temp_deltaS = 0;
                temp_Err = 0;
                if(pointidx < batchidx*n_th)
                {
                    nextsrc_idx = next_src_idx_g( pointidx , i_src );
                }
                else
                {
                    nextsrc_idx = local_info.neighbor3[2];
                }
                
                src_pt_t < f_t > * pt_pointidx;
                src_pt_t < f_t > * pt_nextsrc_idx;
                pt_pointidx = ( (pointidx < batchidx*n_th) ? &pool_margin[i_src * n_point_max + pointidx] : &local_info.new_pts[ pointidx - batchidx*n_th ] );
                pt_nextsrc_idx = ( (nextsrc_idx < batchidx*n_th) ? &pool_margin[i_src * n_point_max + nextsrc_idx] : &local_info.new_pts[ nextsrc_idx - batchidx*n_th ] );
            
                // pt_pointidx = &pool_margin[i_src * n_point_max + pointidx];
                // pt_nextsrc_idx = &pool_margin[i_src * n_point_max + nextsrc_idx];

                for(int j=3;j<5;j++)
                {
                    temp_deltaS += pt_pointidx->deltaS[j];
                    temp_Err += pt_pointidx->deltaS_Err[j];
                }

                for(int j=0;j<3;j++)
                {
                    temp_deltaS += pt_nextsrc_idx->deltaS[j];
                    temp_Err += pt_nextsrc_idx->deltaS_Err[j];
                }
                pool_margin[i_src * n_point_max + pointidx].error_interval = temp_Err;
                pool_margin[i_src * n_point_max + pointidx].area_interval = temp_deltaS;
                deltaS[idx] += temp_deltaS;
                Err[idx] += temp_Err;
            }
            else
            {
                pool_margin[i_src * n_point_max + pointidx].error_interval = 0;
                pool_margin[i_src * n_point_max + pointidx].area_interval = 0;
            }


        }


        __syncthreads();



        while(sum_length<(n_th*2))
        {
            if(((idx & (sum_length-1))==0) && ((idx + sum_length/2) < n_th))    // sum_length = 2^n, idx%sum_length = idx & (sum_length-1)
            {
                deltaS[idx] += deltaS[idx + sum_length/2];
                Err[idx] += Err[idx + sum_length/2];
            }

            sum_length *= 2;
            __syncthreads();
        }
        // if(idx == 0)
        // {
        ret_info.mag = fabs(deltaS[0]) / local_info.src_shape.src_area;
        ret_info.err = Err[0] / local_info.src_shape.src_area;
        ex_info.Area = deltaS[0];
        ex_info.Err_A = Err[0];

        if(Err[0]/fabs(deltaS[0])<=RelTol+RelErrAllowed)
        {
            ex_info.SolveSucceed = true;
        }

        // }
        
    }
}

// __device__ void source_base_t::sum_area_0_g_test  ( local_info_t<f_t>& local_info ) const
// {
//     f_t* deltaS = (f_t*) local_info.deltaS_sum;
//     f_t* Err = (f_t*) local_info.Err_sum;

//     const int batchidx = local_info.batchidx;
//     const int idx = threadIdx.x;
//     const int i_src = blockIdx.x;
//     int nextsrc_idx;
//     int sum_length = 2;
//     const int sum_size = (batchidx+1) * n_th;


//     f_t temp_deltaS;
//     f_t temp_Err;
//     int pointidx;
//     auto & ex_info = local_info.src_ext;
//     auto & ret_info = local_info.src_ret;
//     if(local_info.shared_info->Break)
//     {

//         ex_info.SolveSucceed = false;
//         // if(!muted)
//         // {
//         //     printf("srcidx: %d: fail, break\n",srcidx);
//         // }
        
//         // ret_info.mag = -1;
//         // ret_info.err = -1;
//         // ex_info.Area = -0;
//         // ex_info.Err_A = -0;

//         return;
//     }

//     if(ex_info.SolveSucceed == false){

//         deltaS[idx] = 0;
//         Err[idx] = 0;

//         // new update

//         pointidx = local_info.neighbor3[ 0 ]; 
//         if(pointidx < batchidx*n_th)
//         {
//             for(int j=3;j<5;j++)
//             {
//                 pool_margin[i_src * n_point_max + pointidx].deltaS_Err[j] = local_info.pt_prev.deltaS_Err[j];
//                 pool_margin[i_src * n_point_max + pointidx].deltaS[j] = local_info.pt_prev.deltaS[j];
//             }            
//         }
//         pointidx = local_info.neighbor3[ 1 ];
//         if(pointidx < batchidx*n_th)
//         {
//             printf("srcidx: %d, pointidx: %d\n",blockIdx.x, threadIdx.x);
//         }
//         // for(int j=0;j<5;j++)
//         // {
//         //     pool_margin[i_src * n_point_max + pointidx].deltaS_Err[j] = local_info.new_pts[ pointidx - batchidx*n_th ].deltaS_Err[j];
//         //     pool_margin[i_src * n_point_max + pointidx].deltaS[j] = local_info.new_pts[ pointidx - batchidx*n_th ].deltaS[j];
//         // }
//         pointidx = local_info.neighbor3[ 2 ];
//         if(pointidx<batchidx*n_th)
//         {
//             for(int j=0;j<3;j++)
//             {
//                 pool_margin[i_src * n_point_max + pointidx].deltaS_Err[j] = local_info.pt_next.deltaS_Err[j];
//                 pool_margin[i_src * n_point_max + pointidx].deltaS[j] = local_info.pt_next.deltaS[j];
//             }
//         }
   


//         // __syncthreads();


//         // for(int i=0;i<(batchidx+1);i++)        // ceil(NPOINTS / n_th)
//         // {
//         //     pointidx = i*n_th + idx;
//         //     if(pointidx >= sum_size){continue;}

//         //     if(!pool_margin[i_src * n_point_max + pointidx].skip)
//         //     {
//         //         temp_deltaS = 0;
//         //         temp_Err = 0;
//         //         if(pointidx < batchidx*n_th)
//         //         {
//         //             nextsrc_idx = next_src_idx_g( pointidx , i_src );
//         //         }
//         //         else
//         //         {
//         //             nextsrc_idx = local_info.neighbor3[2];
//         //         }
                
//         //         src_pt_t < f_t > * pt_pointidx;
//         //         src_pt_t < f_t > * pt_nextsrc_idx;
//         //         pt_pointidx = ( (pointidx < batchidx*n_th) ? &pool_margin[i_src * n_point_max + pointidx] : &local_info.new_pts[ pointidx - batchidx*n_th ] );
//         //         pt_nextsrc_idx = ( (nextsrc_idx < batchidx*n_th) ? &pool_margin[i_src * n_point_max + nextsrc_idx] : &local_info.new_pts[ nextsrc_idx - batchidx*n_th ] );
            
//         //         // pt_pointidx = &pool_margin[i_src * n_point_max + pointidx];
//         //         // pt_nextsrc_idx = &pool_margin[i_src * n_point_max + nextsrc_idx];

//         //         for(int j=3;j<5;j++)
//         //         {
//         //             temp_deltaS += pt_pointidx->deltaS[j];
//         //             temp_Err += pt_pointidx->deltaS_Err[j];
//         //         }

//         //         for(int j=0;j<3;j++)
//         //         {
//         //             temp_deltaS += pt_nextsrc_idx->deltaS[j];
//         //             temp_Err += pt_nextsrc_idx->deltaS_Err[j];
//         //         }
//         //         pool_margin[i_src * n_point_max + pointidx].error_interval = temp_Err;
//         //         pool_margin[i_src * n_point_max + pointidx].area_interval = temp_deltaS;
//         //         deltaS[idx] += temp_deltaS;
//         //         Err[idx] += temp_Err;
//         //     }
//         //     else
//         //     {
//         //         pool_margin[i_src * n_point_max + pointidx].error_interval = 0;
//         //         pool_margin[i_src * n_point_max + pointidx].area_interval = 0;
//         //     }


//         // }


//         // __syncthreads();



//         // while(sum_length<(n_th*2))
//         // {
//         //     if(((idx & (sum_length-1))==0) && ((idx + sum_length/2) < n_th))    // sum_length = 2^n, idx%sum_length = idx & (sum_length-1)
//         //     {
//         //         deltaS[idx] += deltaS[idx + sum_length/2];
//         //         Err[idx] += Err[idx + sum_length/2];
//         //     }

//         //     sum_length *= 2;
//         //     __syncthreads();
//         // }
//         // // if(idx == 0)
//         // // {
//         // ret_info.mag = fabs(deltaS[0]) / local_info.src_shape.src_area;
//         // ret_info.err = Err[0] / local_info.src_shape.src_area;
//         // ex_info.Area = deltaS[0];
//         // ex_info.Err_A = Err[0];

//         // if(Err[0]/fabs(deltaS[0])<=RelTol+RelErrAllowed)
//         // {
//         //     ex_info.SolveSucceed = true;
//         // }

//         // }
        
//     }
// }

__device__ void source_base_t::sum_area_3_local ( local_info_t<f_t>& local_info ) const
{
    auto & ex_info = local_info.src_ext;
    auto & ret_info = local_info.src_ret;

    if((local_info.shared_info->Break))
    {
        ex_info.SolveSucceed = false;
        // ret_info.mag = -1;
        // ret_info.err = -1;
        // ex_info.Area = -0;
        // ex_info.Err_A = -0;
        return;
    }  

    const int batchidx = local_info.batchidx;
    if(ex_info.SolveSucceed){return;}
    int sum_length = 2;

    f_t Area,Error;
    f_t tempS;

    f_t* deltaS = (f_t*) local_info.deltaS_sum;
    f_t* Err = (f_t*) local_info.Err_sum;


    for(int loop2=0;loop2<2;loop2++)
    {
        deltaS[ threadIdx.x + loop2*n_th] = 0;
        Err[ threadIdx.x + loop2*n_th] = 0;            
    }

    if(!local_info.new_pts[ threadIdx.x ].skip)
    {
        int idx_list[3]; 

        idx_list[1] = local_info.neighbor3[1];
        idx_list[0] = local_info.neighbor3[0];
        idx_list[2] = local_info.neighbor3[2];

        src_pt_t<f_t>* pt_here = &local_info.new_pts[ threadIdx.x ];
        src_pt_t<f_t>* pt_prev = ((idx_list[0]<batchidx*n_th) ? &local_info.pt_prev : &local_info.new_pts[ idx_list[0] - batchidx*n_th ]);     // 这里必须用 shared 里面的新点，因为新点的面积只有 shared 里面算了
        src_pt_t<f_t>* pt_next = ((idx_list[2]<batchidx*n_th) ? &local_info.pt_next : &local_info.new_pts[ idx_list[2] - batchidx*n_th ]);     // 而且这里只读取一次，不需要额外拷贝了

        // loop2 = 0 /////////////////////////////
        f_t old_err = pt_prev->error_interval;
        f_t old_area = pt_prev->area_interval;
        if(idx_list[0]< batchidx*n_th)
        {
            for(int j=3;j<5;j++)
            {
                deltaS[ threadIdx.x ] += pt_prev->deltaS[j];
                Err[ threadIdx.x ] += pt_prev->deltaS_Err[j];
            }
            for(int j=0;j<3;j++)
            {
                deltaS[ threadIdx.x ] += pt_here->deltaS[j];
                Err[ threadIdx.x ] += pt_here->deltaS_Err[j]; 
            }
            pt_prev->error_interval = Err[ threadIdx.x ];
            pt_prev->area_interval = deltaS[ threadIdx.x ];
            Err[ threadIdx.x ] -= old_err;
            deltaS[ threadIdx.x ] -= old_area;
        }
        // loop2 = 1 /////////////////////////////
        for(int j=3;j<5;j++)
        {
            deltaS[ threadIdx.x + n_th] += pt_here->deltaS[j];
            Err[ threadIdx.x + n_th] += pt_here->deltaS_Err[j];
        }
        for(int j=0;j<3;j++)
        {
            deltaS[ threadIdx.x + n_th] += pt_next->deltaS[j];
            Err[ threadIdx.x + n_th] += pt_next->deltaS_Err[j];    
        }
        pt_here->error_interval = Err[ threadIdx.x + n_th ];
        pt_here->area_interval = deltaS[ threadIdx.x + n_th ];
    }

    // 后半加到前面  
    deltaS[ threadIdx.x ] += deltaS[ threadIdx.x + n_th ];
    Err[ threadIdx.x ] += Err[ threadIdx.x + n_th ];
    __syncthreads();
 

    while(sum_length<(2*n_th))        // 只对前 n_th求和，因为后面的已经加过来了
    {
        if((( threadIdx.x & (sum_length-1))==0) && (( threadIdx.x + sum_length/2) < n_th))    // sum_length = 2^n, idx%sum_length = idx & (sum_length-1)
        {
            deltaS[ threadIdx.x ] += deltaS[ threadIdx.x + sum_length/2];
            Err[ threadIdx.x ] += Err[ threadIdx.x + sum_length/2];

        }
        sum_length *= 2;
        __syncthreads();
    }


    // 对 global 修改 cross 并产生面积变化时，更改这里：
    deltaS[0] += local_info.shared_info->deltaS_cross_global;
    Err[0] += local_info.shared_info->Err_cross_global;
    __syncthreads();


    if(abs(deltaS[0]/ex_info.Area) < RelTol*1e-2 )
    {
        ex_info.SolveSucceed = true;
    }
    Area = deltaS[0] + ex_info.Area;
    Error = Err[0] + ex_info.Err_A;
    tempS = local_info.src_shape.src_area;

    ret_info.mag = fabs(Area) / tempS;
    ret_info.err = Error / tempS;
    ex_info.Area = Area;
    ex_info.Err_A = Error;

    if(Error/fabs(Area)<=RelTol+RelErrAllowed)
    {
        ex_info.SolveSucceed = true;
    } 
}

__device__ void source_base_t::adap_set_g ( local_info_t<f_t>& local_info ) const
{
    const int srcidx = blockIdx.x;
    int* ParentIdx = local_info.parent_idx;
    f_t* ErrorWeight = local_info.adap_sum;                  // size = n_th*(batchidx)
    int batchidx = local_info.batchidx;
    int idx = threadIdx.x;

    for(int i=0;i < batchidx;i++)        // ceil(NPOINTS / n_th)
    {
        ErrorWeight[i*n_th + idx] = 0;
        // 如果该点求解失败，那就没它的事，权重为0（相当于这个点不存在，将由附近的点负责查找）（而且它也压根没算面积误差）
        // 权重是“当前源点和下个源点之间的所有像面积之和”
        if(! pool_margin[ srcidx * n_point_max + i*n_th+idx ].skip)
        {
            ErrorWeight[(i*n_th + idx)] = pool_margin[ srcidx * n_point_max + i*n_th+idx ].error_interval;
            ErrorWeight[(i*n_th + idx)] += pool_margin[ srcidx * n_point_max + i*n_th+idx ].special_Err;
            pool_margin[ srcidx * n_point_max + i*n_th+idx ].special_Err = 0.;  

            ErrorWeight[(i*n_th + idx)] /= fabs(local_info.src_ext.Area);          // 使用相对误差

            ErrorWeight[(i*n_th + idx)] = max_f(0,ErrorWeight[(i*n_th + idx)] - RelTol/(batchidx*n_th));        // 如果自己的部分误差小于均值了，那算作过关       
        }
    }
    // input to shared memory finished
    __syncthreads();

    // 多段求和，每段用二分，再累加
    // 单段长度：n_th; 总长度：n_th*(batchidx)
    for(int i=0;i < batchidx;i++)        // ceil(batchidx*n_th / n_th)
    {
        int piece_length = 2;
        while(piece_length < 2*n_th)
        {
            if((((i*n_th + idx) & (piece_length-1)) >= piece_length/2) && ((i*n_th + idx) < n_th*batchidx))      // piece_length = 2^n, then idx % piece_length = idx & (piece_length-1)
            {
                ErrorWeight[(i*n_th + idx)] += ErrorWeight[(i*n_th + idx) - ((i*n_th + idx) & (piece_length-1)) + piece_length/2 -1];
            }
            __syncthreads();
            piece_length *= 2;
        }
    }
    // 每一段加前一段的尾巴，跳过第0段，从第1段开始，
    for(int i=1;i < batchidx;i++)        // ceil(batchidx*n_th / n_th)
    {
        ErrorWeight[(i*n_th + idx)] += ErrorWeight[(i*n_th - 1)];
        __syncthreads();
    }

    // 归一化
    f_t normalize = ErrorWeight[n_th*batchidx - 1];
    __syncthreads();
    // 如果noramlize比较大，归一化
    if(normalize >= RelErrAllowed)        // 与idx无关的分支，可以
    {
        for(int i=0;i < batchidx;i++)        // ceil(batchidx*n_th / n_th)
        {
            ErrorWeight[(i*n_th + idx)] /= normalize;
        }    
    }
    // 如果normalize非常小，即面积误差累加比阈值只超出非常小，为了规避数值误差（乃至NaN），
    // 这种情况换成均匀分布！由successcheck判定是否要break
    else
    {
        for(int i=0;i<batchidx;i++)        // ceil(batchidx*n_th / n_th)
        {
            ErrorWeight[(i*n_th + idx)] = f_t(i*n_th + idx + 1) / f_t(n_th*batchidx);
        }              
    }     
    __syncthreads();

    // ErrorWeight[i]意味着内存地址 i 的点和内存地址 next_i（不是i+1）单点之间的权重
    // 我只要连接left
    // 反向插值回去，前64个线程负责，以上述权重找64个新点
    // +0.5 避免总是插值到第0个点上
    f_t quantile = (idx+0.5) / f_t(n_th);
    int left_idx = bisearch_left(ErrorWeight,quantile,batchidx*n_th);      // 取值范围：[0,len-1]，即parent节点，相邻的left_idx在内存中连续，但是Q不连续
    int right_idx = next_src_idx_g(left_idx, srcidx);
    // 这里left_idx right_idx存的是左右points的指标
    // quantile_left的指标在points指标（即Q指标）-1，ErrW[0]对应pts[1]的分度值
    f_t left_Q = pool_margin[ srcidx * n_point_max + left_idx ].Q;
    f_t quantile_left = ( (left_idx == 0) ? 0 : ErrorWeight[left_idx-1]);
    f_t right_Q = pool_margin[ srcidx * n_point_max + right_idx ].Q;
    if(right_Q==0)
        right_Q = 1;
    // quantile_right = ErrorWeight[left_idx -1 +1];       // 注意quantile时right和left的差别
    // 不创建变量（塞不下了），直接调用共享内存的数据放在缓存里
    // 线性插值
    f_t Q = left_Q + (right_Q - left_Q) * (quantile-quantile_left)/(ErrorWeight[left_idx -1 +1]-quantile_left);         // 不创建变量（塞不下了），quantile_right = ErrorWeight[left_idx -1 +1]; 
    __syncthreads();

    // 在此往下，不再使用 ErrorWeight，需要覆写原来 ErrorWeight 的部分，因而需要同步
    // pool_margin[ srcidx * n_point_max + idx + batchidx*n_th].Q = Q;
    local_info.new_pts[idx].Q = Q;
    c_t marginloc;
    margin_Q(marginloc, local_info.src_shape, Q);
    local_info.new_pts[idx].position = marginloc;
    local_info.new_pts[idx].Nphys = -1;
    local_info.new_pts[idx].special_Err = 0;
    for(int j=0;j<order-1;j++)
    {
        local_info.new_pts[idx].deltaS[j]=0;
        local_info.new_pts[idx].deltaS_Err[j]=0;
    }

    // 别急，下面还要连接顺序
    ParentIdx[idx] = left_idx;
    __syncthreads();

    // 前半截是“谁连自己”
    // 如果前一个点和自己是姊妹节点
    if((idx >= 1) && (ParentIdx[idx-1]==ParentIdx[idx]))      // 0号的前一个要取模，但负数取模有一定问题，而且反正0号线程肯定连parent
    {
        local_info.new_pts[ idx-1 ].next_src_idx = idx + batchidx*n_th;     // 前一个点后面连自己
        local_info.new_pts[  idx  ].prev_src_idx = idx-1 + batchidx*n_th;
    }
    else
    {
        pool_margin[ srcidx * n_point_max + left_idx].next_src_idx = idx + batchidx*n_th;               // parent节点后面连自己
        local_info.new_pts[idx].prev_src_idx = left_idx;
    }
    // 后半截是“自己连谁”
    // 如果自己是姊妹节点中的老幺
    if(idx==(n_th-1))     // 如果是本batch中的最后一个
    {
        local_info.new_pts[idx].next_src_idx = right_idx;                   // 自己连 right idx
        pool_margin[ srcidx * n_point_max + right_idx].prev_src_idx = idx + batchidx*n_th;
    }
    else{
        if(ParentIdx[idx] != ParentIdx[idx+1])  // 如果是普通情况
        {
            local_info.new_pts[idx].next_src_idx = right_idx;                   // 自己连 right idx
            pool_margin[ srcidx * n_point_max + right_idx].prev_src_idx = idx + batchidx*n_th;
        }
        // else{}   // 前面连好了已经
    }
}

////////////////////////////////////////////////////////////
//

__host__ void source_base_t::run( device_t & f_dev )
{
    auto n_bl = ( n_src + n_th - 1 ) / n_th;
    auto s_sh = 0;


    auto lpar = std::make_tuple
        ( dim3( n_bl ), dim3( n_th ), s_sh );

    f_dev.launch( point_approximation, lpar, stream, ( * this ) );

    n_bl = n_src;
    s_sh = sizeof(int)*n_th + 
    std::max( sizeof(float2_t)*n_point_max, 
              sizeof(float2_t)   * 8*n_th
            + sizeof(shared_info_t< float2_t >) * 1 
            + sizeof(cross_info_t) * n_cross_max * 4
            + sizeof(src_pt_t< float2_t >) * n_th);
    // printf("shared memory size: %d\n",s_sh);
    lpar = std::make_tuple
        ( dim3( n_bl ), dim3( n_th ), s_sh );
    f_dev.launch( solve_extended_sh, lpar, stream, ( * this ) );


    // n_bl = n_src;
    // lpar = std::make_tuple
    //     ( dim3( n_bl ), dim3( n_th ), s_sh );
    // f_dev.launch( boundaryset_solve_connect, lpar, stream, ( * this ) );

    // n_bl = ( n_src + n_th - 1 ) / n_th;
    // lpar = std::make_tuple
    //     ( dim3( n_bl ), dim3( n_th ), s_sh );
    // f_dev.launch( slope_detector, lpar, stream, ( * this ) );

    // n_bl = n_src;
    // lpar = std::make_tuple
    //     ( dim3( n_bl ), dim3( n_th ), s_sh );
    // f_dev.launch( area_normal, lpar, stream, ( * this ) );

    // n_bl = n_src; 
    // // s_sh = 2*sizeof(f_t)* n_th;
    // lpar = std::make_tuple
    //     ( dim3( n_bl ), dim3( n_th ), s_sh );
    // f_dev.launch( sum_area, lpar, stream, ( * this ) );


    return;
}

};                              // namespace twinkle
