#include "source_base.h"
#include "../utils/solve.h"

namespace twinkle
{

using f_t = source_base_t::f_t;
using c_t = source_base_t::c_t;

__device__ void source_base_t::margin_Q
( c_t & marginloc, const src_shape_t<f_t>& shape, const float& Q ) const
{
    f_t theta = Q*2*PI;
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
    if(astrom)
    {
        src.astromX_interval = 0.;
    }
    return;
}


__device__ void source_base_t::dz_wedge
( c_t & dz, f_t & wedge, f_t & jacobian, const img_t & img, const c_t & zeta, const c_t & loc_centre, const f_t rho, const f_t lens_s) const
{
	c_t pzeta_pconjz,ppzeta_ppconjz;
	c_t items[2];
    c_t zc;
    zc.set(img.position.re, -img.position.im);

	items[0] = m1 / (( zc + m2 * lens_s ) * ( zc + m2 * lens_s ));
	items[1] = m2 / (( zc - m1 * lens_s ) * ( zc - m1 * lens_s ));
	pzeta_pconjz = items[0] + items[1];

	jacobian = 1 - pzeta_pconjz.norm2(  );

	items[0] = -2 * items[0] / (zc + m2 * lens_s );
	items[1] = -2 * items[1] / (zc - m1 * lens_s );
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
( img_t * imgs, c_t * dz, f_t * wedge, int & phys_tst, int * temp_j_from_out_j, const c_t & zeta, const c_t & loc_centre, const f_t Rho, const f_t lens_s) const
{
    bool skip = solve_imgs(imgs, zeta, true, lens_s);
    if(skip){ return true; }

    f_t temp_jacobi[ 5 ];
    for(int j=0;j<order-1;j++)
    {   
        dz_wedge( dz[j], wedge[j], temp_jacobi[j], imgs[j], zeta, loc_centre, Rho, lens_s );
        imgs[j].parity = ( (temp_jacobi[j]>0) ? 0 : 1 );
    }   
    is_physical_all( imgs, zeta, phys_tst, lens_s ); 
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

}; // namespace twinkle
