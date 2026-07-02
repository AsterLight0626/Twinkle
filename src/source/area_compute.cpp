#include "source_base.h"

namespace twinkle
{

using f_t = source_base_t::f_t;
using c_t = source_base_t::c_t;

__device__ f_t source_base_t::delta_s_1
( const c_t & z_here, const c_t & z_next) const
{
    return -((z_next.im + z_here.im) * (z_next.re - z_here.re) - (z_next.re + z_here.re) * (z_next.im - z_here.im)) / 4;
}

__device__ void source_base_t::delta_X_1
( f_t& real, f_t&imag, const c_t & z_here, const c_t & z_next) const
{
    // f_t real, imag;
    real = 1./8. * (z_next.re+z_here.re)*(z_next.re+z_here.re)*(z_next.im-z_here.im);
    imag = -1./8. * (z_next.im+z_here.im)*(z_next.im+z_here.im)*(z_next.re-z_here.re);
    return;
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
    if(next_j < 0 || next_j >= 5)
    {
        deltaS = 0;
        Error1234 = 0;
        return;
    }
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

__device__ void source_base_t::deltaS_error_image_beta
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
    if(next_j < 0 || next_j >= 5)
    {
        deltaS = 0;
        Error1234 = 0;
        return;
    }
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
    E_1 = fabs((wedge0 - wedge1) * theta3) / 48;
    deter = theta2 * (dzs0*dzs1).abs(  );
    if(fabs(deter)>2e-12){E_2 = 1.5 * fabs( deltaS_p *  ((zs0-zs1).norm2(  ) / deter -1));}
    else
    {
        E_2=0;
    }        // usually caused by numericial error
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
const src_pt_t < f_t > * pt_here, const src_pt_t < f_t > * pt_other, 
c_t * astromX_out)    const
{
    f_t Qother_src, Qhere;
    f_t theta1, theta2, theta3;
    int next_i;
    f_t deltaS_tp, E_1234;
    // astrometry
    c_t delta_astromX;

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
            if(astrom)
            {
                astromX_out[jghost] = 0.;
            }
            continue;
        }
        next_i = pt_here->next_idx[j];
        if(here_idx==next_i)
        {
            deltaS_error_cross_beta(deltaS_tp, E_1234, j, pt_here->next_j[j], (j>=3), pt_here);      // ghost_direction = bool (j>=3)
            deltaS_out[j] = deltaS_tp;
            deltaS_Err_out[j] = E_1234;
            // // ！！需要补充 cross 的 astromX
            // printf("cross! idx: %d, j: %d\n", threadIdx.x, j);
            if(astrom)
            {
                deltaX_cross( delta_astromX,  j, pt_here->next_j[j], (j>=3), pt_here);   
                astromX_out[j] = delta_astromX;
            }
            continue;
        }
        deltaS_error_image_beta( deltaS_tp, E_1234, j, pt_here, pt_other, theta1, theta2, theta3);

        deltaS_out[j] = deltaS_tp;
        deltaS_Err_out[j] = E_1234; \
        if(astrom)
        {
            deltaX_normal( delta_astromX,  j, pt_here, pt_other, theta3);   
            astromX_out[j] = delta_astromX;
        }

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

__device__ void source_base_t::deltaS_error_cross_beta
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
    f_t deltaS_p1 = (wedge[plus] - wedge[minus]) * theta3 /24;
    f_t deltaS_p2 = wedge_product((zs[1]-zs[0]) , (dzs[1]+dzs[0])) * theta1 / 12;
    f_t deltaS_p = (deltaS_p1 + deltaS_p2) / 2;
    f_t E_1 = fabs(wedge[plus] + wedge[minus]) * theta3 / 48;
    c_t temp0 = (zs[plus]-zs[minus])*(dzs[plus]-dzs[minus]);
    f_t temp1 = 2* (zs[plus]-zs[minus]).abs(  ) * sqrt((dzs[plus]*dzs[minus]).abs(  ));
    temp0 = ( (ghost_direction==1) ? temp0 + temp1 : temp0 - temp1);
    f_t E_2 = 1.5* (temp0).abs(  ) * theta1;
    f_t E_3 = 0.1 * fabs(deltaS_p) * theta2;
    f_t E_4 = fabs((deltaS_p1 - deltaS_p2));
    deltaS = deltaS_t+deltaS_p;
    Error1234 = E_1+E_2+E_3+E_4;
    return;
}


__device__ void source_base_t::deltaX_normal
( c_t & delta_astromX, const int j,
    const src_pt_t < f_t > * pt_here, const src_pt_t < f_t > * pt_other,
    const f_t & theta3) const
{
    // c_t delta_astromX_t, delta_astromX_p;

    f_t temp_real, temp_imag;

    c_t zs0,zs1,dzs0,dzs1;
    f_t wedge0,wedge1;
    int next_j;

    next_j = pt_here->next_j[j];
    if(next_j < 0 || next_j >= 5)
    {
        delta_astromX = 0;
        return;
    }
    zs0 = pt_here->images[j].position;
    zs1 = pt_other->images[next_j].position;
    wedge0 = pt_here->wedge[j];
    wedge1 = pt_other->wedge[next_j];
    dzs0 = pt_here->dz[j];
    dzs1 = pt_other->dz[next_j];

    delta_X_1(temp_real, temp_imag, zs0,zs1);
    c_t delta_astromX_t = c_t(temp_real, temp_imag);

    // avgwedge1 = (xi_1 * wedge_i + xi1_1 * wedge_i1) * mi
    // avgwedge2 = (xi_2 * wedge_i + xi1_2 * wedge_i1) * mi
    // dx2 = dz_i.imag + dz_i1.imag
    // d2x2 = dx2**2
    // dx1 = dz_i.real + dz_i1.real
    // d2x1 = dx1**2

    temp_real = -(-0.125 * (dzs0.re + dzs1.re)*(dzs0.re + dzs1.re) * (dzs0.im + dzs1.im)  - (zs0.re * wedge0 + zs1.re * wedge1) ) * theta3 / 24.;
    temp_imag =  (-0.125 * (dzs0.im + dzs1.im)*(dzs0.im + dzs1.im) * (dzs0.re + dzs1.re)  + (zs0.im * wedge0 + zs1.im * wedge1) ) * theta3 / 24.;

    // printf("threadIdx: %d, j: %d, (%.16f, %.16f), (%.16f, %.16f) \n",threadIdx.x, j, temp_real, temp_imag, delta_astromX_t.re,delta_astromX_t.im);
    // printf("theta3: %.16f, first part: (%.16f), second part: (%.16f)\n", theta3, ((dzs0.re + dzs1.re)*(dzs0.re + dzs1.re) * (dzs0.im + dzs1.im)), (zs0.re * wedge0 + zs1.re * wedge1) );

    c_t delta_astromX_p = c_t(temp_real, temp_imag);
    // delta_astromX_p = 0.;

    delta_astromX = delta_astromX_t + delta_astromX_p;

    return;
}

__device__ void source_base_t::deltaX_cross
( c_t & delta_astromX, const int j0, const int j1, const bool ghost_direction,
    const src_pt_t < f_t > * pt_here) const
{
    f_t temp_real, temp_imag;

    // int next_j;
    // next_j = pt_here->next_j[j];

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

    // f_t deltaS_t = delta_s_1(zs[0],zs[1]);
    f_t theta1 = (zs[1]-zs[0]).abs(  ) / sqrt((dzs[0]*dzs[1]).abs(  ));
    f_t theta2 = theta1*theta1;
    f_t theta3 = theta2*theta1;
    int plus = int(!ghost_direction);
    int minus = int(ghost_direction);

    delta_X_1(temp_real, temp_imag, zs[0],zs[1]);
    c_t delta_astromX_t = c_t(temp_real, temp_imag);

    // avgwedge1 = (xi_1 * wedge_i + xi1_1 * wedge_i1) * mi
    // avgwedge2 = (xi_2 * wedge_i + xi1_2 * wedge_i1) * mi
    // dx2 = dz_i.imag + dz_i1.imag
    // d2x2 = dx2**2
    // dx1 = dz_i.real + dz_i1.real
    // d2x1 = dx1**2

    // temp_real = -(-0.125 * (dzs[plus].re - dzs[minus].re)*(dzs[plus].re - dzs[minus].re) * (dzs[plus].im - dzs[minus].im)  - (zs[plus].re * wedge[plus] - zs[minus].re * wedge[minus]) ) * theta3 / 24.;
    // temp_imag =  (-0.125 * (dzs[plus].im - dzs[minus].im)*(dzs[plus].im - dzs[minus].im) * (dzs[plus].re - dzs[minus].re)  + (zs[plus].im * wedge[plus] - zs[minus].im * wedge[minus]) ) * theta3 / 24.;

    temp_real =  ((dzs[plus].re * dzs[plus].re * dzs[plus].im + zs[plus].re * wedge[plus]) - (dzs[minus].re * dzs[minus].re * dzs[minus].im + zs[minus].re * wedge[minus])) * theta3 / 24.;
    temp_imag = -((dzs[plus].im * dzs[plus].im * dzs[plus].re + zs[plus].im * wedge[plus]) - (dzs[minus].im * dzs[minus].im * dzs[minus].re + zs[minus].im * wedge[minus])) * theta3 / 24.;


    // // printf("threadIdx: %d, j: %d, (%.16f, %.16f), (%.16f, %.16f) \n",threadIdx.x, j, temp_real, temp_imag, delta_astromX_t.re,delta_astromX_t.im);
    // // printf("theta3: %.16f, first part: (%.16f), second part: (%.16f)\n", theta3, ((dzs0.re + dzs1.re)*(dzs0.re + dzs1.re) * (dzs0.im + dzs1.im)), (zs0.re * wedge0 + zs1.re * wedge1) );

    // f_t deltaS_p1 = (wedge[plus] - wedge[minus]) * theta3 /24;
    // f_t deltaS_p2 = wedge_product((zs[1]-zs[0]) , (dzs[1]+dzs[0])) * theta1 / 12;
    // f_t deltaS_p = (deltaS_p1 + deltaS_p2) / 2;

    c_t delta_astromX_p = c_t(temp_real, temp_imag);
    // c_t delta_astromX_p = 0.;

    // printf("j0: %d, j1: %d, zsplus, minus: (%.16f, %.16f), (%.16f, %.16f), delta_astromX_p: (%.16f, %.16f)\n", j0,j1,  zs[plus].re, zs[plus].im, zs[minus].re, zs[minus].im, delta_astromX_p.re/deltaS_p, delta_astromX_p.im/deltaS_p);
    // f_t part1 = -0.125 * (dzs[plus].re - dzs[minus].re)*(dzs[plus].re - dzs[minus].re) * (dzs[plus].im - dzs[minus].im)* theta3 / 24. / deltaS_p;
    // f_t part2 = (zs[plus].re * wedge[plus] - zs[minus].re * wedge[minus]) * theta3 / 24. / deltaS_p;
    // printf("++: %.16f, +-: %.16f, -+: %.16f, --: %.16f\n", part1+part2, part1-part2,-part1+part2, -part1-part2);

    delta_astromX = delta_astromX_t + delta_astromX_p;

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
const int jj, const int idx_here, bool& Break, int& Ncross_all, const int bcidx ) const
{
    Break = false;
    if(pt_here->next_idx[jj] != idx_here)
    {
        return false;
    }
    // else
    // {
    //     if(bcidx==0)
    //     {
    //         atomicAdd(&Ncross_all, 1);            
    //     }
    // }

    // N5↔N3 boundary self-connection from cross-candidate (negative_connect_j012),
    // not a real caustic crossing — skip slope test
    if((pt_here->Nphys == 5 && pt_other->Nphys == 3) ||
       (pt_here->Nphys == 3 && pt_other->Nphys == 5))
        return false;

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

}; // namespace twinkle
