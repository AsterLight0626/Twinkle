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

__device__ void source_base_t::solve_point_approx( f_t coeff_RelTol ) const
{
    const int i_src = threadIdx.x + n_th * blockIdx.x;
    if( i_src >= n_src )
        return;

    auto & src  = pool_center [ i_src ];
    const c_t zeta = src.loc_centre;
    auto  & mag_err = pool_mag [i_src];
    auto  & ext_info = pool_extended [i_src];
    
    c_t temp_astrom_X = 0.;
    c_t* ptr_astrom_Th = nullptr;
    if(astrom)
    {
        ptr_astrom_Th = &pool_astrom_Th [i_src];
    }
    const f_t Rho = src.rho;

    img_t temp_images[order-1];
    int phys_tst = -1;

    const f_t lens_s = pool_lens_s[ i_src ] ;

    if(solve_imgs(temp_images, zeta, false, lens_s , ext_info.roots_point)){return;}      // solve_imgs return "fail"

    is_physical_all( temp_images , zeta , phys_tst, lens_s );

    f_t muPS = 0, muQC = 0;
    for(int j=0;j<order-1;j++)
    {
        if(temp_images[j].physical)
        {
            f_t temp_pt_mag = 1/fabs(jacobian( temp_images[j].position, lens_s ));
            muPS += temp_pt_mag;
            muQC += fabs(mu_qc(temp_images[j].position, lens_s ) * (Rho+1e-3) * (Rho+1e-3));
            if(astrom)
            {
                temp_astrom_X += temp_images[j].position * temp_pt_mag;                
            }
        }
    }

    muQC *= C_Q;
    
    bool additional_test
    = ghost_test_all( temp_images, phys_tst, zeta, Rho, lens_s )
    && safe_distance_test(zeta, Rho, lens_s );

    bool is_point_src
     = (muQC < coeff_RelTol * RelTol * muPS) && additional_test;

    if( is_point_src )
    {
        ext_info.SolveSucceed = true;       // default: false
    }
    if( !additional_test )      // 如果没通过额外测试，也就是 ghost_test 或者 safe_distance 失败了
    {
        // muQC += coeff_RelTol * RelTol * muPS;     
        muQC += muPS;           // 多加些误差
    }
    mag_err.mag = muPS;
    mag_err.err = muQC;

    if(astrom)
    {
        *ptr_astrom_Th = temp_astrom_X / muPS;
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

    const f_t lens_s = local_info.lens_s;

    for(int j=0;j<order-1;j++)
    {
        temp_images[j].position = ex_info.roots_point[j];
    }

    //// calculation
    bool skip = margin_solve_cal (temp_images, dz, wedge, phys_tst, temp_j_from_order_j, zeta, loc_centre, Rho, lens_s );

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

    int next_idx_out[5] = {-1, -1, -1, -1, -1};
    int next_j_out[5] = {-1, -1, -1, -1, -1};
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


__device__ void source_base_t::slope_test_local  (local_info_t<f_t>& local_info, const int bcidx ) const
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
            if(slope_test(pt_prev, pt_here, jj, prev_idx, Break, local_info.shared_info->Ncross_all, bcidx))
            {
                // local_info.shared_info->Break = true;
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
        if(slope_test(pt_here, pt_prev, jj, idx, Break, local_info.shared_info->Ncross_all, bcidx))
        {
            // local_info.shared_info->Break = true;
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
        if(slope_test(pt_here, pt_next, jj, idx, Break, local_info.shared_info->Ncross_all, bcidx))
        {
            // local_info.shared_info->Break = true;
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
            if(slope_test(pt_next, pt_here, jj, next_idx, Break, local_info.shared_info->Ncross_all, bcidx))
            {
                // local_info.shared_info->Break = true;
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
        idx_cross_g = ( (j_test_p[0]<3) ? prev_src_idx_g(idx_cross_p, blockIdx.x) : next_src_idx_g(idx_cross_p, blockIdx.x) );
        // printf("idx_g: %d\n",idx_cross_g);

    
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
            if(pool_margin[srcidx * n_point_max + idx_cross_g].images[2].physical || pool_margin[srcidx * n_point_max + idx_cross_g].images[3].physical)
            {
                local_info.shared_info->Break = true;
                return;
            }

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
                if(j_test_p[1] < 0 || j_test_p[1] >= 5) {
                    local_info.shared_info->Break = true;
                    return;
                }
                if(ghost_direction)
                {
                    idx_cross_p = prev_src_idx_g(idx_cross_p, blockIdx.x);
                }
                else
                {       
                    idx_cross_p = next_src_idx_g(idx_cross_p, blockIdx.x);                  
                }
                // printf("idx_p: %d\n",idx_cross_p);
                if(pool_margin[srcidx * n_point_max + idx_cross_p].Nphys!=5)
                {
                    // 该片段已走到尽头，也没通过
                    test_pass = true;
                    additional_piece = true;
                    break;
                }
                bool j_found = false;
                for(int j=0;j<5;j++)        // 连接：因为nextj是单向的，有一个方向需要遍历查找
                {
                    if(pool_margin[srcidx * n_point_max + idx_cross_p].images[j].parity != ghost_direction)
                    {
                        temp_j = pool_margin[srcidx * n_point_max + idx_cross_p].next_j[j];
                        if((temp_j == j_test_p[0]))
                        {
                            j_test_p[0] = j;
                            j_found = true;
                            break;
                        }                                    
                    }

                }
                if(!j_found)
                {
                    test_pass = true;
                    additional_piece = true;
                    break;
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
                        // printf("idx_p: %d, j0: %d, j1: %d\n",idx_cross_p,j_test_p[0],j_test_p[1]);
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
                            prev_pt_idx = prev_src_idx_g(idx_change[0], blockIdx.x);
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
                            prev_pt_idx = next_src_idx_g(idx_change[0], blockIdx.x);   
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

                    if(j0[iii]!=2 && j0[iii]!=3)
                    {
                        local_info.shared_info->Break = true;
                        // printf("j0 not 2 or 3!\n");
                        return;
                    }
                    if(j1[iii]!=2 && j1[iii]!=3)
                    {
                        local_info.shared_info->Break = true;
                        // printf("j1 not 2 or 3!\n");
                        return;
                    }

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
                // printf("idx_change[0]: %d, j0[0]: %d, j1[0]: %d\n",idx_change[0], j0[0],j1[0]);

                idx_cross_g = (ghost_direction) ? next_src_idx_g(idx_cross_p, blockIdx.x) \
                                                : prev_src_idx_g(idx_cross_p, blockIdx.x) ;
                // printf("idx_p: %d, idx_g inv: %d\n",idx_cross_p,idx_cross_g);
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
                j_test_g[0] = 2;
                j_test_g[1] = 3;
                if(pool_margin[srcidx * n_point_max + idx_cross_g].images[2].physical || pool_margin[srcidx * n_point_max + idx_cross_g].images[3].physical)
                {
                    local_info.shared_info->Break = true;
                    return;
                }

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
                // printf("idx_g: %d, j0: %d, j1: %d\n",idx_cross_g, j_test_g[int(!tempbool)],j_test_g[int( tempbool)]);

                cos2 = (temp*temp) / (norms[iter_i-1]*norms[iter_i]);
                if(cos2<=0.1){test_pass = true;}
                else{
                    if(cos2<0.9)
                    {
                        // printf("unclear! srcidx: %d\n",blockIdx.x);
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
                    idx_cross_g = (ghost_direction) ? next_src_idx_g(idx_cross_p, blockIdx.x) \
                                                    : prev_src_idx_g(idx_cross_p, blockIdx.x) ;
                    // printf("idx_p: %d\n",idx_cross_p);
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
                                // printf("unclear! srcidx: %d\n",blockIdx.x);
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

    // astrometry
    c_t astromX_out[5];


    pt_here = &local_info.new_pts[ threadIdx.x ];
    pt_prev = &local_info.pt_prev;
    pt_next = &local_info.pt_next;


    if(local_info.new_pts[ threadIdx.x ].skip){ return; }

    // loop3 == 0   ///////////////////
    if(idx_prev<batchidx*n_th)
    {
        // positive
        parity = 0;
        deltaS_error_parity(deltaS_out, deltaS_Err_out, parity, idx_prev, idx_here, pt_prev, pt_here, astromX_out);

        for(int output_j=3; output_j<5;output_j++)
        {
            pt_prev->deltaS[output_j] = deltaS_out[output_j];
            pt_prev->deltaS_Err[output_j] = deltaS_Err_out[output_j];
            if(astrom)
            {
                pt_prev -> delta_astromX[output_j]= astromX_out[output_j];
            }
        }
    }

    // loop3 == 1       ////////////////////////

    // negative
    parity = 1;                      
    deltaS_error_parity(deltaS_out, deltaS_Err_out, parity, idx_here, idx_prev, pt_here, pt_prev, astromX_out);                    
    // positive
    parity = 0;
    deltaS_error_parity(deltaS_out, deltaS_Err_out, parity, idx_here, idx_next, pt_here, pt_next, astromX_out);

    for(int output_j=0; output_j<5;output_j++)
    {
        pt_here->deltaS[output_j] = deltaS_out[output_j];
        pt_here->deltaS_Err[output_j] = deltaS_Err_out[output_j];
        if(astrom)
        {
            pt_here -> delta_astromX[output_j]= astromX_out[output_j];
            // printf("threadIdx: %d, j: %d, delta_astromX: %.16f, %.16f\n",threadIdx.x, output_j, astromX_out[output_j].re,astromX_out[output_j].im);
            // printf("%d, %d, %.16f, %.16f\n",threadIdx.x, output_j, pt_here->images[output_j].position.re, pt_here->images[output_j].position.im);
        }
    }

    // loop3 == 2 /////////////////////
    if(idx_next<batchidx*n_th)
    {
        // negative
        parity = 1;                        
        deltaS_error_parity(deltaS_out, deltaS_Err_out, parity, idx_next, idx_here, pt_next, pt_here, astromX_out);

        for(int output_j=0; output_j<3;output_j++)
        {
            pt_next->deltaS[output_j] = deltaS_out[output_j];
            pt_next->deltaS_Err[output_j] = deltaS_Err_out[output_j];
            if(astrom)
            {
                pt_next -> delta_astromX[output_j]= astromX_out[output_j];
            }
        }
    }
}

__device__ void source_base_t::sum_area_0_g  ( local_info_t<f_t>& local_info, f_t coeff_RelTol ) const
{
    f_t* deltaS = (f_t*) local_info.deltaS_sum;
    f_t* Err = (f_t*) local_info.Err_sum;
    c_t* astromX = (c_t*) local_info.astromX_sum;

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
    // astrometry
    c_t temp_astromX;
    auto & astrom_Th_info = local_info.src_astrom_Th;
    
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
        if(astrom)
        {
            astromX[idx] = 0.;
        }

        // new update

        pointidx = local_info.neighbor3[ 0 ]; 
        if(pointidx < batchidx*n_th)
        {
            for(int j=3;j<5;j++)
            {
                pool_margin[i_src * n_point_max + pointidx].deltaS_Err[j] = local_info.pt_prev.deltaS_Err[j];
                pool_margin[i_src * n_point_max + pointidx].deltaS[j] = local_info.pt_prev.deltaS[j];
                if(astrom)
                {
                    pool_margin[i_src * n_point_max + pointidx].delta_astromX[j] = local_info.pt_prev.delta_astromX[j];
                }
            }            
        }
        pointidx = local_info.neighbor3[ 1 ];
        for(int j=0;j<5;j++)
        {
            pool_margin[i_src * n_point_max + pointidx].deltaS_Err[j] = local_info.new_pts[ pointidx - batchidx*n_th ].deltaS_Err[j];
            pool_margin[i_src * n_point_max + pointidx].deltaS[j] = local_info.new_pts[ pointidx - batchidx*n_th ].deltaS[j];
            if(astrom)
            {
                pool_margin[i_src * n_point_max + pointidx].delta_astromX[j] = local_info.new_pts[ pointidx - batchidx*n_th ].delta_astromX[j];
            }
        }
        pointidx = local_info.neighbor3[ 2 ];
        if(pointidx<batchidx*n_th)
        {
            for(int j=0;j<3;j++)
            {
                pool_margin[i_src * n_point_max + pointidx].deltaS_Err[j] = local_info.pt_next.deltaS_Err[j];
                pool_margin[i_src * n_point_max + pointidx].deltaS[j] = local_info.pt_next.deltaS[j];
                if(astrom)
                {
                    pool_margin[i_src * n_point_max + pointidx].delta_astromX[j] = local_info.pt_next.delta_astromX[j];
                }
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
                temp_astromX = 0.;
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
                    if(astrom)
                    {
                        temp_astromX += pt_pointidx->delta_astromX[j];
                    }
                }

                for(int j=0;j<3;j++)
                {
                    temp_deltaS += pt_nextsrc_idx->deltaS[j];
                    temp_Err += pt_nextsrc_idx->deltaS_Err[j];
                    if(astrom)
                    {
                        temp_astromX += pt_nextsrc_idx->delta_astromX[j];
                    }
                }
                pool_margin[i_src * n_point_max + pointidx].error_interval = temp_Err;
                pool_margin[i_src * n_point_max + pointidx].area_interval = temp_deltaS;
                deltaS[idx] += temp_deltaS;
                Err[idx] += temp_Err;
                if(astrom)
                {
                    pool_margin[i_src * n_point_max + pointidx].astromX_interval = temp_astromX;
                    astromX[idx] += temp_astromX;
                    // printf("idx_here: %d, idx_next: %d, astromX, %.16f, %.16f\n", pointidx, nextsrc_idx, temp_astromX.re, temp_astromX.im);
                    // printf("%d, %d, %.16f, %.16f\n", pointidx, nextsrc_idx, temp_astromX.re, temp_astromX.im);
                }
            }
            else
            {
                pool_margin[i_src * n_point_max + pointidx].error_interval = 0;
                pool_margin[i_src * n_point_max + pointidx].area_interval = 0;
                if(astrom)
                {
                    pool_margin[i_src * n_point_max + pointidx].astromX_interval =0.;
                }
            }


        }




        __syncthreads();



        while(sum_length<(n_th*2))
        {
            if(((idx & (sum_length-1))==0) && ((idx + sum_length/2) < n_th))    // sum_length = 2^n, idx%sum_length = idx & (sum_length-1)
            {
                deltaS[idx] += deltaS[idx + sum_length/2];
                Err[idx] += Err[idx + sum_length/2];
                if(astrom)
                {
                    astromX[idx] += astromX[idx + sum_length/2];
                }
            }

            sum_length *= 2;
            __syncthreads();
        }
        // if(idx == 0)
        // {

        // ret_info.mag = fabs(deltaS[0]) / local_info.src_shape.src_area;
        // ret_info.err = Err[0] / local_info.src_shape.src_area;


        if( ((fabs(deltaS[0]) - 3*Err[0]) / local_info.src_shape.src_area) <= sqrt(1+(2/local_info.src_shape.rho)*(2/local_info.src_shape.rho)))
        {
            ret_info.mag = fabs(deltaS[0]) / local_info.src_shape.src_area;
            ret_info.err = Err[0] / local_info.src_shape.src_area;
            if(astrom)
            {
                // printf("astromX[0]: %.16f, %.16f, Area: %.16f\n",astromX[0].re,astromX[0].im, fabs(deltaS[0]));
                astrom_Th_info = astromX[0] / fabs(deltaS[0]);
                // if(threadIdx.x==0)
                // {
                //     printf("astromX: %.16f, %.16f\n", astromX[0].re, astromX[0].im);
                // }
            }
        }
        else        // 如果放大率高于理论上界，且高出的值大于三倍误差，认为此时求解出错，取上一次的结果
        {
            ret_info.mag = min(ret_info.mag, sqrt(1+(2/local_info.src_shape.rho)*(2/local_info.src_shape.rho)));
            ret_info.err = -1;
            ex_info.SolveSucceed = true;
            ex_info.Break = true;               // 应该用 Break，但还在测试，所以吧 Break 和 SolveSucceed 都标记了
            // astrom_Th_info = astrom_Th_info;        // 此时不修改 astromX_info，采用上一次的值。
        }

        ex_info.Area = deltaS[0];
        ex_info.Err_A = Err[0];
        if(astrom)
        {
            ex_info.astromX = astromX[0];            
        }

        if(Err[0]/fabs(deltaS[0])<=coeff_RelTol * RelTol + RelErrAllowed)
        {
            ex_info.SolveSucceed = true;
        }

        // }
        
    }
}

__device__ void source_base_t::sum_area_3_local ( local_info_t<f_t>& local_info, f_t coeff_RelTol ) const
{
    auto & ex_info = local_info.src_ext;
    auto & ret_info = local_info.src_ret;
    // astrometry
    auto & astrom_Th_info = local_info.src_astrom_Th;

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

    // astrometry
    c_t temp_astromX;
    c_t* astromX = (c_t*) local_info.astromX_sum;


    for(int loop2=0;loop2<2;loop2++)
    {
        deltaS[ threadIdx.x + loop2*n_th] = 0;
        Err[ threadIdx.x + loop2*n_th] = 0;
        if(astrom)
        {
            astromX[ threadIdx.x + loop2*n_th] = 0.;
        }
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
        c_t old_astromX = pt_prev->astromX_interval;
        // printf("init: %d, %.16f,%.16f\n",idx_list[0], old_astromX.re, old_astromX.im);
        if(idx_list[0]< batchidx*n_th)
        {
            for(int j=3;j<5;j++)
            {
                deltaS[ threadIdx.x ] += pt_prev->deltaS[j];
                Err[ threadIdx.x ] += pt_prev->deltaS_Err[j];
                if(astrom)
                {
                    astromX[ threadIdx.x ] += pt_prev->delta_astromX[j];
                }
            }
            for(int j=0;j<3;j++)
            {
                deltaS[ threadIdx.x ] += pt_here->deltaS[j];
                Err[ threadIdx.x ] += pt_here->deltaS_Err[j]; 
                if(astrom)
                {
                    astromX[ threadIdx.x ] += pt_here->delta_astromX[j];
                }
            }
            pt_prev->error_interval = Err[ threadIdx.x ];
            pt_prev->area_interval = deltaS[ threadIdx.x ];
            Err[ threadIdx.x ] -= old_err;
            deltaS[ threadIdx.x ] -= old_area;
            if(astrom)
            {
                pt_prev->astromX_interval = astromX[ threadIdx.x ];
                astromX[ threadIdx.x ] -= old_astromX;
                // printf("old: %d, %d, %.16f,%.16f\n",idx_list[0], idx_list[1], old_astromX.re, old_astromX.im);
                // printf("%d, %d, %.16f, %.16f\n", idx_list[0], idx_list[1], pt_prev->astromX_interval.re, pt_prev->astromX_interval.im);
                // printf("%d, %d, %.16f, %.16f\n", idx_list[0], idx_list[1], astromX[ threadIdx.x ].re, astromX[ threadIdx.x ].im);

                // printf("astromX[ %d ]: %.16f, %.16f,   old_astromX: %.16f, %.16f, old_area_err: %.16f, %.16f, idx_prev: %d, %d,%d\n", threadIdx.x, astromX[ threadIdx.x ].re, astromX[ threadIdx.x ].im, old_astromX.re, old_astromX.im, old_area, old_err, idx_list[0], idx_list[1], idx_list[2]);
            }
        }
        // loop2 = 1 /////////////////////////////
        for(int j=3;j<5;j++)
        {
            deltaS[ threadIdx.x + n_th] += pt_here->deltaS[j];
            Err[ threadIdx.x + n_th] += pt_here->deltaS_Err[j];
            if(astrom)
            {
                astromX[ threadIdx.x + n_th] += pt_here->delta_astromX[j];
            }
        }
        for(int j=0;j<3;j++)
        {
            deltaS[ threadIdx.x + n_th] += pt_next->deltaS[j];
            Err[ threadIdx.x + n_th] += pt_next->deltaS_Err[j];
            if(astrom)
            {
                astromX[ threadIdx.x + n_th] += pt_next->delta_astromX[j];
            }
        }
        pt_here->error_interval = Err[ threadIdx.x + n_th ];
        pt_here->area_interval = deltaS[ threadIdx.x + n_th ];
        if(astrom)
        {
            pt_here->astromX_interval = astromX[ threadIdx.x + n_th ];
            // printf("%d, %d, %.16f, %.16f\n", idx_list[1], idx_list[2], pt_here->astromX_interval.re, pt_here->astromX_interval.im);

            // printf("astromX_interval_here[%d]: %.16f, %.16f\n", idx_list[1], astromX[ threadIdx.x + n_th ].re, astromX[ threadIdx.x + n_th ].im);
        }
    }

    // 后半加到前面  
    deltaS[ threadIdx.x ] += deltaS[ threadIdx.x + n_th ];
    Err[ threadIdx.x ] += Err[ threadIdx.x + n_th ];
    if(astrom)
    {
        astromX[ threadIdx.x ]  += astromX[ threadIdx.x + n_th ];

        // printf("astromX[ %d ]: %.16f, %.16f, deltaS, err: %.16f, %.16f \n", threadIdx.x, astromX[ threadIdx.x ].re, astromX[ threadIdx.x ].im, deltaS[ threadIdx.x ], Err[ threadIdx.x ]);
    }
    __syncthreads();
 

    while(sum_length<(2*n_th))        // 只对前 n_th求和，因为后面的已经加过来了
    {
        if((( threadIdx.x & (sum_length-1))==0) && (( threadIdx.x + sum_length/2) < n_th))    // sum_length = 2^n, idx%sum_length = idx & (sum_length-1)
        {
            deltaS[ threadIdx.x ] += deltaS[ threadIdx.x + sum_length/2];
            Err[ threadIdx.x ] += Err[ threadIdx.x + sum_length/2];
            if(astrom)
            {
                astromX[ threadIdx.x ]  += astromX[ threadIdx.x + sum_length/2];
            }
        }
        sum_length *= 2;
        __syncthreads();
    }


    // // 对 global 修改 cross 并产生面积变化时，更改这里：
    // deltaS[0] += local_info.shared_info->deltaS_cross_global;
    // Err[0] += local_info.shared_info->Err_cross_global;
    // // 疑似废弃，未增加 astrometry 部分

    __syncthreads();


    if(abs(deltaS[0]/ex_info.Area) < coeff_RelTol * RelTol * 1e-2 )
    {
        ex_info.SolveSucceed = true;
    }
    __syncthreads();
    Area = deltaS[0] + ex_info.Area;
    Error = Err[0] + ex_info.Err_A;
    tempS = local_info.src_shape.src_area;
    if(astrom)
    {
        temp_astromX = astromX[0] + ex_info.astromX;
        // if(threadIdx.x==0)
        // {
        //     printf("temp_astromX: %.16f, %.16f, astromX[0]: %.16f, %.16f, ex_info.astromX: %.16f, %.16f  \n", temp_astromX.re, temp_astromX.im, astromX[0].re, astromX[0].im, ex_info.astromX.re, ex_info.astromX.im);
        // }
    }

    // ret_info.mag = fabs(Area) / tempS;
    // ret_info.err = Error / tempS;

    if( (fabs(Area)- 3 * Error) / tempS <= sqrt(1+(2/local_info.src_shape.rho)*(2/local_info.src_shape.rho)))
    {
        ret_info.mag = fabs(Area) / tempS;
        ret_info.err = Error / tempS;
        if(astrom)
        {
            astrom_Th_info = temp_astromX / Area;
        }
    }
    else        // 如果放大率高于理论上界，且高出的值大于三倍误差，认为此时求解出错，取上一次的结果
    {
        ret_info.mag = min(ret_info.mag, sqrt(1+(2/local_info.src_shape.rho)*(2/local_info.src_shape.rho)));
        ret_info.err = -1;
        ex_info.SolveSucceed = true;
        ex_info.Break = true;               // 应该用 Break，但还在测试，所以吧 Break 和 SolveSucceed 都标记了
        // astrom_Th_info = astrom_Th_info;        // 此时不修改 astromX_info，采用上一次的值。
    }

    ex_info.Area = Area;
    ex_info.Err_A = Error;
    if(astrom)
    {
        ex_info.astromX = temp_astromX;        
    }


    if(Error/fabs(Area)<= coeff_RelTol * RelTol + RelErrAllowed)
    {
        ex_info.SolveSucceed = true;
    } 
}

__device__ void source_base_t::adap_set_g ( local_info_t<f_t>& local_info, f_t coeff_RelTol ) const
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

            ErrorWeight[(i*n_th + idx)] = max_f(0,ErrorWeight[(i*n_th + idx)] - coeff_RelTol * RelTol /(batchidx*n_th));        // 如果自己的部分误差小于均值了，那算作过关       
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
    local_info.new_pts[idx].skip = false;
    for(int j=0;j<order-1;j++)
    {
        local_info.new_pts[idx].deltaS[j]=0;
        local_info.new_pts[idx].deltaS_Err[j]=0;
        local_info.new_pts[idx].next_idx[j] = -1;
        local_info.new_pts[idx].next_j[j] = -1;
    }
    __threadfence_block();

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

//////////////////////////////////////////////////////////////////////

__device__ void source_base_t::solve_extended_uniform( local_info_t<f_t>& local_info, const int i_src, f_t coeff_RelTol ) const
{

    if( threadIdx.x == 0 )
    {
        local_info.shared_info->Break = local_info.src_ext.Break;
        local_info.shared_info->Ncross_all = 0;
    }
    __syncthreads();

    // calculation
    for(int bcidx=0;bcidx < (n_point_max / n_th) ;bcidx++)
    // for(int bcidx=0;bcidx < 1 ;bcidx++)
    {
        if((local_info.src_ext.SolveSucceed)||(local_info.shared_info->Break)){break;}

        if(local_info.batchidx==0){margin_set_local( local_info );}
        else{adap_set_g( local_info, coeff_RelTol );}
        __syncthreads();
        margin_solve_local( local_info );
        __syncthreads();
        neighbor3_info_g( local_info );         // g 只体现在读入相邻旧点的 src_pt 上, 可以依据 skip 情况直接修改局域的 prev/next_src_idx
        __syncthreads(); 
        connect_next_local( local_info );
        __syncthreads();

        // if(pt_here->next_idx[jj] == idx_here)

        if( threadIdx.x == 0 )
        {
            local_info.shared_info->Ncross = 0;
            local_info.shared_info->deltaS_cross_global = 0;
            local_info.shared_info->Err_cross_global = 0;
        }
        __syncthreads();
        slope_test_local( local_info, bcidx );
        __syncthreads();
        
        if(local_info.shared_info->Ncross > 0  && (!local_info.shared_info->Break))
        {
            {
                // shared neighbor02 to global
                if( local_info.neighbor3[0] < local_info.batchidx * n_th )
                {
                    for(int j=3;j<5;j++)
                    {
                        pool_margin[i_src * n_point_max + local_info.neighbor3[0]].next_idx[j] = local_info.pt_prev.next_idx[j];
                        pool_margin[i_src * n_point_max + local_info.neighbor3[0]].next_j[j] = local_info.pt_prev.next_j[j];
                    }
                    pool_margin[i_src * n_point_max + local_info.neighbor3[0]].next_src_idx = local_info.pt_prev.next_src_idx;
                }
            pool_margin[i_src * n_point_max + local_info.neighbor3[1]] = local_info.new_pts[ threadIdx.x ];
                if( local_info.neighbor3[2] < local_info.batchidx * n_th )
                {
                    for(int j=0;j<3;j++)
                    {
                        pool_margin[i_src * n_point_max + local_info.neighbor3[2]].next_idx[j] = local_info.pt_next.next_idx[j];
                        pool_margin[i_src * n_point_max + local_info.neighbor3[2]].next_j[j] = local_info.pt_next.next_j[j];
                    }
                    pool_margin[i_src * n_point_max + local_info.neighbor3[2]].prev_src_idx = local_info.pt_next.prev_src_idx;
                }
                __syncthreads();

                slope_detector_g( local_info );
                __syncthreads(); 

                local_info.pt_prev = pool_margin[ blockIdx.x * n_point_max + local_info.neighbor3[0] ];
                local_info.pt_next = pool_margin[ blockIdx.x * n_point_max + local_info.neighbor3[2] ];
                local_info.new_pts[ threadIdx.x ] = pool_margin[ blockIdx.x * n_point_max + local_info.neighbor3[1] ];
                __syncthreads();
            }
        }

        area_err_local( local_info );
        __syncthreads();

        if(local_info.batchidx>=4)
        {
            sum_area_3_local( local_info, coeff_RelTol );
            __syncthreads();             
        }
       
        // shared neighbor02 to global
        if(! local_info.new_pts[ threadIdx.x ].skip)
        {
            if( local_info.neighbor3[0] < local_info.batchidx * n_th )
            {
                // {pool_margin[i_src * n_point_max + local_info.neighbor3[0]] = local_info.pt_prev;}
                for(int j=3;j<5;j++)
                {
                    pool_margin[i_src * n_point_max + local_info.neighbor3[0]].deltaS[j] = local_info.pt_prev.deltaS[j];
                    pool_margin[i_src * n_point_max + local_info.neighbor3[0]].deltaS_Err[j] = local_info.pt_prev.deltaS_Err[j];
                    pool_margin[i_src * n_point_max + local_info.neighbor3[0]].next_idx[j] = local_info.pt_prev.next_idx[j];
                    pool_margin[i_src * n_point_max + local_info.neighbor3[0]].next_j[j] = local_info.pt_prev.next_j[j];
                }
                pool_margin[i_src * n_point_max + local_info.neighbor3[0]].area_interval = local_info.pt_prev.area_interval;
                pool_margin[i_src * n_point_max + local_info.neighbor3[0]].error_interval = local_info.pt_prev.error_interval;
                pool_margin[i_src * n_point_max + local_info.neighbor3[0]].astromX_interval = local_info.pt_prev.astromX_interval;
                // pool_margin[i_src * n_point_max + local_info.neighbor3[0]].next_src_idx = local_info.pt_prev.next_src_idx;       // adap 里面已经把 global 的信息改好了
            }
            pool_margin[i_src * n_point_max + local_info.neighbor3[1]] = local_info.new_pts[ threadIdx.x ];
            if( local_info.neighbor3[2] < local_info.batchidx * n_th )
            {
                // {pool_margin[i_src * n_point_max + local_info.neighbor3[2]] = local_info.pt_next;}
                for(int j=0;j<3;j++)
                {
                    pool_margin[i_src * n_point_max + local_info.neighbor3[2]].deltaS[j] = local_info.pt_next.deltaS[j];
                    pool_margin[i_src * n_point_max + local_info.neighbor3[2]].deltaS_Err[j] = local_info.pt_next.deltaS_Err[j];
                    pool_margin[i_src * n_point_max + local_info.neighbor3[2]].next_idx[j] = local_info.pt_next.next_idx[j];
                    pool_margin[i_src * n_point_max + local_info.neighbor3[2]].next_j[j] = local_info.pt_next.next_j[j];
                }
                // pool_margin[i_src * n_point_max + local_info.neighbor3[2]].prev_src_idx = local_info.pt_next.prev_src_idx;
            }
        }
        else
        {
            pool_margin[i_src * n_point_max + local_info.neighbor3[1]] = local_info.new_pts[ threadIdx.x ];
        }
        __syncthreads();
        // shared neighbor02 to global

        if(bcidx<4)
        {
            sum_area_0_g( local_info, coeff_RelTol );
            __syncthreads();
        }
        
        local_info.batchidx += 1;
    }

    __syncthreads();
    // if(threadIdx.x == 0)                // 最后才输出
    // {
    //     // // local_info.src_ext.SolveSucceed = local_info.shared_info->SolveSucceed;
    //     // local_info.src_ext.Break = local_info.shared_info->Break;
    //     // pool_extended[ i_src ] = local_info.src_ext;
    //     // pool_mag[ i_src ] = local_info.src_ret;
    //     printf("(Ncross: %d, mag: %.6f), ",local_info.shared_info->Ncross, local_info.src_ret.mag);
    // }

}

__device__ f_t source_base_t::M0_phi( local_info_t<f_t>& local_info, const int i_src, const f_t phi ) const
{
    f_t coeff_RelTol = 1./3.;
    f_t r2 = 1-phi*phi;
    local_info.batchidx = 0;
    local_info.src_ext.SolveSucceed = false;
    local_info.src_ext.Break = false;
    local_info.src_shape.src_area = pool_center[ i_src ].src_area * r2;
    local_info.src_shape.rho = pool_center [ i_src ].rho * sqrt(r2);
    solve_extended_uniform(local_info, i_src, coeff_RelTol);

    return local_info.src_ret.mag;
}

}; // namespace twinkle
