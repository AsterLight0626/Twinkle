#pragma once

#include "point_img.h"
#include "point_src.h"

namespace twinkle
{

template < class f_T, int n_th = 32 >
struct local_info_t
{
    f_T                        lens_s;
    // char* dat_sh;
    int                      batchidx;
    src_shape_t< f_T >      src_shape;
    src_ext_t  < f_T >        src_ext;
    src_ret_t  < f_T >        src_ret;
    // 这里要加一个和 src_ret 定位一样的，用来记录局部 astrom 求和的量
    complex_t  < f_T >  src_astrom_Th;
    int                  neighbor3[3];
    // int                 *  global_idx;
    src_pt_t   < f_T >        pt_prev;
    src_pt_t   < f_T >        pt_next;
    // src_pt_t   < f_T >        pt_self;

    f_T                 *  deltaS_sum;
    f_T                 *     Err_sum;
    complex_t  < f_T >  * astromX_sum;

    shared_info_t< f_T >* shared_info;
    cross_info_t        *  cross_info;
    // src_pt_t   < f_T >  *  margin_pts;
    src_pt_t   < f_T >  *     new_pts;

    int                 *  parent_idx;
    f_T                 *    adap_sum;




    __device__ void setup_mem
    ( const char* const dat_sh )
    {
        int offset = 0;

        this -> shared_info = (shared_info_t < f_T > *) &dat_sh[offset];    // size = 1
        offset += sizeof(shared_info_t < f_T > ) * 1 ;

        this -> parent_idx = ( int * ) &dat_sh[offset];         // size = n_th, 这部分虽然在adap 里，但是是独立的，不能共用
        offset += sizeof(int) * n_th;

        // 使用 adap 部分时，后面的部分都没有被使用。它们可以共存于同一空间，不同时间
        this -> adap_sum = ( f_T * ) &dat_sh[offset];           // size = n_point_max，与下面的信息共用
        // offset += sizeof(f_T) * n_point_max;

        this -> deltaS_sum  = ( f_T  * ) &dat_sh[offset];       // size = 4*n_th，从 parent_idx 后面开始写
        offset += sizeof(f_T)   * 4*n_th;

        this -> Err_sum  = ( f_T  * ) &dat_sh[offset];       // size = 4*n_th
        offset += sizeof(f_T)   * 4*n_th;  
        
        this -> astromX_sum  = ( complex_t  < f_T >  * ) &dat_sh[offset];       // size = 4*n_th
        offset += sizeof(complex_t  < f_T >)   * 4*n_th;  

        this -> cross_info = ( cross_info_t  *  ) &dat_sh[offset];      // size =  ( 4 * n_th )
        offset += sizeof(cross_info_t) * ( 4 * n_th );
        
        // 这里还能再加，可以给 astrometry

        // 这里还能再加，可以给 astrometry

        this -> new_pts    = ( src_pt_t   < f_T >  * ) &dat_sh[offset]; // size = ( n_th ), medium part of margin_pts


    }


};


};                              // namespace twinkle