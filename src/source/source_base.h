#pragma once

#include "../utils/complex.h"
#include "../device/general.h"
#include "./types/point_img.h"
#include "./types/point_src.h"
#include "./types/extended_src.h"
#include "./types/pool.h"
#include "./types/sh_manager.h"
#include "./types/LimbDarkening.h"
#include <map>
#include <stack>
#include <queue>
// #include <vector>


namespace twinkle
{

// using  float_t = double;
// using float2_t = double;
using  float_t = double;
using float2_t = double;

////////////////////////////////////////////////////////////
//

template< class src_T >
__global__ void point_approximation( const src_T src )
{
    src.init_break_succeed(  );
    return src.solve_point_approx(  );
}


template< class src_T >
__global__ void solve_extended_sh( const src_T src )
{
    __dyn_shared__( char, dat_sh );

    const int i_src = blockIdx.x;
    if( i_src >= src.n_src )
        return;
    if( threadIdx.x >= src.n_th )
        return;

    local_info_t< float2_t > local_info;

    // input
    local_info.setup_mem( dat_sh );
    local_info.src_shape = src.pool_center [ i_src ];
    local_info.src_ext = src.pool_extended [ i_src ];
    local_info.src_ret = src.pool_mag      [ i_src ];
    local_info.lens_s = src.pool_lens_s [ i_src ];

    local_info.batchidx = 0;
    src.solve_extended_uniform(local_info, i_src);


    // 检查 pool_margin, 找到所有自连接次数
    for(int i_batch=0;i_batch<local_info.batchidx;i_batch++)
    {
        for(int j=0;j<5;j++)
        {
            if(src.pool_margin[i_src * src.n_point_max + i_batch*src.n_th + threadIdx.x].next_idx[j] == threadIdx.x + i_batch*src.n_th)
            {
                atomicAdd(&(local_info.shared_info->Ncross_all), 1);
            }
        } 
    }
    __syncthreads();

    if(threadIdx.x == 0)                // 最后才输出
    {
        // local_info.src_ext.SolveSucceed = local_info.shared_info->SolveSucceed;
        local_info.src_ext.Break = local_info.shared_info->Break;
        src.pool_extended[ i_src ] = local_info.src_ext;
        src.pool_mag[ i_src ] = local_info.src_ret;
        // printf("(%d, Ncross: %d, mag: %.6f), ",i_src, local_info.shared_info->Ncross, local_info.src_ret.mag);
        src.pool_Ncross[ i_src ] = local_info.shared_info->Ncross_all;

        // if(local_info.src_ext.Break)
        // {
        //     src.pool_mag[ i_src ].mag = -1.;
        // }
    }

}



template< class src_T >
__global__ void solve_LD_sh( const src_T src )
{
    __dyn_shared__( char, dat_sh );

    const int i_src = blockIdx.x;
    if( i_src >= src.n_src )
        return;
    if( threadIdx.x >= src.n_th )
        return;
    // float2_t phi = phi_a[ i_src ];
    float2_t phi = src.pool_phi[ i_src ];
    if ((phi >= 1) || (phi<0))          // 用这个机制让已经完成计算的点跳过计算
        return;

    float2_t M0;
    local_info_t< float2_t > local_info;

    // input
    local_info.setup_mem( dat_sh );
    local_info.src_shape = src.pool_center [ i_src ];
    local_info.src_ext = src.pool_extended [ i_src ];
    local_info.src_ret = src.pool_mag      [ i_src ];
    local_info.lens_s = src.pool_lens_s [ i_src ];

    // auto M0r2_phi = [&src, &local_info, i_src] __device__ (auto phi) {
    //     return src.M0r2(local_info, i_src, phi);
    // };
    // M0r2 = M0r2_phi(phi);

    // // r=1
    // phi = 0;
    // r2 = 1-phi*phi;
    // local_info.batchidx = 0;
    // local_info.src_ext.SolveSucceed = false;
    // local_info.src_ext.Break = false;
    // local_info.src_shape.src_area = src.pool_center[ i_src ].src_area * r2;
    // local_info.src_shape.rho = src.pool_center [ i_src ].rho * sqrt(r2);
    // src.solve_extended_uniform(local_info, i_src);

    // phi = 0.5;
    M0 = src.M0_phi(local_info, i_src, phi);


    if(threadIdx.x == 0)                // 最后才输出
    {
        // // local_info.src_ext.SolveSucceed = local_info.shared_info->SolveSucceed;
        // local_info.src_ext.Break = local_info.shared_info->Break;
        // src.pool_extended[ i_src ] = local_info.src_ext;
        // src.pool_mag[ i_src ] = local_info.src_ret;
        src.pool_mag[ i_src ].mag = M0;
        src.pool_Ncross[ i_src ] = local_info.shared_info->Ncross_all;
    }

}


////////////////////////////////////////////////////////////
//




////////////////////////////////////////////////////////////
//

static const float C_Q = 6;
static const float C_G = 2;
static const float C_P = 2;
static const float2_t PI = 3.141592653589793;
static const double RelErrAllowed = 1e-10;

class source_base_t
{
    ////////// Types //////////
public:                      // Types
    using f_t =            float2_t;
    using c_t = complex_t< float2_t >;
    using img_t = image_pt_t< f_t >;

    ////////// Host-side initialization //////////
public:
    __host__  source_base_t(  );
    __host__ ~source_base_t(  );
    
    ////////// Source data maintainance //////////
public:                      // Data
    int                           n_cross_max;
    int                           n_point_max;
    pool_t< src_ret_t   < f_t > >    pool_mag;
    pool_t<      cross_info_t   >  pool_cross;
    pool_t< src_shape_t < f_t > > pool_center;
    pool_t< src_pt_t    < f_t > > pool_margin;
    pool_t< src_ext_t   < f_t > > pool_extended;
    pool_t<               f_t   > pool_lens_s;
    pool_t<               f_t   >    pool_phi;
    pool_t<               int   > pool_Ncross;

public:                         // Functions
    __host__ virtual void init( device_t & f_dev );
    // __host__ virtual void set_params_2D( device_t & f_dev, f_t ss, f_t qq, f_t rho, f_t xmax, f_t xmin, f_t ymax, f_t ymin, int Nx, int Ny );
    // __host__ virtual void set_params_1D( device_t & f_dev, f_t ss, f_t qq, f_t rho, f_t xmax, f_t xmin, f_t ymax, f_t ymin, int Nsrc );
    // __host__ virtual void set_same( device_t & f_dev, f_t ss, f_t qq, f_t rho, f_t zeta_x, f_t zeta_y );
    __host__ virtual void free( device_t & f_dev );    

    ////////// Device call parameters //////////

    int             n_src;
    static constexpr int n_th = 32;
    int       streamidx;
protected:                      // Data
    device::stream_t    stream;
    int     size_batch;
    // int          n_src;
    // int           n_th;

    ////////// Lensing parameters //////////
// protected:                      // Data
public:
    static constexpr int order = 6;
    // f_t         s;
    f_t         lens_q;
    f_t    RelTol;
    f_t        m1;
    f_t        m2;
    f_t        a1;          // for linear limb darkening
    // int  batchidx;
    
    ////////// Device-side interfaces //////////
// protected:                      // Functions
public:
    __device__ f_t  yield_lens_coef
    ( c_t   * coef, const c_t & zeta, const f_t lens_s ) const;
    __device__ bool solve_lens_eq
    ( img_t *  img,       c_t * coef,  const bool use_init_roots = false) const;

    __host__ __device__ c_t f_zbar
    ( const c_t & z, const f_t len_s  ) const; 
    __device__ void bubble_arg_sort
    ( int* index, f_t* arr, const int len ) const;        
    __device__ void is_physical_all
    ( img_t * imgs, const c_t & zeta, int & Nphys, const f_t lens_s  ) const; 
    __device__ f_t jacobian
    ( const c_t & z, const f_t lens_s ) const;  
    __device__ f_t mu_qc
    ( const c_t & z, const f_t lens_s ) const;
    __device__ f_t ghost_tst
    ( const c_t & z_G, const c_t & z_hat, const f_t lens_s ) const;   

    __device__ bool solve_imgs
    ( img_t * imgs, const c_t & zeta, const bool use_init_roots, const f_t lens_s,  c_t * roots_point = nullptr) const;
    __device__ f_t point_mag
    ( img_t * imgs, const c_t & zeta, int & Nphys ) const;
    __device__ bool ghost_test_all
    ( const img_t * const imgs, const int & Nphys, const c_t & zeta, const f_t & Rho, const f_t lens_s ) const;
    __device__ bool safe_distance_test
    ( const c_t & zeta, const f_t & Rho, const f_t lens_s ) const;


    __device__ void margin_Q
    ( c_t & marginloc, const src_shape_t<f_t>& shape, const float& Q ) const;

    __device__ void set_skip_info
    ( src_pt_t < f_t > & src ) const;
    __device__ void dz_wedge
    ( c_t & dz, f_t & wedge, f_t & jacobian, const img_t & img, const c_t & zeta, const c_t & loc_centre, const f_t rho, const f_t lens_s ) const;

    __device__ bool polish_pp
    ( img_t * imgs, const f_t * const jacobians) const;
    __device__ bool polish_pp2
    ( img_t * imgs, const f_t * const jacobians, const c_t & zeta) const;    
    __device__ bool check_parity
    ( img_t * imgs, const f_t * const jacobians, const c_t & zeta, const int Nphys) const;

    __device__ bool margin_solve_cal                // return skip
    ( img_t * imgs, c_t * dz, f_t * wedge, int & phys_tst, int * temp_j_from_out_j, const c_t & zeta, const c_t & loc_centre, const f_t Rho, const f_t len_s) const;

    

    __device__ int prev_src_idx_g       // "g" means communicate to global memory
    ( const int idx_here, const int srcidx ) const;    
    __device__ int next_src_idx_g
    ( const int idx_here, const int srcidx ) const; 
    __device__ int prev_src_idx_local
    ( const int idx_here, local_info_t<f_t>& local_info, bool * changed = nullptr ) const;
    __device__ int next_src_idx_local
    ( const int idx_here, local_info_t<f_t>& local_info, bool * changed = nullptr ) const;   
    __device__ int prev_src_idx_local_g
    ( const int idx_here, local_info_t<f_t>& local_info, const int srcidx ) const;
    __device__ int next_src_idx_local_g
    ( const int idx_here, local_info_t<f_t>& local_info, const int srcidx ) const;    


    __device__ void connect_order
    ( const c_t* posi_here, const c_t* posi_other, int len_here, int len_other, int* order) const;
    __device__ int positive_connect_j34     // return j_positive_candidate
    ( int* next_idx_out, int* next_j_out, const int Nphys_here, const int Nphys_next,
     const c_t * posi_here_positive, const c_t * posi_next, const int here_idx, const int next_idx) const;
    __device__ int negative_connect_j012    // return j_negative_candidate
    ( int* next_idx_out, int* next_j_out, const int Nphys_here, const int Nphys_prev,
     const c_t * posi_here_negative, const c_t * posi_prev, const int here_idx, const int prev_idx) const;
    __device__ bool cross_j34
    ( int* next_j_out, const int j_positive_candidate,
     const c_t * posi_here_negative, const c_t * posi_next_negative, const int here_idx) const;
    __device__ bool cross_j012               // 可能修改012
    ( int* next_j_out, const int j_negative_candidate,
    const c_t * posi_here_positive, const c_t * posi_prev_positive, const int here_idx) const;

    __device__ bool connect_prev_j012
    ( int* next_idx_out, int* next_j_out,
    const img_t * imgs_here, const img_t * imgs_prev, const int here_idx, const int prev_idx, 
    const int Nphys_here, const int Nphys_prev) const;
    __device__ bool connect_next_j34
    ( int* next_idx_out, int* next_j_out,
    const img_t * imgs_here, const img_t * imgs_next, const int here_idx, const int next_idx,
    const int Nphys_here, const int Nphys_next) const;


    __device__ f_t delta_s_1
    ( const c_t & z_here, const c_t & z_next) const;
    __device__ f_t wedge_product
    ( const c_t &z1, const c_t &z2 ) const;
    __device__ void deltaS_error_image
    ( f_t & deltaS, f_t & Error1234, const int j,
    const src_pt_t < f_t > * pt_here, const src_pt_t < f_t > * pt_other,
    const f_t & theta1, const f_t & theta2, const f_t & theta3) const;
    __device__ void deltaS_error_image_beta
    ( f_t & deltaS, f_t & Error1234, const int j,
    const src_pt_t < f_t > * pt_here, const src_pt_t < f_t > * pt_other,
    const f_t & theta1, const f_t & theta2, const f_t & theta3) const;
    __device__ void deltaS_error_parity
    ( f_t * deltaS_new_out, f_t * Err_new_out, const bool parity,
    const int here_idx, const int other_idx,
    const src_pt_t < f_t > * pt_here, const src_pt_t < f_t > * pt_other)    const;

    __device__ void deltaS_error_cross
    ( f_t & deltaS, f_t & Error1234, const int j0, const int j1, const bool ghost_direction,
    const src_pt_t < f_t > * pt_here ) const;
    __device__ void deltaS_error_cross_beta
    ( f_t & deltaS, f_t & Error1234, const int j0, const int j1, const bool ghost_direction,
    const src_pt_t < f_t > * pt_here ) const;

    __device__ int bisearch_left
    (const f_t* values_in_order, const f_t to_insert, const int len) const;


    __device__ bool slope_test
    ( src_pt_t < f_t > * pt_here, src_pt_t < f_t > * pt_other,
    const int jj, const int idx_here, bool& Break, int& Ncross_all, const int bcidx ) const;


public:                         // Functions
    __device__ void init_break_succeed(  ) const;
    __device__ void solve_point_approx( f_t coeff_RelTol=1. ) const;
    
    __device__ void margin_set_local  ( local_info_t<f_t>& local_info ) const;
    __device__ void margin_solve_local( local_info_t<f_t>& local_info ) const;

    __device__ void neighbor3_info_g  ( local_info_t<f_t>& local_info ) const;

    __device__ void connect_next_local( local_info_t<f_t>& local_info ) const;


    __device__ void slope_test_local  (local_info_t<f_t>& local_info, const int bcidx ) const;
    __device__ void slope_detector_g  ( local_info_t<f_t>& local_info ) const;

    __device__ void area_err_local    ( local_info_t<f_t>& local_info ) const;

    __device__ void sum_area_0_g      ( local_info_t<f_t>& local_info, f_t coeff_RelTol=1. ) const;
    __device__ void sum_area_3_local  ( local_info_t<f_t>& local_info, f_t coeff_RelTol=1. ) const;
    

    __device__ void adap_set_g        ( local_info_t<f_t>& local_info, f_t coeff_RelTol=1. ) const;

    __device__ void solve_extended_uniform( local_info_t<f_t>& local_info, const int i_src, f_t coeff_RelTol=1. ) const;
    __device__  f_t M0_phi( local_info_t<f_t>& local_info, const int i_src, const f_t phi ) const;

    ////////// Host-side interfaces //////////
public:
    __host__ virtual void run( device_t & f_dev );
    __host__ virtual void run_pt( device_t & f_dev );
    // __host__ virtual void runLD( device_t & f_dev, double LD_a=1, int depth=4 );
    __host__ virtual void runLD_A( device_t & f_dev, double LD_a=1, int max_depth=10 );
    __host__ virtual void runLD_MaxErr( device_t & f_dev, double LD_a=1, int* Nuniform_out=nullptr, int max_depth=30 );
    __host__ virtual void runLD_beta( device_t & f_dev, double LD_a=1, int* Nuniform_out=nullptr, int max_depth=30 );

    __host__ virtual void monotest_single(double val_L, double val_R, double Mag_L, double Mag_R, double phi_self, double& val_self, double& Mag_self);
    __host__ virtual void monotest(twinkle::ErrorUnit_t<double>& eu);
    __host__ virtual void caustic_cal(double lens_s, double lens_q, c_t* caustic_points);
    __host__ virtual void hidden_phi(c_t* caustic_pts, const c_t& loc_center, double rho, double* phi_hidden);


};                              // class source_base_t

};                              // namespace twinkle
