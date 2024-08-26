#pragma once
#include<iostream>
#include<fstream>
#include "../device/general.h"

#include"MacroVaribles.h"
#include"StructComplex.h"
#include"StructSrc.h"
#include"GlobalFuncs_v2.h"
#include"Caustics.h"

// #define cout std::cout
// #define endl std::endl
// using namespace std;
// using namespace chrono;



// #ifndef __HIP_CPU_RT__
// #ifdef        __CUDACC__
// #define __arch__( X ) cuda##X
// #elif defined( __HIPCC__ )
// #define __arch__( X )  hip##X
// #endif // __CUDACC__ or __HIPCC__
// #define __launch__( KER, N_BL, N_TH, S_SH, STRM, ... ) \
//     KER <<< N_BL, N_TH, S_SH, STRM >>> ( __VA_ARGS__ )
// #else  // __HIP_CPU_RT__
// #define __arch__( X )  hip##X
// #define __launch__( KER, N_BL, N_TH, S_SH, STRM, ... ) \
//     hipLaunchKernelGGL                                 \
//     ( KER, N_BL, N_TH, S_SH, STRM, __VA_ARGS__ )
// #endif // __HIP_CPU_RT__


// using stream_t = __arch__(Stream_t);

enum on_dev_t { dev = 0, host = 1, undef = -1 };

template < class f_T, int dim = 1, class i_T = int >
struct array_t
{
    bool on_device;
    f_T     * data;
    i_T       size;
    i_T  n [ dim ];

    ////////// Host-side operations //////////
    template< class idx_T > 
    __host__ void init
    ( const device_t & dev, const idx_T &  shape ,
      const on_dev_t & on_device )
    {
        if constexpr( dim == 1 )
            n[ 0 ] = shape;
        else
            for( int a = 0; a < dim; ++ a )
                n[ a ] = shape[ a ];
        size = 1;        
        for( int a = 0; a < dim; ++ a )        
            size  *= n[ a ];
        this->on_device  = ( on_device == on_dev_t:: dev );
        data = ( this->on_device 
                ? dev.malloc_device< f_T > ( size ) 
                : dev.malloc_host  < f_T > ( size ) );
        return;
    };
    __host__ void free( const device_t & dev )
    {
        return on_device ? dev.free_device( data ) 
                         : dev.free_host  ( data );
    };

    template< class tgt_T >
    __host__ void cp_to
    ( const device_t & dev, tgt_T & tgt, 
      const device::stream_t & stream )
    {
        for( int a = 0; a < dim; ++ a )
            if( tgt.n[ a ] != n[ a ] )
                throw std::runtime_error
                    ( "cp size inconsistent!" );     
        return dev.cp_a( tgt.data, data, size, stream );        // device.h
    };
    
    ////////// Dual-side operations //////////
    template < class idx_T > __forceinline__        // 编译优化，函数强制要求展开
    __host__ __device__ i_T idx( const idx_T & i ) const
    {
        i_T  res( i[ dim - 1 ] );
        for( int a = dim - 2; a >= 0; -- a )
		{
#ifdef __HIP_CPU_RT__
			if( i[ a ] < 0 || i[ a ] >= n[ a ] )
				throw std::runtime_error( "invalid index" );
#endif
            res = res * n[ a ] + i[ a ];
		}
        return res;
    };

    template < class idx_T >
    __forceinline__ __host__ __device__
    f_T & operator[  ] ( const idx_T & i ) const
    {
        // __dyn_shared__( f_T, p_sh );
        // extern __shared__ f_T p_sh[  ];
        f_T * data = const_cast< f_T * > ( this->data );
        if constexpr( dim == 1 )
		{
#ifdef __HIP_CPU_RT__
			if( i >= size )
				throw std::runtime_error( "invalid index" );
#endif
            return data[ i ];
		}
        else
            return data[ this->idx( i ) ];
    };
};


// nsys profile --stats=true ./bin/Twinkle     

// including parameters from MacroVaribles.h
template<isFloating f_T>
class Twinkle
{
    public:
    int Nsrcs = NSRCS;
    // int DEV_NUM;
    device_t dev;
    device::stream_t stream;
    // src_params_t<f_T>* h_params;
    // src_params_t<f_T>* d_params;
    array_t<src_ext_t<f_T>> h_srcs,d_srcs;
    // array_t<complex_t<f_T>> h_crit,d_crit;
    // array_t<complex_t<f_T>> h_caus,d_caus;
    // array_t<int> h_connect,d_connect;
    array_t<CC_t<f_T>> h_CC,d_CC;
    // int* h_Nfail;        //  initially 0
    // int* d_Nfail;
    // array_t<params_back_t<f_T>> h_params_back,d_params_back;
    array_t<int> h_srcidx,d_srcidx;
    int Ncal;           //  initally NSRCS
    // int* h_Ncal;         //  initally NSRCS
    // int* d_Ncal;

    array_t<src_params_t<f_T>> h_params,d_params;
    array_t<int> h_Nfail,d_Nfail,h_Ncal,d_Ncal;

    array_t<f_T> h_Mag,d_Mag;
    array_t<f_T> h_Err,d_Err;

    array_t<f_T> h_time,d_time;

    



    // public:x

    bool muted = false;      // if true, no warning information
    bool critcaus = false;

    Twinkle(){};
    void malloc(int device_num = 0);
    void free();
    void init();

    void set_time();
    void set_time(f_T* time);

    void set_params();
    void set_params(src_params_t<f_T>* params);

    void set_path();                        // set path is always necessary, even set time and params in previous
    void set_path(f_T* time);
    void set_path(src_params_t<f_T>* params);
    void set_path(f_T* time, src_params_t<f_T>* params);

    void set_array(int NaxisX=30, int NaxisYf_T=30, f_T xmin=-3e-5, f_T xmax=3e-5, f_T ymin=1e-4, f_T ymax=1e-4);

    void set_rho_loc(f_T* rhos, complex_t<f_T>* zetas);

    void solve();
    void cp_back();
    void writeto(std::string path, std::string suffix = "");

    void cp_back_ME();
    void writeto_ME(std::string path, std::string suffix = "");

    void cp_back_caus();
    void solve_caus();
    void writeto_caus(std::string path, std::string suffix = "");

};


template<isFloating f_T>
_Global void SrcsInit(src_ext_t<f_T>* srcs);

_Global void ResTest();