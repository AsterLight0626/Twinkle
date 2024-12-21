#ifdef __HIPCC__

#include "hip.h"
#include <hip/hip_runtime.h>

namespace device
{
inline void hip_api( const hipError_t & e )
{
    if( e != hipSuccess )
        throw std::runtime_error( hipGetErrorString( e ) );
}

static const size_t const_size( 60 * 1024 );

#ifndef   __HIP_CPU_RT__
__constant__
#endif // __HIP_CPU_RT__
char const_pool[ const_size ];

////////////////////////////////////////////////////////////
// Con-/destructor

hip_t::hip_t(  )
{
    f_malloc       = [ & ] ( size_t size ) -> void *
    {
        void * p;
        hip_api( hipMalloc    ( &   p, size ) );
        return p;
    };
    f_malloc_host  = [ & ] ( size_t size ) -> void *
    {
        void * p;
        hip_api( hipHostMalloc( &   p, size ) );
        return p;
    };
    f_free         = [ & ]  ( void * p )
    {
        hip_api( hipFree    (      p ) );
    };
    f_free_host    = [ & ]  ( void * p )
    {
        hip_api( hipHostFree(      p ) );
    };
    f_cp = [ & ] ( void * t, const void * s, size_t n )
    {
        hip_api( hipMemcpy( t, s, n, hipMemcpyDefault ) );
    };
    a_cp = [ & ] ( void * t, const void * s, size_t n,
                   const stream_t & strm )
    {
        hip_api( hipMemcpyAsync
               ( t, s, n, hipMemcpyDefault,  strm ) ) ;
    };
    f_cc =    [ & ] ( void * t, const void * s, size_t n )
    {
        int di( intptr_t( t ) - intptr_t( head_const ) );
#ifdef   __HIP_CPU_RT__
        memcpy( const_pool + di, s, n );
#else
       // hip_api( hipMemcpyToSymbol( head_const, s, n, di ) );
       hip_api( hipMemcpyToSymbol
              ( const_pool, s, n, int( di ) ) );
#endif // __HIP_CPU_RT__
    };
    f_const = [ & ] ( void * t, const void * s, size_t n )
    {
#ifndef __HIP_CPU_RT__
        const
#endif
        void * p = t;
        hip_api( hipMemcpyToSymbol( p, s, n ) );
    };
        
#ifndef __HIP_CPU_RT__          // GPU mode
    f_launch = [ & ]( const void * ker,  const dim3 & n_bl ,
                      const dim3 & n_th, const int  & s_sh ,
                      const void * strm, void  **     args )
    {
        const auto stream = ( hipStream_t ) ( strm     );
        hip_api( hipLaunchKernel
               ( ker, n_bl, n_th, args, s_sh, stream ) );
    };
    f_mset = [ & ]( void         *    p, const int  &  val ,
                    const size_t & size, const void * strm )
    {
        hip_api( hipMemsetAsync
               ( p, val, size, hipStream_t( strm ) ) );
    };
    max_streams = 2048;
#else                           // CPU mode
    head_const = const_pool;
    f_exec = [ & ]( const dim3 & n_bl, const dim3 & n_th,
                    const int  & s_sh, const void * strm,
                    std::function< void(  ) > ker )
    {
        const auto stream = ( hipStream_t ) ( strm );
        hipLaunchKernelGGL( ker, n_bl, n_th, s_sh, stream );
    };
    f_mset = [ & ]( void         *    p, const int  &  val ,
                    const size_t & size, const void * strm )
    {
        launch_host( ( hipStream_t ) ( strm ), [ = ] (  )
        {  return std::memset( p, val, size ); }   ) ;
    };
    max_streams = 244;    
#endif // __HIP_CPU_RT__
    size_const  = const_size;
    idx_streams =          0;
    rng_hd      =    nullptr;    
    return;
}

////////////////////////////////////////////////////////////
// Initialization

void hip_t::init( const int & rank )
{
    base_t::init ( rank );
    if( idx_device <  0 )
    {
        hip_api( hipGetDeviceCount( & num_device ) );
        idx_device = rank  % num_device  ;
    }
    this->prepare(      );

#ifndef   __HIP_CPU_RT__
    hip_api( hipDeviceSetSharedMemConfig
             ( hipSharedMemBankSizeEightByte ) );
    hip_api( hipDeviceSetCacheConfig
             ( hipFuncCachePreferShared    ) );
#endif // __HIP_CPU_RT__
    return;
}

void hip_t::prepare (  )
{
    hip_api( hipSetDevice( idx_device ) );
#ifndef   __HIP_CPU_RT__
    hip_api( hipGetSymbolAddress
           ( & head_const, const_pool ) );
#endif // __HIP_CPU_RT__
    size_const           = const_size;    
    return;
}

void hip_t::finalize(  )
{
    delete_all_streams(  );
    return;
}

////////////////////////////////////////////////////////////
// Streams

void hip_t::  pin( void * p, const size_t & s )
{
#ifndef   __HIP_CPU_RT__
    return hip_api
    ( hipHostRegister( p, s, hipHostRegisterDefault ) );
#endif // __HIP_CPU_RT__
}

void hip_t::unpin( void * p )
{
#ifndef   __HIP_CPU_RT__
    return hip_api( hipHostUnregister( p ) );
#endif // __HIP_CPU_RT__
}

stream_t hip_t::yield_stream( const int & pri )
{
    if( streams.size(  ) < max_streams )
    {
        streams.emplace_back(  );
        hip_api( hipStreamCreateWithFlags
        ( & streams.back(  ), hipStreamNonBlocking ) );
    }
    auto   res = streams[ idx_streams ];
    idx_streams += 1 ;
    idx_streams  = idx_streams % max_streams;
    return  res ;
}

void hip_t::sync_stream_base( const stream_t & stream )
{
    return hip_api( hipStreamSynchronize( stream ) );
}

void hip_t::sync_all_streams(  )
{
    hip_api( hipDeviceSynchronize(  ) );
    return base_t::sync_all_streams(  );
}

void hip_t::delete_all_streams(  )
{
    // sync_all_streams(  );
    for( auto & p :  streams )
    {
        hip_api( hipStreamSynchronize( p ) );
        hip_api( hipStreamDestroy    ( p ) );
    }
    return streams.clear(  );
}

////////////////////////////////////////////////////////////
// Events

event_t hip_t::event_create(  )
{
    event_t e;
    hip_api( hipEventCreate( & e ) );
    return  e;
}

void hip_t::event_destroy( const event_t & event )
{
    return hip_api( hipEventDestroy( event ) );
}

void hip_t::event_record ( const event_t  &  event ,
                           const stream_t & stream )
{
#ifdef    __HIP_CPU_RT__
    event_m[ event ] = stream;
#else  // __HIP_CPU_RT__
    return hip_api( hipEventRecord( event, stream ) );
#endif // __HIP_CPU_RT__
}

void hip_t::event_wait   ( const stream_t & stream ,
                           const event_t  &  event )
{
#ifdef    __HIP_CPU_RT__
    event_sync (  event );
    sync_stream( stream );
#else  // __HIP_CPU_RT__
   return hip_api( hipStreamWaitEvent( stream, event, 0 ) );
#endif // __HIP_CPU_RT__
}

void hip_t::event_sync   ( const event_t & event )
{
#ifdef    __HIP_CPU_RT__
    if( event_m.find( event ) != event_m.end(  ) )
    {
        sync_stream_base( event_m[ event ] );
        event_m.erase            ( event )  ;
    }
#else  // __HIP_CPU_RT__
    return hip_api( hipEventSynchronize  ( event ) );
#endif // __HIP_CPU_RT__
}

////////////////////////////////////////////////////////////
// Callbacks

void hip_t::f_launch_host
( const stream_t & s, const f_cb_t f, void * p )
{
#ifdef    __HIP_CPU_RT__
    f_exec( 1, 1, 0, s, [ = ] (  )
    {
        ( * reinterpret_cast< func_t * >( p ) ) (  );
    }   ) ;
#else  // __HIP_CPU_RT__
    return hip_api( hipLaunchHostFunc( s, f, p ) );
#endif // __HIP_CPU_RT__
}

////////////////////////////////////////////////////////////
// Device-side RNG

#ifndef     __HIP_CPU_RT__
__constant__  hip_t::rand_t * rng_c;
__constant__ size_t         n_rng_c;

__global__ void rng_init_ker( int seed, hip_t::rand_t * p )
{
    // const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    // if( idx < n_rng_c )
    //     hiprand_init( seed, threadIdx.x, idx, & p[ idx ] );
    return;
}
#endif // __HIP_CPU_RT__

void hip_t::prep_rng( const size_t & n_rng )
{
#ifdef    __HIP_CPU_RT__
    srand( time( NULL ) );
#else  // __HIP_CPU_RT__
    throw std::runtime_error( "RNG not ready for AMD GPU!" );
    // hipMalloc( ( void ** ) & rng_hd,
    //             n_rng * sizeof( rand_t ) );
    // const int seed( 77 );
    // const int n_th( 64 );
    // const int n_bl( ( n_rng + n_th - 1 ) / n_th );
    // hip_api( hipMemcpyToSymbol
    //          ( n_rng_c, & n_rng, sizeof( int   ) ) );
    // rng_init_ker<<< n_bl, n_th >>>( seed, rng_hd );
    // hip_api( hipMemcpyToSymbol
    //          ( rng_c, & rng_hd, sizeof( rng_hd ) ) );
    // return hip_api( hipDeviceSynchronize (  ) );
#endif // __HIP_CPU_RT__
}

__device__ float rand_dev(  )
{
#ifdef    __HIP_CPU_RT__
    return std::rand(  ) / ( RAND_MAX + 1.f );    
#else  // __HIP_CPU_RT__
    return 0.f;
    // const auto i_th = threadIdx.x + blockDim.x * blockIdx.x;
    // float i_rng = int( float( i_th )  / n_rng_c );
    // i_rng       = i_th - i_rng * n_rng_c;
    // return curand_uniform( rng_c + int( i_rng ) );
#endif // __HIP_CPU_RT__
}

};                              // namespace device

#endif // __HIPCC__
