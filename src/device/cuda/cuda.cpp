#ifdef __CUDACC__

#include <stdio.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>

#include "cuda.h"

namespace device
{
////////////////////////////////////////////////////////////
// Con-/destructor

inline void cuda_api( const cudaError_t & e )
{
    if( e != cudaSuccess )
        throw std::runtime_error( cudaGetErrorString( e ) );
}


static const size_t const_size( 60  * 1024 );
__constant__ char   const_pool[ const_size ];

cuda_t::cuda_t(  )
{
    f_malloc       = [ & ]   ( size_t size ) -> void *
    {
        void * p;
        cuda_api( cudaMalloc( & p, size ) );
        return p;
    };
    f_malloc_host  = [ & ]   ( size_t size ) -> void *
    {
        void * p;
        cuda_api( cudaMallocHost( & p, size ) );
        return p;
    };
    f_free           = [ & ] ( void * p )
    {
        cuda_api( cudaFree( p ) );
    };
    f_free_host      = [ & ] ( void * p )
    {
        cuda_api( cudaFreeHost( p ) );
    };
    f_cp = [ & ] ( void * t, const void * s, size_t n )
    {
        cuda_api
            ( cudaMemcpy( t, s, n, cudaMemcpyDefault ) );
    };
    a_cp = [ & ] ( void * t, const void * s, size_t n ,
                   const stream_t & strm )
    {
        cuda_api( cudaMemcpyAsync
                  ( t, s, n, cudaMemcpyDefault, strm ) );
    };
    f_cc    = [ & ] ( void * t, const void * s, size_t n )
    {
        auto di( intptr_t( t ) - intptr_t( head_const ) );
        cuda_api( cudaMemcpyToSymbol
                  ( const_pool, s, n, int( di ) ) );
    };
    f_const = [ & ] ( void * t, const void * s, size_t n )
    {
        cuda_api( cudaMemcpyToSymbol
                  ( ( const void * ) t, s, n, 0,
                    cudaMemcpyHostToDevice ) );
    };
    f_launch = [ & ]( const void * ker,  const dim3 & n_bl ,
                      const dim3 & n_th, const int  & s_sh ,
                      const void * strm, void  **     args )
    {
        const auto stream = ( cudaStream_t ) ( strm );
        cuda_api( cudaLaunchKernel
                  ( ker, n_bl, n_th, args, s_sh, stream ) );
    };
    f_mset = [ & ]( void         *    p, const int  &  val ,
                    const size_t & size, const void * strm )
    {
        cuda_api( cudaMemsetAsync
                  ( p, val, size, stream_t( strm ) ) );
    };
    idx_streams =          0;
    max_streams =         64;
    rng_hd      =    nullptr;
    return;
}

////////////////////////////////////////////////////////////
// Initialization

void cuda_t::init( const int & rank )
{
    base_t::init ( rank );
    if( idx_device <  0 )
    {
        cuda_api( cudaGetDeviceCount( & num_device ) );
        idx_device =  rank  % num_device  ;
    }
    this->prepare(      );

    cuda_api( cudaDeviceSetCacheConfig
              ( cudaFuncCachePreferShared ) );
    return;
}

void cuda_t::prepare (  )
{
    cuda_api( cudaSetDevice ( idx_device ) );
    cuda_api( cudaGetSymbolAddress
              ( & head_const, const_pool ) );
    size_const  = const_size;
    return;
}

void cuda_t::finalize(  )
{
    if( rng_hd != nullptr )
        cuda_api( cudaFree( rng_hd ) );
    return delete_all_streams(  );
}

////////////////////////////////////////////////////////////
// Streams

void cuda_t::  pin( void * p, const size_t & s )
{
    return cuda_api( cudaHostRegister
                     ( p, s, cudaHostRegisterDefault ) );
}

void cuda_t::unpin( void * p )
{
    return cuda_api( cudaHostUnregister( p ) );
}

stream_t cuda_t::yield_stream( const int & pri )
{
    if( streams.size(  ) < max_streams )
    {
        streams.emplace_back(  );
        cuda_api( cudaStreamCreateWithPriority
        ( & streams.back( ), cudaStreamNonBlocking, pri ) );
    }
    auto    res  = streams[ idx_streams ];
    idx_streams += 1 ;
    idx_streams  = idx_streams % max_streams;
    return  res ;
}

void cuda_t::sync_stream_base( const stream_t & stream )
{
    return cuda_api( cudaStreamSynchronize( stream ) );
}

void cuda_t::sync_all_streams(  )
{
    cuda_api( cudaDeviceSynchronize(  ) );
    return base_t::sync_all_streams(  );
}

void cuda_t::delete_all_streams(  )
{
    sync_all_streams(        );
    for( auto & p :  streams )
        cuda_api( cudaStreamDestroy( p ) );
    return streams. clear(   );
}

////////////////////////////////////////////////////////////
// Events

event_t cuda_t::event_create(  )
{
    event_t e ;
    cuda_api( cudaEventCreateWithFlags
            ( & e, cudaEventDisableTiming ) );
    return  e ;
}

void cuda_t::event_destroy( const event_t & event )
{
    return cuda_api( cudaEventDestroy( event ) );
}

void cuda_t::event_record ( const  event_t &  event ,
                            const stream_t & stream )
{
    return cuda_api( cudaEventRecord( event, stream ) );
}

void cuda_t::event_wait   ( const stream_t & stream ,
                            const  event_t &  event )
{
    return cuda_api( cudaStreamWaitEvent( stream, event ) );
}

void cuda_t::event_sync   ( const event_t & event )
{
    return cuda_api( cudaEventSynchronize( event ) );
}

////////////////////////////////////////////////////////////
// Callbacks

void cuda_t::f_launch_host
( const stream_t & s, const f_cb_t f, void * p )
{
    return cuda_api( cudaLaunchHostFunc( s, f, p ) );
}

////////////////////////////////////////////////////////////
// Device-side RNG

// __constant__ cuda_t::rand_t * rng_c;
// __constant__ size_t         n_rng_c;

__global__ void rng_init_ker( int seed, cuda_t::rand_t * p )
{
    // const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    // if( idx < int( n_rng_c ) )
    //     curand_init( seed, threadIdx.x, idx, & p[ idx ] );
    // return;
}

void cuda_t::prep_rng( const size_t & n_rng )
{
    // cuda_api( cudaMalloc( ( void ** ) & rng_hd,
    //                       n_rng * sizeof( rand_t ) ) );
    // const int seed( 77 );
    // const int n_th( 64 );
    // const int n_bl( ( n_rng + n_th - 1 ) / n_th );
    // cuda_api( cudaMemcpyToSymbol
    //           ( n_rng_c, & n_rng, sizeof( n_rng ) ) );
    // rng_init_ker<<< n_bl, n_th >>>  ( seed, rng_hd  );
    // cuda_api( cudaMemcpyToSymbol
    //           ( rng_c, & rng_hd, sizeof( rng_hd ) ) );
    // return cuda_api( cudaDeviceSynchronize (  ) );
}

__device__ float rand_dev(  )
{
    // const auto i_th = threadIdx.x + blockDim.x * blockIdx.x;
    // float i_rng = int( float( i_th )  / n_rng_c );
    // i_rng       = i_th - i_rng * n_rng_c;
    // return curand_uniform( rng_c + int( i_rng ) );
    return float( 1 );
}

};                              // namespace device

#endif // __CUDACC__
