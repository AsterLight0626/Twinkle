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

static const size_t const_size  ( 56  * 1024 );
__constant__ char   const_pool  [ const_size ];

cuda_t::cuda_t(  )
{
    f_malloc       = [ & ]   ( size_t size ) -> void *
    {
        void * p;
        cudaMalloc     ( & p,        size );
        cudaCheckErrors( "cudaMalloc"     );
        return p;
    };
    f_malloc_host  = [ & ]   ( size_t size ) -> void *
    {
        void * p;
        cudaMallocHost ( & p,        size );
        cudaCheckErrors( "cudaMallocHost" );
        return p;
    };
    f_free           = [ & ] ( void * p )
    {
        cudaFree       (              p );
        cudaCheckErrors( "cudaFree"     );
    };
    f_free_host      = [ & ] ( void * p )
    {
        cudaFreeHost   (              p );
        cudaCheckErrors( "cudaFreeHost" );
    };
    f_cp = [ & ] ( void * t, const void * s, size_t n )
    {
        cudaMemcpy( t, s, n, cudaMemcpyDefault );
        cudaCheckErrors(    "cudaMemcpy"       );
    };
    a_cp = [ & ] ( void * t, const void * s, size_t n ,
                   const stream_t & strm )
    {
        cudaMemcpyAsync( t, s, n, cudaMemcpyDefault, strm );
        cudaCheckErrors(         "cudaMemcpyAsync"        );
    };
    f_cc    = [ & ] ( void * t, const void * s, size_t n )
    {
        auto di( intptr_t( t ) - intptr_t( head_const ) );
        cudaMemcpyToSymbol( const_pool, s, n, int( di ) );
        cudaCheckErrors   ( "cudaMemcpyToSymbol"        );
    };
    f_const = [ & ] ( void * t, const void * s, size_t n )
    {
        cudaMemcpyToSymbol( ( const void * ) t, s, n, 0,
                               cudaMemcpyHostToDevice );
        cudaCheckErrors   (   "cudaMemcpyToSymbol"    );
    };    
    f_launch = [ & ]( const void * ker,  const dim3 & n_bl ,
                      const dim3 & n_th, const int  & s_sh ,
                      const void * strm, void  **     args )
    {
        const auto stream = ( cudaStream_t ) ( strm );
        cudaLaunchKernel
              ( ker, n_bl, n_th, args, s_sh, stream );
        cudaCheckErrors(         "cudaLaunchKernel" );
    };
    f_mset = [ & ]( void         *    p, const int  &  val ,
                    const size_t & size, const void * strm )
    {
        cudaMemsetAsync( p, val, size, stream_t( strm ) );
        cudaCheckErrors(              "cudaMemsetAsync" );
    };
    max_streams =       2048;
    idx_streams =          0;
    return;
}

////////////////////////////////////////////////////////////
// Initialization

void cuda_t::init( const int & rank )
{
    base_t::init ( rank );
    if( idx_device <  0 )
    {
        cudaGetDeviceCount( & num_device );
        idx_device =  rank  % num_device  ;
    }
    this->prepare(      );

    cudaDeviceSetSharedMemConfig
            ( cudaSharedMemBankSizeEightByte );
    cudaDeviceSetCacheConfig
            ( cudaFuncCachePreferShared );
    return;
}

void cuda_t::prepare (  )
{
    cudaSetDevice        (               idx_device );
    cudaGetSymbolAddress ( & head_const, const_pool );
    cudaCheckErrors      ( "Preparations of device" );
    size_const           = const_size;
    return;
}

void cuda_t::finalize(  )
{
    delete_all_streams(  );
    return;
}

////////////////////////////////////////////////////////////
// Streams

void cuda_t::  pin( void * p, const size_t & s )
{
    cudaHostRegister( p, s, cudaHostRegisterDefault );
    cudaCheckErrors (      "cudaHostRegister"        );
}

void cuda_t::unpin( void * p )
{
    cudaHostUnregister( p );
    cudaCheckErrors   ( "cudaHostUnregister" );
}

stream_t cuda_t::yield_streams( const int & pri )
{
    if( streams.size(  ) < max_streams )
    {
        streams.emplace_back(  );
        cudaStreamCreateWithPriority
        ( & streams.back(  ), cudaStreamNonBlocking, pri );
        cudaCheckErrors( "cudaStreamCreate" );
    }
    auto    res  = streams[ idx_streams ];
    idx_streams += 1 ;
    idx_streams  = idx_streams % max_streams;
    return  res ;
}

void cuda_t::sync_stream_base( const stream_t & stream )
{
    cudaStreamSynchronize            ( stream );
    cudaCheckErrors ( "cudaStreamSynchronize" );
}

void cuda_t::sync_all_streams(  )
{
    cudaDeviceSynchronize(  );
    cudaCheckErrors( "cudaDeviceSynchronize" );
    return base_t::sync_all_streams(  );
}

void cuda_t::delete_all_streams(  )
{
    sync_all_streams(        );
    for( auto & p :  streams )
        cudaStreamDestroy( p );
    return streams. clear(   );
}

////////////////////////////////////////////////////////////
// Events

event_t cuda_t::event_create(  )
{
    event_t e;
    cudaEventCreateWithFlags( & e, cudaEventDisableTiming );
    cudaCheckErrors( "cudaEventCreate" );
    return  e;
}

void cuda_t::event_destroy( const event_t & event )
{
    cudaEventDestroy(      event         );
    cudaCheckErrors ( "cudaEventDestroy" );
    return;
}

void cuda_t::event_record ( const  event_t &  event ,
                            const stream_t & stream )
{
    cudaEventRecord( event,     stream );
    cudaCheckErrors( "cudaEventRecord" );
    return;
}

void cuda_t::event_wait   ( const stream_t & stream ,
                            const  event_t &  event )
{
    cudaStreamWaitEvent( stream,    event   );
    cudaCheckErrors( "cudaStreamWaitEvent" );
    return;
}

void cuda_t::event_sync   ( const event_t & event )
{
    cudaEventSynchronize(      event             );
    cudaCheckErrors     ( "cudaEventSynchronize" );
    return;
}

////////////////////////////////////////////////////////////
// Callbacks

void cuda_t::f_launch_host
( const stream_t & s, const f_cb_t f, void * p )
{
    cudaLaunchHostFunc( s, f, p );
    cudaCheckErrors   ( "cudaLaunchHostFunc" );
    return;
}

};                              // namespace device

#endif // __CUDACC__
