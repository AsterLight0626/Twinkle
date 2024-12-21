#ifdef __MUSACC__

#include <stdio.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>

#include "musa.h"

namespace device
{
////////////////////////////////////////////////////////////
// Con-/destructor

static const size_t const_size( 20  * 1024 );
__constant__ char   const_pool[ const_size ];

musa_t::musa_t(  )
{
    f_malloc       = [ & ]   ( size_t size ) -> void *
    {
        void * p;
        musaMalloc     ( & p,        size );
        musaCheckErrors( "musaMalloc"     );
        return p;
    };
    f_malloc_host  = [ & ]   ( size_t size ) -> void *
    {
        void * p;
        musaMallocHost ( & p,        size );
        musaCheckErrors( "musaMallocHost" );
        return p;
    };
    f_free           = [ & ] ( void * p )
    {
        musaFree       (              p );
        musaCheckErrors( "musaFree"     );
    };
    f_free_host      = [ & ] ( void * p )
    {
        musaFreeHost   (              p );
        musaCheckErrors( "musaFreeHost" );
    };
    f_cp = [ & ] ( void * t, const void * s, size_t n )
    {
        const auto type( t > ( void * )( 0x8000000000 ) &&
                         t < ( void * )( 0x100000000000 )
                         ? musaMemcpyHostToDevice
                         : musaMemcpyDeviceToHost );
        musaMemcpy     ( t, s, n, type );
        musaCheckErrors(  "musaMemcpy" );
    };
    a_cp = [ & ] ( void * t, const void * s, size_t n ,
                   const stream_t & strm )
    {
        const auto type( t > ( void * )( 0x8000000000 ) &&
                         t < ( void * )( 0x100000000000 )
                         ? musaMemcpyHostToDevice
                         : musaMemcpyDeviceToHost );        
        musaMemcpyAsync( t, s, n, type, strm );
        musaCheckErrors( "musaMemcpyAsync"   );
    };
    f_cc    = [ & ] ( void * t, const void * s, size_t n )
    {
        auto di( intptr_t( t ) - intptr_t( head_const ) );
        musaMemcpyToSymbol( const_pool, s, n, int( di ) );
        musaCheckErrors   ( "musaMemcpyToSymbol"        );
    };
    f_const = [ & ] ( void * t, const void * s, size_t n )
    {
        musaMemcpyToSymbol( ( const void * ) t, s, n, 0,
                               musaMemcpyHostToDevice );
        musaCheckErrors   (   "musaMemcpyToSymbol"    );
    };    
    f_launch = [ & ]( const void * ker,  const dim3 & n_bl ,
                      const dim3 & n_th, const int  & s_sh ,
                      const void * strm, void  **     args )
    {
        const auto stream = ( musaStream_t ) ( strm );
        musaLaunchKernel
              ( ker, n_bl, n_th, args, s_sh, stream );
        musaCheckErrors(         "musaLaunchKernel" );
    };
    f_mset = [ & ]( void         *    p, const int  &  val ,
                    const size_t & size, const void * strm )
    {
        musaMemsetAsync( p, val, size, stream_t( strm ) );
        musaCheckErrors(              "musaMemsetAsync" );
    };
    max_streams =       2048;
    idx_streams =          0;
    return;
}

////////////////////////////////////////////////////////////
// Initialization

void musa_t::init( const int & rank )
{
    base_t::init ( args );
    if( idx_device <  0 )
    {
        musaGetDeviceCount( & num_device );
        idx_device =  rank  % num_device  ;
    }
    this->prepare(      );

    musaDeviceSetCacheConfig
        ( musaFuncCachePreferShared );
    return;
}

void musa_t::prepare (  )
{
    musaSetDevice        (               idx_device );
    musaGetSymbolAddress ( & head_const, const_pool );
    musaCheckErrors      ( "Preparations of device" );
    size_const           = const_size;
    return;
}

void musa_t::finalize(  )
{
    delete_all_streams(  );
    return;
}

////////////////////////////////////////////////////////////
// Streams

void musa_t::  pin( void * p, const size_t & s )
{
    musaHostRegister( p, s, musaHostRegisterDefault );
    musaCheckErrors (      "musaHostRegister"        );
}

void musa_t::unpin( void * p )
{
    musaHostUnregister( p );
    musaCheckErrors   ( "musaHostUnregister" );
}

stream_t musa_t::yield_stream( const int & pri )
{
    if( streams.size(  ) < max_streams )
    {
        streams.emplace_back(  );
        musaStreamCreateWithPriority
        ( & streams.back(  ), musaStreamNonBlocking, pri );
        musaCheckErrors( "musaStreamCreate" );
    }
    auto    res  = streams[ idx_streams ];
    idx_streams += 1 ;
    idx_streams  = idx_streams % max_streams;
    return  res ;
}

void musa_t::sync_stream_base( const stream_t & stream )
{
    musaStreamSynchronize            ( stream );
    musaCheckErrors ( "musaStreamSynchronize" );
}

void musa_t::sync_all_streams(  )
{
    musaDeviceSynchronize(  );
    musaCheckErrors( "musaDeviceSynchronize" );
    return base_t::sync_all_streams(  );
}

void musa_t::delete_all_streams(  )
{
    sync_all_streams(        );
    for( auto & p :  streams )
        musaStreamDestroy( p );
    return streams. clear(   );
}

////////////////////////////////////////////////////////////
// Events

event_t musa_t::event_create(  )
{
    event_t e;
    musaEventCreateWithFlags( & e, musaEventDisableTiming );
    musaCheckErrors( "musaEventCreate" );
    return  e;
}

void musa_t::event_destroy( const event_t & event )
{
    musaEventDestroy(      event         );
    musaCheckErrors ( "musaEventDestroy" );
    return;
}

void musa_t::event_record ( const  event_t &  event ,
                            const stream_t & stream )
{
    musaEventRecord( event,     stream );
    musaCheckErrors( "musaEventRecord" );
    return;
}

void musa_t::event_wait   ( const stream_t & stream ,
                            const  event_t &  event )
{
    musaStreamWaitEvent( stream,    event   );
    musaCheckErrors( "musaStreamWaitEvent" );
    return;
}

void musa_t::event_sync   ( const event_t & event )
{
    musaEventSynchronize(      event             );
    musaCheckErrors     ( "musaEventSynchronize" );
    return;
}

////////////////////////////////////////////////////////////
// Callbacks

void musa_t::f_launch_host
( const stream_t & s, const f_cb_t f, void * p )
{
    musaLaunchHostFunc( s, f, p );
    musaCheckErrors   ( "musaLaunchHostFunc" );
    return;
}

};                              // namespace device

#endif // __MUSACC__
