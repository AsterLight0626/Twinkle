#include <cstdlib>
#include <cstring>
#include <thread>

#include "device.h"

namespace device
{
////////////////////////////////////////////////////////////
// Con-/destructor

base_t:: base_t(  ) : idx_device( -1 ), num_device( -1 )
{
    f_malloc = [ & ] ( size_t s ) -> void *
               { return std::malloc( s ); };
    f_mset   = [ & ] ( void         * p, const int  & v,
                       const size_t & s, const void * s_ )
               { std::memset( p,  v, s ); };
    f_free   = [ & ] ( void * p )  { std::free( p ); };
    f_cp     = [ & ] ( void * t, const void * s, size_t n )
               { std::memcpy( t, s, n ); };

    f_malloc_host  = f_malloc;
    f_free_host    =   f_free;
    f_cc           =     f_cp;
    max_streams    =     1024;
    idx_streams    =        0;
    head_const     =  nullptr;
    size_const     =        0;
    offset_const   =        0;
    f_malloc_const = [ & ]   ( size_t size ) -> void *
    {
        static const int ali ( 8 );
        const  int s_new( ( size + ali - 1 ) / ali * ali );
        const  int offset = const_offset_get( s_new );
        if( offset + s_new > size_const )
            throw std::runtime_error( "const oversize" );
        void * p = ( ( char * ) head_const ) + offset;
        return p;
    };
    return;
}

////////////////////////////////////////////////////////////
// Initialization

void base_t::init( const int & rank )
{
    // max_streams = args.get< int >
    //     ( "device", "max_streams", max_streams );
    // idx_device  = args.get< int >
    //     ( "device", "idx_device",           -1 );
    num_device  = 0;
    return;
}

////////////////////////////////////////////////////////////
// Callback

void base_t::callback_std( void * p )
{
    return ( * static_cast< func_t * > ( p ) )(  );
}

void base_t::callback_clear( const stream_t & stream )
{
    return cb_map[ stream ].clear(   );
}

void base_t::launch_host( const stream_t & s, func_t && f )
{
    cb_map[ s ].push_back( f );
    auto it( -- cb_map[ s ].end(  ) );
    return f_launch_host( s, callback_std, & ( * it ) );
}

////////////////////////////////////////////////////////////
// Stream and events

void base_t::sync_stream( const stream_t & stream )
{
    if( !  stream_m.empty(  ) )
    {
        auto it  =  stream_m.find( stream );
        if(  it !=  stream_m. end(      ) )
        {
            for( auto & event : it->second )
                event_sync( event );
            it->second.clear(     );
        }
    }
    if( ! sstream_m.empty(  ) )
    {
        auto it  = sstream_m.find( stream );
        if(  it != sstream_m. end(      ) )
        {
            for( auto & stream : it->second )
                sync_stream_base   ( stream );
            it->second.clear(  );
        }
    }
    sync_stream_base     ( stream );
    return callback_clear( stream );
}

void base_t::sync_all_streams(  )
{
    for( auto & [ s, e ] : stream_m )
    {
        for( auto & event  : e )
             event_sync( event );
        e.clear(  );
    }
    for( auto & [ e, s ] :  event_m )
        sync_stream_base ( s );
    event_m.clear(   );
    for( auto & [ s, v ] :   cb_map )
        callback_clear ( s );    
    return;
}

void base_t::stream_record( const stream_t & stream ,
                            const event_t  &  event )
{
    return   stream_m[ stream ].push_back  (  event );
}

void base_t::stream_record( const stream_t & str_w ,
                            const stream_t & str_s )
{
    return  sstream_m [ str_w ]. push_back ( str_s );
}

////////////////////////////////////////////////////////////
// Constant memory pool

int base_t::const_offset_get( const int & add )
{
    offset_const += add ;
    return offset_const - add;
}

};                              // namespace device
