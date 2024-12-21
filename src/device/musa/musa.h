#pragma once

#include "../device.h"
#include "musa_utils.h"

namespace device
{
////////////////////////////////////////////////////////////
// MUSA

class musa_t : public base_t
{
    ////////// Con-/destructor //////////
public:
    musa_t(  );

    ////////// Initialization //////////
public:                         // Functions
    virtual void prepare (  );
    virtual void finalize(  );
    virtual void init    ( const int & rank = 0 );

    ////////// Memory //////////
public:                         // Functions
    virtual void   pin( void * p, const size_t & s );
    virtual void unpin( void * p                   );

    ////////// Streams //////////
public:
    virtual stream_t  yield_stream ( const int  & p = 0 );
    virtual void delete_all_streams(                    );
    virtual void   sync_all_streams(                    );
    virtual void   sync_stream_base( const stream_t & s );

    ////////// Events //////////
public:                         // Functions
    virtual event_t event_create (  );
    virtual void    event_destroy( const event_t  & e );
    virtual void    event_record ( const event_t  & e ,
                                   const stream_t & s );
    virtual void    event_wait   ( const stream_t & s ,
                                   const event_t  & e );
    virtual void    event_sync   ( const event_t  & e );

    ////////// Callback //////////
public:                         // Functions
    virtual void f_launch_host
    ( const stream_t & s, const f_cb_t f, void * p );
};
};                              // namespace device
