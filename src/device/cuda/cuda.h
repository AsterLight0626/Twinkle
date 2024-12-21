#pragma once

#include <curand_kernel.h>
#include "../device.h"

namespace device
{
////////////////////////////////////////////////////////////
// CUDA

class cuda_t : public base_t
{
    ////////// Con-/destructor //////////
public:
    cuda_t(  );

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

    ////////// RNG //////////
public:                         // Type
    using rand_t = curandStateXORWOW_t;
protected:                      // Data
    rand_t * rng_hd;
public:
    virtual void prep_rng( const size_t & n_rng = 1024 );
};

__device__ float rand_dev(  );
};                              // namespace device
