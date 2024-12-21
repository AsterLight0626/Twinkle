#pragma once

#include "device_macros.h"

#include <iostream>
#include <functional>
#include <type_traits>
#include <vector>
#include <atomic>
#include <list>
#include <map>

namespace device
{
////////////////////////////////////////////////////////////
// Default base class: CPU-side operations

class base_t
{
    ////////// Con-/destructor //////////
public:                         // Functions
     base_t(  ) ;
    ~base_t(  ) {  };

    ////////// Initialization //////////
protected:                      // Data
    int            idx_device;
    int            num_device;
public:                         // Functions
    virtual void prepare (  ) {  };
    virtual void finalize(  ) {  };
    virtual void init    ( const int & rank = 0 );

    ////////// Malloc/free //////////
protected:                      // Types
    typedef std::function
    < void( void         *, const  int  & ,
            const size_t &, const  void * ) > f_memset_t;
    typedef std::function< void *( size_t ) > f_malloc_t;
    typedef std::function< void  ( void * ) >   f_free_t;
public:                         // Modules
    f_malloc_t f_malloc;
    f_malloc_t f_malloc_host;
    f_malloc_t f_malloc_const;
    f_memset_t f_mset;
    f_free_t   f_free;
    f_free_t   f_free_host;
public:                         // Functions
    template< typename T  > T * malloc_device
    ( const size_t & n    = 1    )        const;
    template< typename T  > T * malloc_const
    ( const size_t & n    = 1    )        const;    
    template< typename T  > T * malloc_host
    ( const size_t & n    = 1    )        const;
    template< typename T  > T * malloc
    ( const size_t & n    = 1,
      const bool   & host = true )        const;
    template< typename T  >
    void  free_device( T  * & p )         const;
    template< typename T  >
    void  free_host  ( T  * & p )         const;
    template< typename T  > void   free
    ( T * & p, const bool & host = true ) const;
    virtual void    pin( void * p , const size_t & s ) {  };
    virtual void  unpin( void * p )                    {  };

    ////////// Streams //////////
protected:                      // Data
    int                 idx_streams;
    int                 max_streams;
    std::vector< stream_t > streams;
public:
    virtual stream_t  yield_stream
    ( const int & priority = 0 ) = 0;
    virtual void delete_all_streams(                  ) = 0;
    virtual void   sync_all_streams(                  ) ;
    virtual void   sync_stream_base( const stream_t & ) = 0;
    virtual void   sync_stream     ( const stream_t & ) ;

    ////////// Events //////////
protected:                      // Data
    std::map< event_t,               stream_t   >   event_m;
    std::map< stream_t, std::vector<  event_t > >  stream_m;
    std::map< stream_t, std::vector< stream_t > > sstream_m;
public:                         // Functions
    virtual event_t event_create (  ) = 0;
    virtual void    event_destroy( const event_t  & e ) = 0;
    virtual void    event_record ( const event_t  & e ,
                                   const stream_t & s ) = 0;
    virtual void    event_wait   ( const stream_t & s ,
                                   const event_t  & e ) = 0;
    virtual void    event_sync   ( const event_t  & e ) = 0;
    virtual void    stream_record( const stream_t & s ,
                                   const event_t  & e );
    virtual void    stream_record( const stream_t & s ,
                                   const stream_t & r );

    ////////// Copy of memory //////////
protected:                      // Types
    typedef std::function
    < void( void *, const void *, size_t ) >  f_cp_t;
    typedef std::function
    < void( void *, const void *, size_t ,
            const   stream_t           & ) >  a_cp_t;
protected:                      // Data
    void *   head_const;        // __constant__ head loc.
    size_t   size_const;        // __constant__ size max.
    size_t offset_const;        // __constant__ offset
    virtual  int  const_offset_get( const int & );
public:                         // Modules
    f_cp_t f_cp;
    a_cp_t a_cp;
    f_cp_t f_cc;                // copy to __constant__ pool
    f_cp_t f_const;             // copy to __constant__ 
public:                         // Functions
    template< typename T, std::enable_if_t
    < ! std::is_pointer< T >::value, int > = 0 > void cp
    ( T & tgt, const T & src, const size_t & n = 1 ) const;
    template< typename T > void cp
    ( T * tgt, const T * src, const size_t & n = 1 ) const;
    template< typename T > void cp_a
    ( T * tgt, const T * src, const size_t & n = 1,
       const   stream_t & str  = 0 )                 const;
    template< typename T > void cp_const
    ( T & tgt, const T & src, const size_t & n = 1 ) const;

    ////////// Launch functions (kernels) //////////
protected:                      // Types
    using f_lau_t = std::function // Generic device mode
    < void( const void *, const dim3 &,  const dim3 &,
            const int  &, const void *,  void ** ) > ;
    using f_exe_t = std::function // CPU-mode
    < void( const dim3 &, const dim3 &,  const int & ,
            const void *, std::function< void(   ) > ) >;
    using lpar_t = std::tuple< dim3, dim3, int   > ;
protected:
    f_lau_t f_launch;
    f_exe_t f_exec;
protected:                      // Function templates
    template< size_t n, typename ... args_T,
              typename std::enable_if_t
            < n == sizeof ... ( args_T ), int > = 0 >
    void parse_args
    ( void  *                          args_arr ,
      const std::tuple< args_T ... > & args_tup ) const;
    template< size_t n, typename ... args_T,
              typename std::enable_if_t
            < n != sizeof ... ( args_T ), int > = 0 >
    void parse_args
    ( void  **                         args_arr ,
      const std::tuple< args_T ... > & args_tup ) const;
public:                         // Interface
   template< typename ... args_T,
              typename fun_T = void( * ) ( args_T ... ) >
    __forceinline__ void __launch
    ( const fun_T & ker,    const dim3   &     n_bl ,
      const dim3  & n_th,   const int    &     s_sh ,
      const void  * stream, const args_T & ... args ) const;
    template< typename ... args_T,
              typename fun_T = void( * ) ( args_T ... ) >
    __forceinline__ void launch
    ( const fun_T & ker,    const lpar_t &     lpar ,
      const void  * stream, const args_T & ... args ) const;

    ////////// Callbacks on the host side //////////
protected:                      // Types
    using func_t = std::function< void(  ) >;
    using f_cb_t = void ( * )   ( void   * );
protected:                      // Datakj
    std::map< stream_t, std::list< func_t > > cb_map;
protected:                      // Function
    static  void callback_std  ( void           * p ) ;
public:                         // Function    
    virtual void callback_clear( const stream_t & s ) ;
    virtual void  f_launch_host
    ( const stream_t &, const f_cb_t, void * ) = 0;
    virtual void launch_host
    ( const stream_t & s, func_t && f );

    ////////// RNG //////////
public:                         // Functions
    virtual void prep_rng( const size_t & n_rng = 1024 ){ };
};

////////////////////////////////////////////////////////////
// Malloc/free

template< typename T >
T * base_t::malloc_host  ( const size_t & n ) const
{
    return ( T * ) ( f_malloc_host( n * sizeof( T ) ) );
}

template< typename T >
T * base_t::malloc_device( const size_t & n ) const
{
    return ( T * ) ( f_malloc     ( n * sizeof( T ) ) );
}

template< typename T >
T * base_t::malloc_const ( const size_t & n ) const
{
    return ( T * ) ( f_malloc_const( n * sizeof( T ) ) );
}

template< typename T > T * base_t::malloc
( const size_t & n, const bool & host ) const
{
    if( host )
        return malloc_host  < T > ( n );
    else
        return malloc_device< T > ( n );
}

template< typename T >
void base_t::free_host  ( T * & p ) const
{
    f_free_host( p );
    p = nullptr;
    return;
}

template< typename T >
void base_t::free_device( T * & p ) const
{
    f_free     ( p );
    p = nullptr;
    return;
}

template< typename T > void base_t::free
( T * & p, const bool & host ) const
{
    if( host )
        free_host  ( p );
    else
        free_device( p );
    return;
}

////////////////////////////////////////////////////////////
// Memory copy interfaces

template< typename T > void base_t::cp
( T * tgt, const T * src, const size_t & n ) const
{
    f_cp( ( void * )     tgt,   ( const void * )     src,
          n * sizeof( T ) );
    return;
}

template< typename T > void base_t::cp_a
( T * tgt, const T *  src , const size_t & n,
  const stream_t & stream ) const
{
    a_cp( ( void * )     tgt,   ( const void * )     src,
          n * sizeof( T ),              stream );
    return;
}

template< typename T, std::enable_if_t
        < ! std::is_pointer< T >::value, int > >
void base_t::cp
( T & tgt, const T & src, const size_t & n ) const
{
    f_cp( ( void * ) ( & tgt ), ( const void * ) ( & src ),
          n * sizeof( T ) );
    return;
}

template< typename T > void base_t::cp_const
( T & tgt, const T & src, const size_t & n ) const
{
    f_const
    ( ( void * ) ( & tgt ), ( const void * ) ( & src ),
      n * sizeof( T ) );
    return;
};

////////////////////////////////////////////////////////////
// Kernel launcher

template< size_t n, typename ... args_T,
          typename std::enable_if_t
        < n == sizeof ... ( args_T ), int > >
void base_t::parse_args
( void  *                          args_arr,
  const std::tuple< args_T ... > & args_tup ) const
{
    return;
}

template< size_t n, typename ... args_T,
          typename std::enable_if_t
        < n != sizeof ... ( args_T ), int> >
void base_t::parse_args
( void  **                         args_arr,
  const std::tuple< args_T ... > & args_tup ) const
{
    args_arr[ n ] = ( void * )
                  ( & std::get< n >( args_tup ) );
    return parse_args< n + 1 >( args_arr, args_tup );
}

template< typename ... args_T, typename   fun_T >
__forceinline__ void base_t::__launch
( const fun_T & ker,    const dim3   &     n_bl ,
  const dim3  & n_th,   const int    &     s_sh ,
  const void  * stream, const args_T & ... args ) const
{
    if( f_launch )
    {
        void * args_arr[ sizeof ... ( args )  ];        
        parse_args< 0 >( args_arr,
                         std::forward_as_tuple( args... ) );
        f_launch( reinterpret_cast< const void * >  ( ker ),
                  n_bl, n_th, s_sh, stream,      args_arr );
    }
    else if( f_exec )           // CPU mode
        f_exec( n_bl, n_th, s_sh, stream, [ = ] (  )
        {  ker( args ... );  } );
    else
        throw std::runtime_error( "Kernel exec undefined" );
    return;
}

template< typename ... args_T, typename   fun_T >
__forceinline__ void base_t::launch
( const fun_T & ker,    const lpar_t &     lpar ,
  const void  * stream, const args_T & ... args ) const
{
    __launch( ker,                   std::get< 0 >( lpar ),
              std::get< 1 >( lpar ), std::get< 2 >( lpar ),
              stream,  args  ...  );
    return;
}

};                              // namespace device
