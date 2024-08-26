#pragma once

#include <type_traits>
#include <iostream>
#include <cfloat>
#include <cmath>

////////////////////////////////////////////////////////////
// Precision of floating points

#ifndef    PRECISION
#define    PRECISION 1
#endif  // PRECISION

namespace type
{
#if   PRECISION == 0
typedef float   float_t;
typedef float  float2_t;
#elif PRECISION == 1
typedef float   float_t;
typedef double float2_t;
#elif PRECISION == 2
typedef double  float_t;
typedef double float2_t;
#endif

typedef int             int_t;
typedef unsigned int   uint_t;
typedef bool           bool_t;
typedef signed   char schar_t;
typedef unsigned char uchar_t;

static const int size_f ( sizeof(  float_t ) );
static const int size_f2( sizeof( float2_t ) );
static const int size_i ( sizeof(    int_t ) );
};

////////////////////////////////////////////////////////////
// Vector types

namespace type
{

template< typename T >
struct v3_base_t
{
    T x[ 3 ];

    __host__ __device__ __forceinline__
    static constexpr int size(  ) {  return 3; };

    template < typename  ...  L > __host__ __device__
    v3_base_t( const L & ... ls ) : x { ( T ) ls ... } {  };

    __host__ __device__
    v3_base_t( const T & t ) : x { t, t, t } {  };

    __host__ __device__ __forceinline__
    T       & operator [  ] ( const int_t   & i )
    {
        return x[ i ];
    };
    __host__ __device__ __forceinline__
    const T & operator [  ] ( const int_t   & i ) const
    {
        return x[ i ];
    };
    __host__ __device__ __forceinline__
    bool operator != ( const v3_base_t< T > & r ) const
    {
        bool res( false );
#pragma unroll
        for( int i = 0; i < 3; ++ i )
            res |= ( x[ i ] != r.x[ i ] );
        return res;
    };
    __host__ __device__ __forceinline__
    bool operator == ( const v3_base_t< T > & r ) const
    {
        return ! ( this->operator != ( r ) );
    };
    __host__ __device__ __forceinline__ T prod(  ) const
    {
        T res( 1 );
        for( int i = 0; i < 3; ++ i )
            res *= x[ i ];
        return res;
    };
    __host__ __device__ __forceinline__ T norm(  ) const
    {
        T res( 0 );
        for( int i = 0; i < 3; ++ i )
            res += x[ i ] * x[ i ];
        return sqrtf( res );
    };
    __host__ static constexpr v3_base_t< T > null(  )
    {
        return v3_base_t< T >
             ( { ( T ) 0, ( T ) 0, ( T ) 0 } );
    };

    template< typename F, std::enable_if_t
            < std::is_fundamental< F >::value, int > = 0 >
    __host__ __device__ __forceinline__
    v3_base_t< T > & operator = ( const F & f )
    {
#pragma unroll
        for( int a = 0; a < 3; ++ a )
            x[ a ] = f;
        return ( * this );
    };
};

typedef v3_base_t<    int_t >    idx_t;
typedef v3_base_t<   bool_t >   bvec_t;
typedef v3_base_t<  float_t >  coord_t;
typedef v3_base_t< float2_t > coord2_t;

};
