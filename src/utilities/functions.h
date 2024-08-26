#pragma once

#include "types.h"

namespace utils
{
////////////////////////////////////////////////////////////
// For HIP-CPU compatibility

#if defined(__CUDA_ARCH__) || ( ! defined(__HIP_CPU_RT__) )
#define __isnan       isnan
#define __isinf       isinf
#else
#define __isnan std::isnan
#define __isinf std::isinf
#endif // __CUDA_ARCH__ || ( ! __HIP_CPU_RT__)

template< class F > __device__ __host__
F max( const F & a, const F & b )
{
    return ( a > b ? a : b );
}

template< class F > __device__ __host__
F min( const F & a, const F & b )
{
    return ( a < b ? a : b );
}

template< class F > __device__ __host__
F abs( const F & a )
{
    return ( a > 0 ? a : -a );
}

////////////////////////////////////////////////////////////
// Precision-dependent functions

//////////////////////////////////////////////////
#if PRECISION == 2

static const type::float_t  __epsilon_fl( DBL_EPSILON );

#ifdef __HIP_CPU_RT__
inline
#else
__forceinline__
#endif
__device__ __host__ double dmax( double a, double b )
{
    return ( a > b ? a : b );
}

#ifdef __HIP_CPU_RT__
inline
#else
__forceinline__
#endif
__device__ __host__ double dmin( double a, double b )
{
    return ( a < b ? a : b );
}

#define __max1   utils::dmax
#define __min1   utils::dmin
#define __pow1   pow
#define __exp1   exp
#define __abs1  fabs

#if defined(__CUDA_ARCH__) || ( ! defined(__HIP_CPU_RT__) )
#define __sqrt1 __dsqrt_rn
#else
#define __sqrt1 sqrt
#endif  // __CUDA_ARCH__ || ( ! __HIP_CPU_RT__ )

#else   // PRECISION == 0 || PRECISION == 1

static const type::float_t __epsilon_fl( FLT_EPSILON );

#define __max1  fmaxf
#define __min1  fminf

#if defined( __MUSACC__ )
#define __pow1  __powf
#define __exp1  __expf
#define __sqrt1  sqrtf
#define __abs1    fabs
#elif defined(__CUDA_ARCH__) || ( !defined(__HIP_CPU_RT__) )
#define __pow1  __powf
#define __exp1  __expf
#define __sqrt1 __fsqrt_rn
#define __abs1    fabs
#else
#define __sqrt1 sqrtf
#define __exp1   expf
#define __pow1   powf
#define __abs1   fabs
#endif  // __CUDA_ARCH__ || ( ! __HIP_CPU_RT__ )

#endif  // PRECISION == 2

//////////////////////////////////////////////////
#if PRECISION == 0

static const type::float2_t __epsilon_f2( FLT_EPSILON );

#define __max2  fmaxf
#define __min2  fminf
#define __abs2  fabsf

#if defined(__CUDA_ARCH__) || ( ! defined(__HIP_CPU_RT__) )
#define __pow2  __powf
#define __sqrt2 __fsqrt_rn
#else
#define __pow2   powf
#define __sqrt2 sqrtf
#endif  // __CUDA_ARCH__ || ( ! __HIP_CPU_RT__ )

#else   // PRECISION != 0

static const type::float2_t __epsilon_f2( DBL_EPSILON );

#define __max2  fmax
#define __min2  fmin
#define __pow2   pow
#define __abs2  fabs

#if defined(__CUDA_ARCH__) || ( ! defined(__HIP_CPU_RT__) )
#define __sqrt2 __fsqrt_rn
#else
#define __sqrt2 sqrt
#endif  // __CUDA_ARCH__ || ( ! __HIP_CPU_RT__ )

#endif  // PRECISION == 0

////////////////////////////////////////////////////////////
// Functions

template< class i_T >__device__ __forceinline__
int idiv( const i_T & n, const i_T & m )
{
#ifdef __MUSACC__
    if( n < 1e5f )
        return type::float_t( n ) / m + 1e-5f;
#endif    
    return ( n / m );
};

template< class f_T > __forceinline__ __host__ __device__
bool neg_prod( const f_T & a, const f_T & b )
{
    return ( a < 0 ) ^ ( b < 0 );
}

template< class f_T > __forceinline__ __host__ __device__
const f_T minmod( const f_T & a, const f_T & b )
{
    return ( neg_prod( a, b ) ? 0 :
             ( __abs1 ( a ) < __abs1( b ) ? a : b ) );
}

template< class f_T > __forceinline__ __host__ __device__
const f_T maxmod( const f_T & a, const f_T & b )
{
    return ( neg_prod( a, b ) ? 0 :
             ( __abs1 ( a ) > __abs1( b ) ? a : b ) );
}

template< class f_T > __forceinline__ __host__ __device__
const f_T superbee( const f_T & a, const f_T & b )
{
    if( neg_prod( a, b ) )
        return 0;
    const auto & a_2 ( a * 2 );
    const auto & b_2 ( b * 2 );
    const f_T  & minmod_2a = __abs1( a_2 ) < __abs1( b )
                           ? a_2 : b;
    const f_T  & minmod_2b = __abs1( b_2 ) < __abs1( a )
                           ? b_2 : a;
    return __abs1( minmod_2a ) > __abs1( minmod_2b ) ?
                   minmod_2a   :         minmod_2b ;
}

template< class T > __forceinline__ __host__ __device__
int  sgn( const T & a )
{
    return ( a > 0 ? 1 : ( a < 0 ? -1 : 0 ) );
}

template< class U, class V, class W >
__forceinline__ __host__ __device__
void cross( U * x, const V & a, const W & b )
{
    x[ 0 ] = a[ 1 ] * b[ 2 ] - a[ 2 ] * b[ 1 ];
    x[ 1 ] = a[ 2 ] * b[ 0 ] - a[ 0 ] * b[ 2 ];
    x[ 2 ] = a[ 0 ] * b[ 1 ] - a[ 1 ] * b[ 0 ];    
    return;
}

template< class I , class J > __forceinline__ __host__
__device__ I shift( const I & i, const J & j )
{
    return j >= 0 ? ( i << j ) : ( i >> ( -j ) );
}
};                              // namespace utils
