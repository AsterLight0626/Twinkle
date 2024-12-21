#pragma once

////////////////////////////////////////////////////////////
//

#if   defined (__CUDACC__)

#define __dyn_shared__(__TYPE__, __NAME__) \
 extern     __shared__ __TYPE__  __NAME__  [  ]

namespace device
{
using stream_t = cudaStream_t;
using  event_t =  cudaEvent_t;
};                              // namespace device

#elif defined (__MUSACC__)

#define __dyn_shared__(__TYPE__, __NAME__) \
 extern     __shared__ __TYPE__  __NAME__  [  ]

namespace device
{
using stream_t = musaStream_t;
using  event_t =  musaEvent_t;
};                              // namespace device

#elif defined (__HIPCC__)

#include <cmath>
#include <hip/hip_runtime.h>
#include <hip/hip_runtime_api.h>

namespace device
{
using stream_t =  hipStream_t;
using  event_t =   hipEvent_t;
};                              // namespace device

#define __dyn_shared__ HIP_DYNAMIC_SHARED

#ifdef __HIP_CPU_RT__

__device__ static inline int __float_as_int( float x )
{
    static_assert( sizeof( int ) == sizeof( float ), "" );
    int tmp;
    __builtin_memcpy( & tmp, & x, sizeof( tmp ) );
    return tmp;
}

__device__ static inline float __int_as_float( int x )
{
    static_assert( sizeof( float ) == sizeof( int ), "" );
    float tmp;
    __builtin_memcpy( & tmp, & x, sizeof( tmp ) );
    return tmp;
}

template  < class f_T > __device__ __forceinline__
bool isinf( const f_T & f )
{
    return std::isinf ( f );
}

template  < class f_T > __device__ __forceinline__
bool isnan( const f_T & f )
{
    return std::isnan ( f );
}

__device__ static inline long long int
__double_as_longlong( double x )
{
    static_assert
        ( sizeof( long long ) == sizeof( double ), "" );
    long long tmp;
    __builtin_memcpy( & tmp, & x, sizeof( tmp ) );
    return tmp;
}

#undef  __forceinline__
#define __forceinline__ inline

#endif // __HIP_CPU_RT__

#endif // __CUDACC__ v.s. __HIPCC__
