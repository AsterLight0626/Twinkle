#pragma once

#ifdef    __CUDACC__
#include "cuda/cuda.h"
using device_t = device::cuda_t;
#elif defined (__HIPCC__)
#include "hip/hip.h"
using device_t = device:: hip_t;
#elif defined (__MUSACC__)
#include "musa/musa.h"
using device_t = device::musa_t;
#endif // __CUDACC__
