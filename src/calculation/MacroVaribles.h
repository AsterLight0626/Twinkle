#pragma once

#define ONGPU
// CPU/GPU调试用
#ifdef ONGPU
#define _Host __host__
#define _Device __device__
#define _Host_Device __host__ __device__ 
#define _Global __global__
#else
#define _Device
#define _Host_Device
#define _Global
#define __host__
#define __device__
#endif

// NPOINTS >= MAX_ADAPTIVE
// #define NBODY 2
// #define LENGTH NBODY*NBODY+2
// #define NPOINTS 1024
// #define NSRCS 10240
// #define PI 3.141592653589793
// #define SOLVEMAX 200
// #define POINTS_PER_BLOCK 64
// #define MAX_ADAPTIVE 1024
// #define NEWTON_STEPS 10
// #define MCMC_STEPS 1000
// #define NPARAMS 7
// #define RAND_SEED 212626
// #define NWALKERS 1
// #define THRESHOLD_REAL_DOUBLE 1e-5
// #define THRESHOLD_REAL_FLOAT 3e-3

static const int BATCH_SIZE = 64;            // must be the power of 2, like 32, 64, 128...           // 边缘取点，每次取 BATCH_SIZE 个；BATCH_SIZE <= NPOINTS, NPOINTS % BATCH_SIZE = 0  
static const float RELTOL = 1e-7;            // relative tolerance
static const int NCRIT = 256;               // must be 2^n and < 1024
static const int NCROSS_CANDIDATE = 128;     // needn't to change unless warning information let you do so
static const int NPOINTS = 4096;             // maximum sampling points on the edge, set memory space. max = 6144 (shared memory)
static const int NSRCS = 1024;                // set memory space, must >= Nsrcs in Twinkle

#define DetailedRecord false