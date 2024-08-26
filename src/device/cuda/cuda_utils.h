#pragma once

////////////////////////////////////////////////////////////
//
#define cudaCheckErrors(msg)                               \
    do                                                     \
    {                                                      \
        cudaError_t __err = cudaGetLastError(  );          \
        if( __err != cudaSuccess )                         \
        {                                                  \
            printf( "Fatal error: %s (%s at %s:%d)\n",     \
                    msg, cudaGetErrorString( __err ),      \
                    __FILE__, __LINE__ );                  \
            printf( "*** FAILED - ABORTING\n" );           \
            throw std::runtime_error( "cudaCheckErrors" ); \
        }                                                  \
    } while( 0 );
