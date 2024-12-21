#pragma once

////////////////////////////////////////////////////////////
//
#define musaCheckErrors(msg)                               \
    do                                                     \
    {                                                      \
        musaError_t __err = musaGetLastError(  );          \
        if( __err != musaSuccess )                         \
        {                                                  \
            printf( "Fatal error: %s (%s at %s:%d)\n",     \
                    msg, musaGetErrorString( __err ),      \
                    __FILE__, __LINE__ );                  \
            printf( "*** FAILED - ABORTING\n" );           \
            throw std::runtime_error( "musaCheckErrors" ); \
        }                                                  \
    } while( 0 );
