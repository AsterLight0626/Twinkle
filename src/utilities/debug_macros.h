#pragma once

////////////////////////////////////////////////////////////
// Macros

#define __check_size__(SIZE,LIM,IERR)                     \
    do{                                                   \
        if( unsigned( LIM ) < unsigned( SIZE ) )          \
        {                                                 \
            std::cerr << "[ "#LIM" = " << LIM << " ] < [ "\
                      << #SIZE" = " << SIZE << " ]"       \
                      << std::endl;                       \
            IERR = true;                                  \
        }                                                 \
    }while( 0 );

#define PRINT_ARR(A,NNN)                        \
    do{                                         \
        printf( #A":\n" );                      \
        for( int III = 0; III < NNN; ++ III )   \
            printf( "%+0.7e ", A[ III ] );      \
        printf( "\n" );                         \
    }while( 0 );

#define PRINT_ARR_SP(A,NNN,MMM)                 \
    do{                                         \
        printf( #A":" );                        \
        for( int II = 0; II < NNN; ++ II )      \
        {                                       \
            printf( "\nline %d: ", II );        \
            for( int JJ = 0; JJ < MMM; ++ JJ )  \
                printf( "%+0.7e ",              \
                        A[ II * MMM + JJ ] );   \
        }                                       \
        printf( "\n" );                         \
    }while( 0 );

#define PRINT_IARR(A,NNN)                       \
    do{                                         \
        printf( #A":\n" );                      \
        for( int III = 0; III < NNN; ++ III )   \
            printf( "%d ", A[ III ] );          \
        printf( "\n" );                         \
    }while( 0 );

#define PRINT_IAR3(A)                           \
    printf( #A": [ %d %d %d ]\n\n", A[ 0 ],     \
            A[ 1 ], A[ 2 ] );                   

#define PRINT_AR3(A)                            \
    printf( #A": [ %g %g %g ]\n\n", A[ 0 ],     \
            A[ 1 ], A[ 2 ] );                   

#define PRINT_DIM3(A)                                   \
    do{                                                 \
        const int a[ 3 ] = { int( A.x ), int( A.y ),    \
                             int( A.z ) };              \
        printf( #A": [ %d %d %d ]\n\n", a[ 0 ],         \
                a[ 1 ], a[ 2 ] );                       \
    }while( false );

#define PRINT_MAT(A,MMM,NNN)                            \
do{                                                     \
    printf( #A":\n[" );                                 \
    for( int III = 0; III < MMM; ++ III )               \
    {                                                   \
        printf( "[" );                                  \
        for( int JJJ = 0; JJJ < NNN; ++ JJJ )           \
            printf( "%+0.6e, ", A( III, JJJ ) );        \
        printf( "],\n" );                               \
    }                                                   \
    printf( "]\n\n" );                                  \
}while( 0 );

#define COMPARE(A,B,IERR)                                 \
    do{                                                   \
        if( unsigned( A ) < unsigned( B ) )               \
        {                                                 \
            std::cerr << "[ "#A" = " << A << " ] < [ "    \
                      << #B" = " << B << " ]"             \
                      << std::endl;                       \
            IERR = true;                                  \
        }                                                 \
    }while( 0 );
    
#define TST_PRINT( A )     __syncthreads(  );             \
    if( threadIdx.x == 0 && blockIdx.x == 0 &&            \
         blockIdx.y == 0 && blockIdx.z == 0 && A )
#define TST_ERR throw std::runtime_error( "TST_ERR" );
#define TST std::cerr << __FILE__ << ": " \
                      << __LINE__ << '\n' ;
#define __TEST__
#define __O0 __attribute__((optimize("O0")))
