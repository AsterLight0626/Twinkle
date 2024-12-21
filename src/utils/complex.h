#pragma once

#include <cmath>

namespace twinkle
{
////////////////////////////////////////////////////////////
//

template < class f_T >
struct complex_t
{
    ////////// Data and types //////////
    using this_t = complex_t< f_T >;
    using    f_t = f_T  ; 
    f_T   re, im ;

    ////////// Constructors //////////
    __host__ __device__ complex_t(  )
              : re( 0 ), im( 0 ) {  };

    template < class r_T, class i_T > __device__ __host__
    complex_t( const r_T & re, const i_T & im )
    {
        this->re = re;
        this->im = im;
        return;
    };
    template < class r_T >
    __device__ __host__ complex_t( const r_T & re )
    {
        this->re = re;
        this->im =  0;
        return;
    };

    ////////// Assign values //////////
    template < class r_T, class i_T > __device__ __host__
    void set( const r_T & re, const i_T & im )
    {
        this->re = re;
        this->im = im;
        return;
    };
    template < class r_T > __device__ __host__
    void  operator = ( const r_T & re )
    {
        this->re = re;
        this->im =  0;
        return;
    };
    __device__ __host__
    void operator = ( const this_t &  s )
    {
        this->re = s.re;
        this->im = s.im;
        return;
    };

    ////////// Uniary operators //////////    
    __device__ __host__ __forceinline__ f_T norm2(  ) const
    {
        return re * re + im * im;
    };
    __device__ __host__ __forceinline__ f_T abs  (  ) const
    {
        return std::sqrt( this->norm2(  ) );
    };
    __device__ __host__ __forceinline__
    complex_t < f_T > conj(  ) const
    {
        return complex_t < f_T > ( this->re, -this->im );
    };
    __device__ __host__ complex_t < f_T > sqrt(  ) const
    {
        complex_t< f_T > res;
        const auto abs_old = this->abs(    );
        if( abs_old <= 0 )
            return res;
        const auto abs_new = std::sqrt( abs_old );
        res.re = abs_new
            * std::sqrt( fabs( ( 1 + re / abs_old ) / 2 ) )
            * ( im >= 0 ? 1 : -1 ) ;
        res.im = abs_new
            * std::sqrt( fabs( ( 1 - re / abs_old ) / 2 ) );
        return res;
    };    

    __device__ __host__ __forceinline__
    void operator += ( const this_t & s )
    {
        this->re += s.re;
        this->im += s.im;
        return;
    };
    template < class r_T > __device__ __host__ 
    void operator += ( const r_T & s )
    {
        this->re += s;
        return;
    };
    __device__ __host__ __forceinline__
    void operator -= ( const this_t & s )
    {
        this->re -= s.re;
        this->im -= s.im;
        return;
    };
    template < class r_T > __device__ __host__ 
    void operator -= ( const r_T & s )
    {
        this->re -= s;
        return;
    };    
    __device__ __host__ __forceinline__
    void operator *= ( const this_t & s )
    {
        this_t previous;
        previous.set( this->re, this->im );
        this->re = previous.re * s.re - previous.im * s.im;
        this->im = previous.re * s.im + previous.im * s.re;
        return;
    };
    template < class r_T > __device__ __host__
    void operator *= ( const r_T & s )
    {
        this->re *= s;
        this->im *= s;        
        return;
    };
    __device__ __host__ __forceinline__
    void operator /= ( const this_t & s )
    {
        this_t previous;
        previous.set( this->re, this->im );
        this->re =   previous.re * s.re + previous.im * s.im;
        this->im = - previous.re * s.im + previous.im * s.re;
        this->re /= s.norm2(  );
        this->im /= s.norm2(  );
        return;
    };
    template < class r_T > __device__ __host__
    void operator /= ( const r_T & s )
    {
        this->re /= s;
        this->im /= s;        
        return;
    };
    __device__ __host__ friend 
    complex_t< f_T > operator -
    ( const complex_t< f_T > & s )
    {
        return complex_t< f_T > ( -s.re, -s.im );
    };

    ////////// Binary operators //////////
    __device__ __host__ __forceinline__ friend 
    complex_t< f_T > operator +
    ( const complex_t< f_T > & l ,  const f_T & r )
    {
        return complex_t< f_T > ( l.re + r, l.im );
    };
    __device__ __host__ __forceinline__ friend 
    complex_t< f_T > operator +
    ( const f_T & l, const complex_t< f_T > & r )
    {
        return r + l;
    };
    __device__ __host__ __forceinline__ friend 
    complex_t< f_T > operator +
    ( const complex_t< f_T > & l,
      const complex_t< f_T > & r )
    {
        return complex_t< f_T >( l.re + r.re, l.im + r.im );
    };

    __device__ __host__ __forceinline__ friend 
    complex_t< f_T > operator -
    ( const complex_t< f_T > & l ,  const f_T & r )
    {
        return complex_t< f_T > ( l.re - r, l.im );        
    };
    __device__ __host__ __forceinline__ friend 
    complex_t< f_T > operator -
    ( const f_T & l, const complex_t< f_T > & r )
    {
        return complex_t< f_T > ( l - r.re, -r.im );
    };
    __device__ __host__ __forceinline__ friend 
    complex_t< f_T > operator -
    ( const complex_t< f_T > & l,
      const complex_t< f_T > & r )
    {
        return complex_t< f_T >( l.re - r.re, l.im - r.im );
    };

    __device__ __host__ __forceinline__ friend 
    complex_t< f_T > operator *
    ( const complex_t< f_T > & l ,  const f_T & r )
    {
        return complex_t< f_T >( l.re * r, l.im * r );
    };    
    __device__ __host__ __forceinline__ friend 
    complex_t< f_T > operator *
    ( const f_T & l, const complex_t< f_T > & r )
    {
        return r * l;
    };
    __device__ __host__ __forceinline__ friend 
    complex_t< f_T > operator *
    ( const complex_t< f_T > & l,
      const complex_t< f_T > & r )
    {
        auto res( l );
        res *=   r;
        return res;
    };

    __device__ __host__ __forceinline__ friend 
    complex_t< f_T > operator /
    ( const complex_t< f_T > & l ,  const f_T & r )
    {
        return complex_t< f_T >( l.re / r, l.im / r );
    };    
    __device__ __host__ __forceinline__ friend 
    complex_t< f_T > operator /
    ( const f_T & l, const complex_t< f_T > & r )
    {
        complex_t< f_T > res ( r.conj(  ) );
        res *= l / r.norm2(  );
        return  res ;
    };
    __device__ __host__ __forceinline__ friend 
    complex_t< f_T > operator /
    ( const complex_t< f_T > & l,
      const complex_t< f_T > & r )
    {
        complex_t< f_T > res ( l );
        res *=   ( f_T ( 1 ) / r );
        return res ;
    };

    __device__ __host__ __forceinline__ friend 
    bool operator == 
    ( const complex_t< f_T > & l ,
      const complex_t< f_T > & r )
    {
        return ( l.re == r.re ) && ( l.im == r.im );
    };
    __device__ __host__ __forceinline__ friend 
    bool operator != 
    ( const complex_t< f_T > & l ,
      const complex_t< f_T > & r )
    {
        return ! ( l == r );
    };
};

};                              // namespace twinkle

