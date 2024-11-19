#pragma once

#include<iostream>
// #include<concepts>
#include"MacroVaribles.h"
#include"../utilities/functions.h"


static const int NBODY=2;
static const int LENGTH=NBODY*NBODY+2;
static const double PI = 3.141592653589793;


// // 概念：复数
// template <typename f_T>
// concept isComplex = requires(){f_T::re;f_T::im;};    // 有re和im这两个属性 // requires可以看作一个type的函数，只要没报错，就通过
// 概念：实数（其实是float和double）
// template <typename f_T>
// concept isFloating = std::is_floating_point_v<f_T>;
// // 概念：任意
// template <typename T>
// concept Any = true;


#define isFloating typename


// 结构体：复数

template < isFloating f_T >
struct complex_t {
    f_T re, im;
    
    template< class r_T, class i_T >
    _Host_Device complex_t( const r_T & re_, const i_T & im_ )  : re(re_), im(im_)   {};
    template< class r_T >
    _Host_Device complex_t( const r_T & re_ ) : re( re_ ), im( 0 ) {};
    // _Host_Device complex_t()                 : complex_t(0. , 0.) {};

    // _Host_Device complex_t(f_T _re, f_T _im);
    // _Host_Device complex_t(f_T _re);
    _Host_Device complex_t() {  };

    // // 转换构造，本质上还是构造函数
    // template<isFloating f2_T>
    // _Host_Device complex_t(const complex_t<f2_T> & z);
    // {
    //     this->re = (f_T) z.re;
    //     this->im = (f_T) z.im;
    // } 

    template <isFloating f2_T>
    _Host_Device void operator=(const complex_t<f2_T> & z);
    // {
    //     this->re = (f_T) z.re;      // 转化为f_T类型
    //     this->im = (f_T) z.im;
    // };

    // template <isFloating f2_T>
    // _Host_Device void operator=(const f2_T & z);
    // {
    //     this->re = (f_T) z;
    //     this->im = (f_T) 0;
    // };
};

// 函数声明：复数运算

// Unary

template <isFloating f_T>
_Host_Device f_T norm(const complex_t<f_T> &z);

// real 与 imag 直接通过 .re .im访问
// 好吧，也写一份
template <isFloating f_T>
_Host_Device f_T real(const complex_t<f_T> &z);

template <isFloating f_T>
_Host_Device f_T imag(const complex_t<f_T> &z);

template <isFloating f_T>
_Host_Device complex_t<f_T> conj(const complex_t<f_T> &z);

// 负号
template<isFloating f_T>
_Host_Device complex_t<f_T> operator-(const complex_t<f_T> z);

// _Host float abs_c(const complex_t<float> &z);
// _Host double abs_c(const complex_t<double> &z);
_Device float abs_c(const complex_t<float> &z);
_Device double abs_c(const complex_t<double> &z);

// _Host complex_t<float> sqrt(const complex_t<float> &x);
// _Host complex_t<double> sqrt(const complex_t<double> &x);
_Device complex_t<float> sqrt(const complex_t<float> &x);
_Device complex_t<double> sqrt(const complex_t<double> &x);

// 支持cout输出
template <isFloating f_T>
std::ostream & operator<<(std::ostream &out, const complex_t<f_T> &A);


template <isFloating f_T>
_Host_Device bool operator==(const complex_t<f_T> &z1, const complex_t<f_T> &z2);

template <isFloating f_T>
_Host_Device bool operator!=(const complex_t<f_T> &z1, const complex_t<f_T> &z2);

// same type
template <isFloating f_T>
_Host_Device complex_t<f_T> operator+(const complex_t<f_T> &z1,const complex_t<f_T> &z2);

template <isFloating f_T>
_Host_Device complex_t<f_T> operator-(const complex_t<f_T> &z1,const complex_t<f_T> &z2);

template <isFloating f_T>
_Host_Device complex_t<f_T> operator*(const complex_t<f_T> &z1,const complex_t<f_T> &z2);
// _Host_Device complex_t<double> operator*(const complex_t<double> &z1,const complex_t<double> &z2);

template <isFloating f_T>
_Host_Device complex_t<f_T> operator/(const complex_t<f_T> &z1,const complex_t<f_T> &z2);

template <isFloating f_T>
_Host_Device void operator+=(complex_t<f_T> &z1,const complex_t<f_T> &z2);
// template <isFloating f_T>
// _Host_Device void operator*=(const complex_t<f_T> &z1,const complex_t<f_T> &z2);

// Real & Complex 只有同精度的才可以直接四则运算
template <isFloating f_T>
_Host_Device complex_t<f_T> operator+(const f_T &f1,const complex_t<f_T> &z2);

template <isFloating f_T>
_Host_Device complex_t<f_T> operator+(const complex_t<f_T> &z1,const f_T &f2);

template <isFloating f_T>
_Host_Device complex_t<f_T> operator-(const f_T &f1,const complex_t<f_T> &z2);

template <isFloating f_T>
_Host_Device complex_t<f_T> operator-(const complex_t<f_T> &z1,const f_T &f2);

template <isFloating f_T>
_Host_Device complex_t<f_T> operator*(const f_T &f1,const complex_t<f_T> &z2);

template <isFloating f_T>
_Host_Device complex_t<f_T> operator*(const complex_t<f_T> &z1,const f_T &f2);
// _Host_Device complex_t<double> operator*(const complex_t<double> &z1,const double &f2);

template <isFloating f_T>
_Host_Device complex_t<f_T> operator/(const f_T &f1,const complex_t<f_T> &z2);

template <isFloating f_T>
_Host_Device complex_t<f_T> operator/(const complex_t<f_T> &z1,const f_T &f2);

// template <typename f_T>
// _Host_Device f_T abs(const f_T &f);

template <typename f_T>
_Host_Device f_T min(const f_T &f1, const f_T &f2);

template <typename f_T>
_Host_Device f_T max(const f_T &f1, const f_T &f2);



template<typename f_T>
_Host_Device void ADD(complex_t<f_T> &z, const complex_t<f_T> &x, const complex_t<f_T> &y);
template<typename f_T>
_Host_Device void SUB(complex_t<f_T> &z, const complex_t<f_T> &x, const complex_t<f_T> &y);
template<typename f_T>
_Host_Device void MUL(complex_t<f_T> &z, const complex_t<f_T> &x, const complex_t<f_T> &y);
template<typename f_T>
_Host_Device void DIV(complex_t<f_T> &z, const complex_t<f_T> &x, const complex_t<f_T> &y);


_Host_Device complex_t<double> INV(const complex_t<double>& z);


template<typename f_T>
_Host_Device void SUB(complex_t<f_T> &z, const f_T &x, const complex_t<f_T> &y);

template<typename f_T>
_Host_Device void conj(complex_t<f_T> &z, const complex_t<f_T> &x);

template<typename f_T>
_Host_Device void MUL(complex_t<f_T> &z, const f_T &x, const complex_t<f_T> &y);

_Device void sqrt(complex_t<float> &z, const complex_t<float> &x);
_Device void sqrt(complex_t<double> &z, const complex_t<double> &x);

_Host_Device void SUB(complex_t<double> &z, const complex_t<double> &zD, const complex_t<float> &zF);
_Host_Device void SUB(complex_t<double> &z, const complex_t<float> &zF, const complex_t<double> &zD);
_Host_Device void SUB(complex_t<float> &z, const complex_t<double> &zD1, const complex_t<double> &zD2);

_Host_Device void MUL(complex_t<double> &z, const complex_t<float> &x, const complex_t<double> &y);

template <isFloating f_T>
_Host_Device void operator-=(complex_t<f_T> &z1,const complex_t<f_T> &z2);

template<typename f_T>
_Host_Device f_T wedge_product(const complex_t<f_T> &z1, const complex_t<f_T> &z2);