#pragma once

#include"StructComplex.h"
#include"MacroVaribles.h"

template<isFloating f_T>
_Host_Device void ParamsReal2Complex(const f_T* params_lens, complex_t<f_T>* params);


// 多项式相乘
// 如要数乘，对一个数取址& 传入，degree = 0
// 仅限数乘时可以把计算结果存回原地
template<typename f_T>
_Host_Device void PolyMultiply(const f_T* coeff1, const int degree1, const f_T* coeff2, const int degree2, f_T* coeff_res);

// 数乘，把数放在前一个输入
template<typename f_T>
_Host_Device void PolyMultiply(const f_T coeff0, const int degree0, const f_T* coeff2, const int degree2, f_T* coeff_res);

// 多项式相加
// 注意：没有初始化
// 可以存在原地
template<typename f_T>
_Host_Device void PolyAdd(const f_T* coeff1, const int degree1, const f_T* coeff2, const int degree2, f_T* coeff_res);


// GPU版和CPU版稍有不同，将loc_src从params中单独拎了出来
template<isFloating f_T>
_Host_Device void PolyCoeff(const complex_t<f_T>* params,const complex_t<f_T>& loc_src, complex_t<f_T>* coeff_res);

// f(z_bar)，注意输入的是 z 不是 zbar
template<isFloating f_T>
_Host_Device complex_t<f_T> f_zbar(const complex_t<f_T>& z, const complex_t<f_T>* params);
// mass centre frame
template<isFloating f_T>
_Host_Device complex_t<f_T> f_zbar(const complex_t<f_T>& z, const f_T s, const f_T m1, const f_T m2);

// d(f(z_bar)) / d(z_bar)，注意对谁求导
template<isFloating f_T>
_Host_Device complex_t<f_T> D_f_zbar(const complex_t<f_T>& z, const complex_t<f_T>* params);

// Jacobi 放大率的倒数
template <isFloating f_T>
_Host_Device f_T Jacobi(const complex_t<f_T>& z, const complex_t<f_T>* params);

template <typename f_T>
_Host_Device f_T PolyValue(const f_T* coeff, const f_T& z, int degree);

template <typename f_T>
_Host_Device void PolyDer(const f_T* coeff,f_T* res,const int prev_degree);

template <isFloating f_T>
_Host_Device f_T PolyCoeffGC(const f_T s, const f_T q, const complex_t<f_T>& zeta_M, complex_t<f_T>* coeffs);

_Host_Device float PolyCoeffGC(const float s, const float q, const complex_t<float>& zeta_M, complex_t<float>* coeffs);

template <isFloating f_T>
_Host_Device f_T PolyCoeffGCOld(const f_T s, const f_T q, const complex_t<f_T>& zeta_M, complex_t<f_T>* coeffs);

template <isFloating f_T>
_Host_Device f_T PolyCoeffVB(const f_T s, const f_T q, const complex_t<f_T>& zeta_M, complex_t<f_T>* coeffs_out);

template <isFloating f_T>
_Host_Device f_T PolyCoeffVB(const f_T s, const f_T m1, const f_T m2, const complex_t<f_T>& zeta_M, complex_t<f_T>* coeffs_out);

template <isFloating f_T>
_Host_Device f_T PolyCoeffVB1(const f_T s, const f_T q, const complex_t<f_T>& zeta_M, complex_t<f_T>* coeffs_out);

template <isFloating f_T>
_Host_Device f_T PolyCoeffVB1(const f_T s, const f_T m1, const f_T m2, const complex_t<f_T>& zeta_M, complex_t<f_T>* coeffs_out);

template <isFloating f_T>
_Host_Device f_T PolyCoeffTwinkle(const f_T s, const f_T m1, const f_T m2, const complex_t<f_T>& zeta_M, complex_t<f_T>* coeffs_out);