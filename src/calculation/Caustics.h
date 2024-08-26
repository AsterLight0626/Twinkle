#pragma once

#include"StructComplex.h"
#include"MacroVaribles.h"
#include"CoeffGenerator.h"
#include"StructSrc.h"
#include"RootsFinder.h"
#include"Operators.h"
#include"../device/device_macros.h"

// 2 body only
template<isFloating f_T>
_Device void CritCoeff_2b(const complex_t<f_T>* params, const f_T& psi, complex_t<f_T>* res, bool print=false);

// template<isFloating f_T>
// _Global void SolveCrit(const src_params_t<f_T>* src_params,complex_t<f_T>* CritLoc);

// template<isFloating f_T>
// _Global void SolveCritCaus(const src_params_t<f_T>* src_params,complex_t<f_T>* CritLoc, complex_t<f_T>* CausLoc);

template<isFloating f_T>
_Device bool CritAddRoot(const complex_t<f_T>& root, const complex_t<f_T>* params);



// // size  <<<(NWALKERS,1,1),(NCRIT)>>>
// template<isFloating f_T>
// _Global void CNext(complex_t<f_T>* caus, int* next_j);

template<isFloating f_T>
_Host_Device f_T TriAreaD(f_T x0, f_T y0, f_T x1, f_T y1, f_T x2, f_T y2);

// template<isFloating f_T>
// _Global void CollideTest(int* srclist, src_ext_t<f_T>* srcs, const complex_t<f_T>* caus, const int* next_j, int batchidx);

// template<isFloating f_T>
// _Global void CollideTest2(int* srclist, src_ext_t<f_T>* srcs, const complex_t<f_T>* caus, const int* next_j, int batchidx);

template <isFloating f_T>
_Global void SolveCritCaus(const src_params_t<f_T>* src_params, CC_t<f_T>* CC, int batchidx);

// template<isFloating f_T>
// _Global void CNext(CC_t<f_T>* CC, int batchidx);

template<isFloating f_T>
_Global void CNext2(CC_t<f_T>* CC, int batchidx);

template<isFloating f_T>
_Global void psiLoc(CC_t<f_T>* CC, int batchidx);

template<isFloating f_T>
_Global void CollideTest(int* srclist, src_ext_t<f_T>* srcs, const CC_t<f_T>* CC, int batchidx);

template<isFloating f_T>
_Global void CollideTest2(int* srclist, src_ext_t<f_T>* srcs, const CC_t<f_T>* CC, int batchidx);
