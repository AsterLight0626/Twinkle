#pragma once

#include"StructComplex.h"
#include"StructSrc.h"
#include"MacroVaribles.h"
#include"CoeffGenerator.h"
#include"../utilities/functions.h"
#include"RootsFinder.h"

static const double THRESHOLD_REAL_DOUBLE=1e-5;
static const float THRESHOLD_REAL_FLOAT=3e-3;


// // 使用方法：image_pt_t.Physical = isPhysical(image_pt_t.position,params,src_pt_t.position)
_Device bool isPhysical(const complex_t<double>& lambda, const complex_t<double>* params, const complex_t<double>& loc_src, double& physical_err);
_Device bool isPhysical(const complex_t<float>& lambda, const complex_t<float>* params, const complex_t<float>& loc_src, float& physical_err);

template<isFloating f_T>
_Device bool isPhysicalAll(image_pt_t<f_T>* imgs, const f_T s, const f_T m1, const f_T m2, const complex_t<f_T>& loc_src, int* Nphys, f_T* errors_out);

template <isFloating f_T>
_Device f_T DeltaS1(const complex_t<f_T>& z_here,const complex_t<f_T>& z_next);

template<isFloating f_T>
_Device void BubbleArgSort(int* index, f_T* arr, const int len);

template<isFloating f_T>
_Device int BisearchLeft(const f_T* values_in_order, const f_T to_insert, const int len=BATCH_SIZE);


_Device complex_t<float> Margin_Q(const src_shape_t<float>& shape,const float& q);
_Device complex_t<double> Margin_Q(const src_shape_t<double>& shape,const double& q);

template<isFloating f_T>
_Device f_T Area_src(const src_shape_t<f_T>& shape);


_Device void PathLoc(const src_params_t<float>& move_params, const float time, complex_t<float>& loc);
_Device void PathLoc(const src_params_t<double>& move_params, const double time, complex_t<double>& loc);

template<isFloating f_T>
_Device bool PolishPP(image_pt_t<f_T>* images, const f_T* jacobi,  bool muted=true);

template<isFloating f_T>
_Device f_T Mu_QC(const complex_t<f_T>* params, const complex_t<f_T>& z);

template<isFloating f_T>
_Device f_T GhostTst(const complex_t<f_T>* params, const complex_t<f_T>& z_G, const complex_t<f_T>& z_hat);

template<isFloating f_T>
_Device f_T WedgeDxDDx(const complex_t<f_T>& z, const complex_t<f_T>& zeta, const complex_t<f_T>* params, const f_T& rho, const complex_t<f_T>& loc_centre, complex_t<f_T>& dz_output);

_Device float sqrt_f_T(const float& f);
_Device double sqrt_f_T(const double& f);


template<isFloating f_T>
_Device void ConnectOrder(complex_t<f_T>* posi_here, complex_t<f_T>* posi_other, int len_here, int len_other, int* order, bool print=false);

template<isFloating f_T>
_Device int PrevSrcIdx(const src_pt_t<f_T>* margin_pts, const int& idx_here);
template<isFloating f_T>
_Device int NextSrcIdx(const src_pt_t<f_T>* margin_pts, const int& idx_here);

template<isFloating f_T>
_Host_Device f_T Jacobi_M(const complex_t<f_T>& z, const f_T s, const f_T q);

template<isFloating f_T>
_Device bool PolishPP2(image_pt_t<f_T>* imgs, f_T* jacobis, const complex_t<f_T>& zeta,int idx=-1);