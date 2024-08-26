#pragma once

#include "MacroVaribles.h"
#include "StructComplex.h"
#include "StructSrc.h"
#include "Operators.h"
#include "CoeffGenerator.h"
#include "RootsFinder.h"
#include "../device/device_macros.h"

static const double RelErrAllowed = 1e-10;

template<isFloating f_T>
_Global void MarginBatch(int const * srclist, src_ext_t<f_T>* srcs);



template<isFloating f_T>
_Global void SolveBatch(int const * srclist,src_ext_t<f_T>*srcs,const src_params_t<f_T>* src_params, int batchidx, bool muted = false, bool check=true);

// size: <<<(NWALKERS*NSRCS,1, 1),BATCH_SIZE>>>
template<isFloating f_T>
_Global void FindNextBatch3(int* srclist, src_ext_t<f_T>* srcs, int batch_idx, bool muted = false);

template <isFloating f_T>
_Global void AreaErrBatchSrc4(int* srclist, src_ext_t<f_T>* srcs, int batch_idx, bool muted);

// size <<<[(NSRCS-1) // 64 +1,1,1],64>>>
template<isFloating f_T>
_Global void AreaErrBatchSrcCross(int* srclist, src_ext_t<f_T>* srcs,const src_params_t<f_T>* src_params, int batch_idx, int prev_Ncal, bool muted = false);

template <isFloating f_T>
_Global void SumArea0(int* srclist, src_ext_t<f_T>* srcs, int batchidx, bool muted = false);

template <isFloating f_T>
_Global void SumArea3(int* srclist, src_ext_t<f_T>* srcs, int batchidx, bool muted = false);


static const float C_Q = 6;
static const float C_G = 2;
static const float C_P = 2;

template<isFloating f_T>
_Global void QuadruTst(src_ext_t<f_T>* srcs, const src_params_t<f_T>* src_params, const int Nsrcs=NSRCS);

template <isFloating f_T>
_Global void SetPath(src_ext_t<f_T>*srcs, const f_T* time, const src_params_t<f_T>* move_params, int Nsrcs=NSRCS);

// size:<<<NaxisY,NaxisX>>>; 
template <isFloating f_T>
_Global void SetArray(src_ext_t<f_T>*srcs, const src_params_t<f_T>* move_params, const int NaxisX, const int NaxisY, const f_T xmin,const f_T xmax,const f_T ymin,const f_T ymax, int Nsrcs=NSRCS);

// size:<<<NaxisY,NaxisX>>>; 
template <isFloating f_T>
_Global void SetRhoLoc(src_ext_t<f_T>*srcs, const f_T* rhos, const complex_t<f_T>* zetas, int Nsrcs=NSRCS);

// size:<<<1,1>>>; 
template <isFloating f_T>
_Global void SetPoint(src_ext_t<f_T>*srcs, const src_params_t<f_T>* move_params, const f_T x, const f_T y);

template <isFloating f_T>
_Global void AdaptiveLocLoop(int* srclist, src_ext_t<f_T>*srcs, int batchidx, bool muted = false);

// size <<<[(NWALKERS*NSRCS-1) // 64 +1,1,1],64>>>
template<isFloating f_T>
_Global void PhysicalTest(int* srclist, src_ext_t<f_T>*srcs, int batch_idx, int prev_Ncal, bool muted = false);


// size: <<<(Ncal-1) / 64+1 , min(64,Ncal)>>>
template<isFloating f_T>
_Global void SuccessCheck(int* Ncal, int* srclist, src_ext_t<f_T>* srcs, int prev_Ncal, bool muted = false);
