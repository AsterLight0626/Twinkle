#pragma once

#include"StructComplex.h"
#include"MacroVaribles.h"
#include"CoeffGenerator.h"
#include"StructSrc.h"

static const double ROUND_OFF_D=2e-15;
static const float ROUND_OFF_F=1e-6;

static const int MAXIT = 80;
static const int FRAC_JUMP_EVERY = 16;      // must be 2^n
static const double dlmin = 1e-4;

template<isFloating f_T>
_Device bool cmplx_roots_gen(complex_t<f_T> *roots, complex_t<f_T> *poly, int degree, bool polish_roots_after, bool use_roots_as_starting_points);
template<isFloating f_T>
_Device bool cmplx_roots_gen(image_pt_t<f_T> *imgs, complex_t<f_T> *poly, int degree, bool polish_roots_after, bool use_roots_as_starting_points);

template<isFloating f_T>
_Device void solve_quadratic_eq(complex_t<f_T> &x0, complex_t<f_T> &x1, complex_t<f_T> *poly);

template<isFloating f_T>
_Device void cmplx_newton_spec(complex_t<f_T> *poly, int degree, complex_t<f_T> *root, int &iter, bool &success);

template<isFloating f_T>
_Device void cmplx_laguerre(complex_t<f_T> *poly, int degree, complex_t<f_T> *root, int &iter, bool &success);

template<isFloating f_T>
_Device void cmplx_laguerre2newton(complex_t<f_T> *poly, int degree, complex_t<f_T> *root, int &iter, bool &success, int starting_mode);

