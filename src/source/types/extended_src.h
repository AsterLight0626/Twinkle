#pragma once

#include "point_img.h"
#include "point_src.h"

namespace twinkle
{

template < class f_T >
struct src_shape_t
{
	f_T                     rho;
	complex_t< f_T > loc_centre;
	f_T                src_area;
};


template< class f_T>
struct src_ext_t
{
	bool                SolveSucceed;
	bool                       Break;    
	// f_T                     src_area;
	// f_T                          Mag;
	// f_T                          Err;
	f_T                         Area;
	f_T                        Err_A;
	// f_T                        Mag_P;
	// f_T                    Err_Mag_P;
	// int                  points_used;	
	// int                       Ncross;
	complex_t< f_T > roots_point[5];
 
};

struct cross_info_t
{
	int          j1;
    int     j_cross;
    int   idx_cross;
    bool additional;
};

template< class f_T >
struct    src_ret_t
{
    f_T mag;
    f_T err;
};


template< class f_T>
struct    shared_info_t
{
	// bool                SolveSucceed;
	bool                       Break;   
	int Ncross;
	// int in_global;
	f_T deltaS_cross_global;
	f_T    Err_cross_global;
};


};                              // namespace twinkle