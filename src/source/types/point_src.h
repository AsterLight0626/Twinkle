#pragma once

#include "point_img.h"

namespace twinkle
{
////////////////////////////////////////////////////////////
//

template< class f_T, int n_body = 2 >
struct src_pt_t
{
    ////////// Max number of images //////////
    static constexpr int n_img = n_body * n_body + 1 ;
    
    ////////// Data //////////
    complex_t < f_T >        position;
	complex_t < f_T >     dz[ n_img ];    
	f_T                wedge[ n_img ];
    image_pt_t< f_T > images[ n_img ];

	f_T deltaS              [ n_img ];
	f_T deltaS_Err          [ n_img ];
	// f_T deltaS_new          [ n_img ];
	// f_T Err_new             [ n_img ];	

	int next_idx            [ n_img ];
	int next_j              [ n_img ];
	int                         Nphys;
	f_T                   special_Err;

	int                  next_src_idx;
	int                  prev_src_idx;
	bool                         skip;
	f_T                             Q;

	f_T                error_interval;
	f_T                 area_interval;

};

};                              // namespace twinkle


