#pragma once

namespace twinkle
{
////////////////////////////////////////////////////////////
//

template< class f_T >
struct image_pt_t
{
    ////////// 
    complex_t< f_T > position;
    bool             physical;
    bool               parity;

	// __device__ image_pt_t(  ) {  };
	// __device__ image_pt_t( const image_pt_t< f_T > & s )
    // {
    //     this->position = s.position;
    //     this->physical = s.physical;
    //     this->parity   = s.  parity;
    //     return;
    // };
};

};                              // namespace twinkle
