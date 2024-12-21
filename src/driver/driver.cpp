#include "driver.h"

namespace twinkle
{
////////////////////////////////////////////////////////////
//

driver_t::driver_t(  )
{
    #warning "TEST"
    vp_sol.resize( 1 );
    return;
}

driver_t::~driver_t(  )
{
    p_dev->finalize(  );
    return;
}

void driver_t::init(  )
{
    p_dev = std::make_shared< device_t > (  );
    p_dev-> init   ( 0 );
    p_dev-> prepare(   ); 
    for( auto & p_sol : vp_sol )
    {
        p_sol = std::make_shared< source_base_t > (  );
        p_sol->init( * p_dev );
    }
    return;
}

void driver_t::set_params( double ss, double qq, double rho, double xmax, double xmin, double ymax, double ymin )
{
    for( auto & p_sol : vp_sol )
    {
        p_sol->set_params( * p_dev, ss, qq, rho, xmax, xmin, ymax, ymin );
    }
    return;
}

void driver_t::free(  )
{
    for( auto & p_sol : vp_sol )    
        p_sol->free( * p_dev );
    return;
}

////////////////////////////////////////////////////////////
//

void driver_t::run(  )
{
    for( auto & p_sol : vp_sol )    
        p_sol-> run ( *  p_dev );
    // p_dev->sync_all_streams(   );
    return;
}

};
 
