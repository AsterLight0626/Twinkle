#include "driver.h"

namespace twinkle
{
////////////////////////////////////////////////////////////
//


driver_t::driver_t( const int n_stream )
{
    // #warning "TEST"
    vp_sol.resize( n_stream );
    return;
}

driver_t::~driver_t(  )
{
    p_dev->finalize(  );
    return;
}

void driver_t::init( const int n_srcs, const int device_num )
{
    n_srcs_all = n_srcs;
    p_dev = std::make_shared< device_t > (  );
    p_dev-> init   ( device_num );
    p_dev-> prepare(   ); 

    int n_src_per_stream = (n_srcs-1) / vp_sol.size() + 1;  // vp_sol.size = nstream
    // for( auto & p_sol : vp_sol )
    for(int i=0;i<vp_sol.size()-1;i++)
    {
        auto & p_sol = vp_sol[i];
        p_sol = std::make_shared< source_base_t > (  );
        p_sol->n_src = n_src_per_stream;
        p_sol->init( * p_dev );
    }
    auto & p_sol = vp_sol[vp_sol.size()-1];
    p_sol = std::make_shared< source_base_t > (  );
    p_sol->n_src = n_srcs - n_src_per_stream * ( vp_sol.size()-1 );
    p_sol->init( * p_dev );

    return;
}

void driver_t::set_params_2D( double ss, double qq, double rho, double xmax, double xmin, double ymax, double ymin, int Nx, int Ny )
{
    for( auto & p_sol : vp_sol )
    {
        p_sol->set_params_2D( * p_dev, ss, qq, rho, xmax, xmin, ymax, ymin, Nx, Ny );
    }
    return;
}

void driver_t::set_params_1D( double ss, double qq, double rho, double xmax, double xmin, double ymax, double ymin, int Nsrc )
{
    for( auto & p_sol : vp_sol )
    {
        p_sol->set_params_1D( * p_dev, ss, qq, rho, xmax, xmin, ymax, ymin, Nsrc );
    }
    return;
}

void driver_t::set_params( double ss, double qq, double rho, double RELTOL, double* xs, double* ys )
{
    for( auto & p_sol : vp_sol )
    {
        p_sol->s = ss;
        p_sol->q = qq;
        p_sol->RelTol    =   RELTOL;
        p_sol->m1 =  1 / ( qq + 1 );
        p_sol->m2 = qq / ( qq + 1 );
    }
    int n_src_per_stream = (n_srcs_all-1) / vp_sol.size() + 1;  // vp_sol.size = nstream
    for(int idx=0;idx<n_srcs_all;idx++)
    {
        auto & shape  = vp_sol[ idx/n_src_per_stream ]->pool_center.dat_h[ idx%n_src_per_stream ];
        shape.rho = rho;
        shape.loc_centre.re = xs[idx];
        shape.loc_centre.im = ys[idx];            
    }
    for( auto & p_sol : vp_sol )
    {
        p_sol->pool_center.cp_h2d( * p_dev );
    }

    return;
}

void driver_t::return_mag_to( double* mag )
{
    int n_src_per_stream = (n_srcs_all-1) / vp_sol.size() + 1;  // vp_sol.size = nstream
    for( auto & p_sol : vp_sol )
    {
        p_sol->pool_mag.cp_d2h( * p_dev );
    }
    for(int idx=0;idx<n_srcs_all;idx++)
    {
        mag[idx] = vp_sol[ idx/n_src_per_stream ]->pool_mag.dat_h[ idx%n_src_per_stream ].mag;
    } 
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
 
