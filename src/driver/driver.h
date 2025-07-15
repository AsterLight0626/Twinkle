#pragma once

#include <memory>
#include <vector>
#include "../device/general.h"
#include "../source/source_base.h"

namespace twinkle
{
////////////////////////////////////////////////////////////
//

class driver_t
{
    ////////// Con-/destructor and initialization //////////
public:
     driver_t(  );
    ~driver_t(  );
    virtual void init( const int n_srcs, const int device_num=0, const int n_stream = 1 );
    // virtual void set_params_2D( double ss, double qq, double rho, double xmax, double xmin, double ymax, double ymin, int Nx, int Ny );
    // virtual void set_params_1D( double ss, double qq, double rho, double xmax, double xmin, double ymax, double ymin, int Nsrc );
    virtual void set_params( double ss, double qq, double rho, double RELTOL, double* xs, double* ys );
    virtual void set_params( const double* ss, double qq, double rho, double RELTOL, double* xs, double* ys );
    virtual void return_mag_to( double* mag );
    virtual void free(  );    

    int n_srcs_all;

    ////////// Device module //////////
public:
    std::shared_ptr<      device_t > p_dev;
    
    ////////// Solver module //////////
public:
    std::vector< std::shared_ptr< source_base_t > > vp_sol;
    
    ////////// Overall interface //////////
public:
    virtual void run(  );
};

};
