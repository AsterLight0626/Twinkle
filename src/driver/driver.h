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
    virtual void init(  );
    virtual void set_params( double ss, double qq, double rho, double xmax, double xmin, double ymax, double ymin );
    virtual void free(  );    

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
