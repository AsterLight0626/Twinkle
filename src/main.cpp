#include <iostream>
#include <chrono>
#include <string>
#include "driver/driver.h"
#include <fstream>
#include "source/types/point_img.h"
#include "source/types/point_src.h"
#include "source/types/extended_src.h"
// #include "source/types/sh_manager.h"

////////////////////////////////////////////////////////////
// 

int main()
{
    static const int Nsrcs = 1000;

    // coordinates of sources center
    double xs[Nsrcs];
    double ys[Nsrcs];
    for(int idx=0;idx<Nsrcs;idx++)
    {
        xs[idx] = double(idx) / double(Nsrcs) * 1e-4 - 1.5001 + 0.11*1e-3;
        ys[idx] = -0.0035;
    }
    double ss = 0.5;
    double qq = 1e-6;
    double rho = 1e-4;
    double RelTol = 1e-4;

    // set number of streams
    int n_stream = 10;
    twinkle::driver_t driver( n_stream );

    // set device, allocate memory space
    int device_num = 0;
    driver.init( Nsrcs, device_num );

    // set parameters
    driver.set_params(ss,qq,rho,RelTol,xs,ys);

    // solve magnification
    driver.run (  );
    driver.p_dev->sync_all_streams(  );        

    double magnification[Nsrcs];
    driver.return_mag_to(magnification);

    driver.free();

    return 0;
}
