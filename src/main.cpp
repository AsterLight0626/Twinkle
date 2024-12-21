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
    using std::chrono::high_resolution_clock;
    

    twinkle::driver_t driver;

    double ss = 0.9361015017623882;
    double qq = 0.0001;
    double xmax = -0.08184367955441649;
    double xmin = -0.1796960807612752;
    double ymax = 0.02655526837492009;
    double ymin = 0.008455501382753586;

    driver.init(  );
    driver.set_params( ss, qq, 1e-4, xmax, xmin, ymax, ymin );
    driver.run (  );
    driver.p_dev->sync_all_streams(  );

    auto & pool_mag( driver.vp_sol[ 0 ]->pool_mag );
    pool_mag.cp_d2h( * driver.p_dev );

    // ofs_mag<<pool_mag.dat_h[iii].mag<<",";
    // ofs_mag<<pool_mag.dat_h[iii].err<<",";

    driver.free();

    return 0;
}
