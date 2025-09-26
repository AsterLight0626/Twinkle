#include <iostream>
#include <chrono>
#include <string>
#include "driver/driver.h"
#include <fstream>

////////////////////////////////////////////////////////////
// 

int main()
{
    static const int Nsrcs = 1024;

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


    int n_stream = 2;
    // twinkle::driver_t driver( n_stream );
    twinkle::driver_t driver;

    int device_num = 2;
    driver.init( Nsrcs, device_num, n_stream, RelTol );
    // driver.set_params(ss,qq,rho,RelTol,xs,ys);
    driver.set_params(ss,qq,rho,xs,ys);

    driver.runLD ( 1.0 );
    // driver.p_dev->sync_all_streams(  );        

    double magnification[Nsrcs];
    driver.return_mag_to(magnification);

    // printf("mag[2]: %.16f\n",magnification[2]);
    // printf("mag[12]: %.16f\n",magnification[12]);

    std::ofstream outfile("magnitude_data.txt");
    if (outfile.is_open()) {
        outfile.precision(16); // 设置16位小数精度
        outfile << std::fixed; // 固定小数点格式
        for(int i_src=0;i_src<Nsrcs;i_src++)
        {
            outfile << magnification[i_src] << std::endl;
        }
        
        outfile.close();
    }

    driver.free();

    return 0;
}
