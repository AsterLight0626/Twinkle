// #include <iostream>
// #include <chrono>
// #include <string>
// #include "driver/driver.h"
// #include <fstream>

// ////////////////////////////////////////////////////////////
// // 

// int main()
// {
//     // static const int Nsrcs = 1;

//     // double xs[Nsrcs];
//     // double ys[Nsrcs];
//     // for(int idx=0;idx<Nsrcs;idx++)
//     // {
//     //     xs[idx] = 0.02;
//     //     ys[idx] = 0.01;
//     // }
//     // double ss = 1;
//     // double qq = 1e-4;
//     // double rho = 1e-2;
//     // double RelTol = 1e-6;

//     static const int Nsrcs = 1024;

//     double xs[Nsrcs];
//     double ys[Nsrcs];
//     for(int idx=0;idx<Nsrcs;idx++)
//     {
//         xs[idx] = double(idx) / double(Nsrcs) * 1e-4 - 1.5001 + 0.11*1e-3;
//         ys[idx] = -0.0035;
//     }
//     double ss = 0.5;
//     double qq = 1e-6;
//     double rho = 1e-4;
//     double RelTol = 1e-6;


//     int n_stream = 1;
//     // twinkle::driver_t driver( n_stream );
//     twinkle::driver_t driver;

//     int device_num = 6;
//     bool astrom = true;
//     driver.init( Nsrcs, device_num, n_stream, RelTol, astrom );
//     // driver.set_params(ss,qq,rho,RelTol,xs,ys);
//     driver.set_params(ss,qq,rho,xs,ys);

//     // driver.run_pt(  );
//     // driver.run(  );
//     double LD_a = 1.0;          // linear limb darkening coefficient
//     driver.runLD ( LD_a );
//     // driver.p_dev->sync_all_streams(  );        

//     double magnification[Nsrcs];
//     driver.return_mag_to(magnification);

//     twinkle::complex_t<double> astrom_Th[Nsrcs];
//     driver.return_astrom_to(astrom_Th);

//     // std::ofstream outfile("magnitude_data.txt");
//     // if (outfile.is_open()) {
//     //     outfile.precision(16);
//     //     outfile << std::fixed;
//     //     for(int i_src=0;i_src<Nsrcs;i_src++)
//     //     {
//     //         outfile << magnification[i_src] << std::endl;
//     //     }
        
//     //     outfile.close();
//     // }
//     printf("Mag: %.16f\n",magnification[0]);
//     for(int i=0;i<64;i++)
//         printf("astrom_Th: %.16f, %.16f\n",astrom_Th[i].re,astrom_Th[i].im);

//     driver.free();

//     return 0;
// }
