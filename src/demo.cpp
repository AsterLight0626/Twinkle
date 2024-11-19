#include<iostream>

#include"calculation/init.h"

#define cout std::cout
#define endl std::endl


// #include <vector>
#include <fstream>

// // nsys profile --stats=true ./bin/Twinkle 

int main()
{
    ////// SET MEMORY
    Twinkle<double> LittleStar;
    LittleStar.malloc(0);       // set GPU device number
    // LittleStar.critcaus = true;


    ////// INPUT 
    // input observation time
    double time[NSRCS];
    for(int iii=0;iii<NSRCS;iii++)
    {
        // time[iii] = double(iii) / double(NSRCS) * 1e-1 - 5e-2;
        time[iii] = double(iii) / double(NSRCS) * 1e-4 - 1.5001 + 0.11*1e-3;
    }
    // LittleStar.set_time(time);

    // input parameters
    src_params_t<double> params;
    params.t_0 = 0;
    params.t_E = 1;
    params.u_0 = -0.0035;
    params.alpha = 0;
    params.shape.rho = (1e-4);
    params.q = 1e-6;
    params.s = 0.5;


    // LittleStar.set_params(&params);

    LittleStar.set_path(time,&params);


    ////// SOLVE
    LittleStar.solve();


    ////// COPY BACK MAGNIFICATION AND ERROR
    LittleStar.cp_back_ME();   


    // data in LittleStar.h_Mag (Magnification) and LittleStar.h_Err (Error)
    // .writeto(path) method is recommended
        // std::ofstream ofs;     // 创建流对象
        // ofs.precision(16);    
        // ofs.open("Magnification.txt", std::ios::out);
        // for(int i=0;i<LittleStar.Nsrcs;i++)
        // {
        //     ofs<<LittleStar.h_Mag[i]<<endl;
        // }


    ////// COPY BACK ALL DATA
    LittleStar.cp_back();
    LittleStar.writeto("./");
    

    // FREE
    LittleStar.free();
    return 0;
}


