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
    LittleStar.malloc(0);       // set GPU device number 0
    // LittleStar.critcaus = true;


    ////// INPUT 
    // input observation time
    double time[NSRCS];
    for(int iii=0;iii<NSRCS;iii++)
    {
        time[iii] = double(iii) / double(NSRCS) * 1e-1 - 5e-2;
    }
    // LittleStar.set_time(time);

    // input parameters
    src_params_t<double> params;
    params.t_0 = 0;
    params.t_E = 1;
    params.u_0 = 0.035;
    params.alpha = 0;
    params.shape.rho = (1e-4);
    params.q = 1e-5;
    params.s = 1;


    // LittleStar.set_params(&params);

    LittleStar.set_path(time,&params);


    ////// SOLVE
    LittleStar.solve();


    ////// COPY BACK MAGNIFICATION AND ERROR
    LittleStar.cp_back_ME();   



    ////// COPY BACK ALL DATA
    LittleStar.cp_back();
    LittleStar.writeto("./");
    

    // FREE
    LittleStar.free();
    return 0;
}


