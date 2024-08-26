// #include<iostream>

// #include"calculation/init.h"

// #define cout std::cout
// #define endl std::endl


// // #include <vector>
// #include <fstream>

// // // nsys profile --stats=true ./bin/Twinkle 

// int main()
// {
//     ////// SET MEMORY
//     Twinkle<double> LittleStar;
//     LittleStar.critcaus = true;     // must be ahead of malloc


//     LittleStar.malloc(0);       // set GPU device number

//     ////// INPUT 
//     // input parameters
//     src_params_t<double> params;
//     params.t_0 = 0;
//     params.t_E = 0.218999842;
//     params.u_0 = 0.01;
//     params.alpha = 0.53;
//     params.shape.rho = 1.939366667;
//     params.q = 0.125;
//     params.s = 0.437613333;

//     LittleStar.set_params(&params);

//     LittleStar.solve_caus();
//     LittleStar.cp_back_caus();
//     LittleStar.writeto_caus("./");


//     // FREE
//     LittleStar.free();
//     return 0;
// }


