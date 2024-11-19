// #include<iostream>

// #include"calculation/init.h"
// #include"MCMC/MCMC.h"
// #include<sstream>
// #include<fstream>

// #define cout std::cout
// #define endl std::endl


// // #include <vector>
// // #include <sstream>

// // nsys profile --stats=true ./bin/Twinkle 

// int main()
// {
//     ////// INPUT OBSERVATIONAL MAGNIFICATION
//     double obsMag[NSRCS], obsErr[NSRCS];
//     std::ifstream inFile("Magnification.txt");

//     if (!inFile.is_open()) {
//     std::cerr << "Error: Input file 'Magnification.txt' does not exist." << endl;
//     return 1;
//     }

//     std::string line;
//     for(int i=0;i<NSRCS;i++)        // NSRCS data in file
//     {
//         std::getline(inFile, line);
//         obsMag[i] = std::stod(line);
//         obsErr[i] = 0.1;
//     }


//     ////// SET MEMORY
//     Twinkle<double> LittleStar;
//     LittleStar.malloc(0);       // set GPU device number


//     ////// INPUT 
//     // input observation time
//     double time[NSRCS];
//     for(int iii=0;iii<NSRCS;iii++)
//     {
//         time[iii] = double(iii) / double(NSRCS) * 1e-1 - 5e-2;
//     }
//     LittleStar.set_time(time);

//     // input parameters
//     int MCMCsteps = 100;
//     src_params_t<double> params_list[MCMCsteps];
//     src_params_t<double>& params0 = params_list[0];
//     // Real Parameters
//     // params0.t_0 = 0;
//     // params0.t_E = 1;
//     // params0.u_0 = 0.035;
//     // params0.alpha = 0;
//     // params0.shape.rho = (1e-4);
//     // params0.q = 1e-5;
//     // params0.s = 1;
//     params0.t_0 = 0 + 0.005;
//     params0.t_E = 1+ 0.005;
//     params0.u_0 = 0.035+ 0.005;
//     params0.alpha = 0+ 0.005;
//     params0.shape.rho = (1e-4);
//     params0.q = 1e-5;
//     params0.s = 1+ 0.005;
//     // LittleStar.set_params(&params);

//     LittleStar.set_path(&params0);


//     ////// SOLVE
//     LittleStar.solve();


//     ////// COPY BACK
//     LittleStar.cp_back_ME();

//     ////// MCMC START
//     double LogProbOld, LogProbNew, alpha, u;
//     std::default_random_engine generator(626212);
//     std::uniform_real_distribution<double> u01(0.,1.);


//     LogProbOld = LogProbGaussian(obsMag,LittleStar.h_Mag.data,obsErr,NSRCS);

//     for(int i=1;i<MCMCsteps;i++)
//     {
//         q_Sampling(params_list[i],params_list[i-1],generator);
//         LittleStar.set_path(&params_list[i]);
//         LittleStar.solve();
//         LittleStar.cp_back_ME();
//         LogProbNew = LogProbGaussian(obsMag,LittleStar.h_Mag.data,obsErr,NSRCS);
//         alpha = Alpha(LogProbOld,LogProbNew);
//         u = u01(generator);
//         // cout<<"i: "<<i<<", u: "<<u<<", alpha: "<<alpha<<endl;
//         // cout<<"LPOld: "<<LogProbOld<<", LPNew: "<<LogProbNew<<endl;
//         if(u>alpha)     // if u > alpha, keep old parameters
//         {
//             params_list[i] = params_list[i-1];
//         }
//         else            // accept new parameters
//         {
//             LogProbOld = LogProbNew;
//         }
//     }

//     std::ofstream ofs;     // 创建流对象
//     ofs.precision(12);    
//     ofs.open("Parameters.txt", std::ios::out);
//     for(int i=0;i<MCMCsteps;i++)
//     {
//         ofs<<params_list[i].t_0<<",";
//         ofs<<params_list[i].t_E<<",";
//         ofs<<params_list[i].u_0<<",";
//         ofs<<params_list[i].alpha<<",";
//         ofs<<params_list[i].shape.rho<<",";
//         ofs<<params_list[i].s<<",";
//         ofs<<params_list[i].q<<endl;
//     }
    



//     // FREE
//     LittleStar.free();
//     return 0;
// }


