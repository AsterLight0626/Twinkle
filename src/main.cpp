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

#include <vector>
#include <sstream>

std::vector<double> readNthLine(const std::string& filename, int n) {
    std::ifstream file(filename);
    std::string line;
    int currentLine = 0;

    while (std::getline(file, line)) {
        if (currentLine == n) {
            std::vector<double> numbers;
            std::stringstream ss(line);
            std::string number;

            while (std::getline(ss, number, ',')) {
                numbers.push_back(std::stod(number));
            }

            return numbers;
        }
        currentLine++;
    }

    return {};  // Return empty vector if the line is not found
}


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
    // driver.free();


    // 读入参数用
    std::string filename = "/home/suweiwang/Twinkle/normal_range.txt";
    std::vector<double> numbers;    

    // const int repeat = 1;
    int Ns=100, Nq = 100;
    double temptime;


    // std::string path = "/home/suweiwang/Twinkle/normal_range_res/";
    std::string path = "./";
    std::string name = "time_2D_Txh_1219";
    std::string txt = ".txt";
    name = path + name + txt;

    std::ofstream ofs;     // 创建流对象

    ofs.precision(16);    
    ofs.open(name, std::ios::out);      // 以 write 模式打开 test.txt



    double time_write;
    // double xmin,xmax,ymin,ymax,deltax,deltay;
    // int NaxisX = 30, NaxisY = 30;


    for(int i=0;i<Ns;i++)
    // for(int i=0;i<1;i++)
    {
        for(int j=0;j<Nq;j++)
        // for(int j=25;j<26;j++)
        {
            {

                numbers = readNthLine(filename, i*Nq+j);
                // numbers = readNthLine(filename, 11*100+97);

                // if (numbers.empty()) {
                //     std::cout << "Line " << i*Nq+j << " not found." << std::endl;
                // }
                // else
                {
                    double ss = numbers[0];
                    double qq = numbers[1];
                    double xmax = numbers[2];
                    double xmin = numbers[3];
                    double ymax = numbers[4];
                    double ymin = numbers[5];
                    // double ss = 0.9361015017623882;
                    // double qq = 0.0001;
                    // double xmax = -0.08184367955441649;
                    // double xmin = -0.1796960807612752;
                    // double ymax = 0.02655526837492009;
                    // double ymin = 0.008455501382753586;

                    driver.set_params( ss, qq, 1e-4, xmax, xmin, ymax, ymin );
                    driver.run (  );
                    driver.p_dev->sync_all_streams(  );
                    for(int repeati = 0;repeati<3;repeati++)
                    {   
                        
                        auto start = std::chrono::system_clock::now();
                        driver.run (  );
                        driver.p_dev->sync_all_streams(  );                            
                        auto end = std::chrono::system_clock::now();
                        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

                        temptime = double(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;
                        if(repeati==0)
                        {
                            time_write = temptime;
                        }
                        else
                        {
                            if(temptime < time_write)
                            {
                                time_write = temptime;
                            }
                        }
                    }
                    
                    auto & pool_mag( driver.vp_sol[ 0 ]->pool_mag );
                    pool_mag.cp_d2h( * driver.p_dev );
                    auto & pool_center( driver.vp_sol[ 0 ]->pool_center );
                    pool_center.cp_d2h( * driver.p_dev );
                    auto & pool_extended( driver.vp_sol[ 0 ]->pool_extended );
                    pool_extended.cp_d2h( * driver.p_dev );
                    auto & pool_margin( driver.vp_sol[ 0 ]->pool_margin );
                    pool_margin.cp_d2h( * driver.p_dev );






                    ofs<<ss<<","<<qq<<","<<1e-4<<","<<time_write<<","<<std::endl;
                    // printf("time_write: %.16f\n",time_write);

                    std::ofstream ofs_mag;     // 创建流对象
                    ofs_mag.precision(16);

                    std::string name = "src_ext";
                    std::string txt = ".txt";
                    std::string suffix = "_sq"+std::to_string(i*Nq+j);
                    std::string path = "/home/suweiwang/Twinkle/normal_range_res/";
                    name = path + name + suffix + txt;
                    {
                        ofs_mag.open(name, std::ios::out);      // 以 write 模式打开 test.txt
                        for(int iii=0;iii<driver.vp_sol[ 0 ]->n_src;iii++)
                        {

                            ofs_mag<<pool_center.dat_h[iii].loc_centre.re<<",";
                            ofs_mag<<pool_center.dat_h[iii].loc_centre.im<<",";
                            ofs_mag<<pool_center.dat_h[iii].rho<<",";
                            ofs_mag<<pool_mag.dat_h[iii].mag<<",";
                            ofs_mag<<pool_mag.dat_h[iii].err<<",";
                            ofs_mag<<pool_extended.dat_h[iii].SolveSucceed<<",";
                            ofs_mag<<pool_extended.dat_h[iii].Break<<",";
                            ofs_mag<<std::endl;
                            
                        }
                        ofs_mag.close();
                        // std::cout<<"print file: "<<name<<std::endl;
                    }


                    // int i_src = 313;
                    // int n_point_max = driver.vp_sol[ 0 ]->n_point_max;
                    // {
                    //     std::ofstream ofs_img;     // 创建流对象
                    //     name= "img_pt" + std::to_string(i_src);
                    //     name = path + name + suffix + txt;

                    //     ofs_img.open(name, std::ios::out);      // 以 write 模式打开 test.txt
                    //     ofs_img.precision(16);

                    //     // auto & test = pool_margin.dat_h[ i_src * n_point_max + 2 ].images[1].position.re;
                    //     // printf("test: %.16f\n",test);
                    //     // 第i个src，第j个root，实部，虚部
                    //     for (int j = 0; j < 5; j++) {
                    //         for (int i = 0; i < n_point_max; i++) {
                    //             ofs_img << i << ",";
                    //             ofs_img << j << ",";
                    //             ofs_img<<pool_margin.dat_h[ i_src * n_point_max + i ].images[j].position.re<<",";
                    //             ofs_img<<pool_margin.dat_h[ i_src * n_point_max + i ].images[j].position.im<<",";
                    //             ofs_img<<pool_margin.dat_h[ i_src * n_point_max + i ].images[j].physical << ",";
                    //             ofs_img<<pool_margin.dat_h[ i_src * n_point_max + i ].images[j].parity << "," ;
                    //             ofs_img<<pool_margin.dat_h[ i_src * n_point_max + i ].deltaS[j]<<",";
                    //             ofs_img<<pool_margin.dat_h[ i_src * n_point_max + i ].next_idx[j]<<",";
                    //             ofs_img<<pool_margin.dat_h[ i_src * n_point_max + i ].next_j[j]<<",";
                    //             ofs_img<<pool_margin.dat_h[ i_src * n_point_max + i ].deltaS_Err[j]<<",";
                    //             // ofs_img<<pool_margin.dat_h[ i_src * n_point_max + i ].deltaS_new[j]<<",";
                    //             // ofs_img<<pool_margin.dat_h[ i_src * n_point_max + i ].Err_new[j]<<",";
                    //             ofs_img<<pool_margin.dat_h[ i_src * n_point_max + i ].dz[j].re<<",";
                    //             ofs_img<<pool_margin.dat_h[ i_src * n_point_max + i ].dz[j].im<<",";
                    //             ofs_img<<pool_margin.dat_h[ i_src * n_point_max + i ].wedge[j]<<",";

                    //             ofs_img << std::endl;
                    //         }        
                    //     }
                    //     ofs_img.close();
                    //     std::cout<<"print file: "<<name<<std::endl;
                    // }

                    // {
                    //     std::ofstream ofs_src;     // 创建流对象
                    //     name= "src_pt" + std::to_string(i_src);
                    //     name = path + name + suffix + txt;
                    //     ofs_src.open(name, std::ios::out);      // 以 write 模式打开 test.txt
                    //     ofs_src.precision(16);
                        
                    //     for(int i = 0;i<n_point_max;i++)
                    //     {
                    //         ofs_src<<pool_margin.dat_h[ i_src * n_point_max + i ].position.re<<",";
                    //         ofs_src<<pool_margin.dat_h[ i_src * n_point_max + i ].position.im<<",";
                    //         ofs_src<<pool_margin.dat_h[ i_src * n_point_max + i ].next_src_idx<<",";
                    //         ofs_src<<pool_margin.dat_h[ i_src * n_point_max + i ].Q<<",";
                    //         // ofs_src<<h_srcs[iii].margin_pts[i].quantile;
                    //         ofs_src<<pool_margin.dat_h[ i_src * n_point_max + i ].prev_src_idx<<",";
                    //         ofs_src<<pool_margin.dat_h[ i_src * n_point_max + i ].skip<<",";
                    //         ofs_src<<pool_margin.dat_h[ i_src * n_point_max + i ].Nphys<<",";
                    //         ofs_src<<pool_margin.dat_h[ i_src * n_point_max + i ].error_interval<<",";
                    //         ofs_src<<pool_margin.dat_h[ i_src * n_point_max + i ].area_interval<<",";

                    //         ofs_src<<std::endl;
                    //     }
                    //     ofs_src.close();
                    //     std::cout<<"print file: "<<name<<std::endl;
                    // }



                    printf("qj: %d\n",j);
    
                }
            }
            
        }
        printf("si: %d\n",i);
    }



    ofs.close();
    std::cout<<"print file: "<<name<<std::endl;



    // twinkle::driver_t driver;

    // double ss = 0.9361015017623882;
    // double qq = 0.0001;
    // double xmax = -0.08184367955441649;
    // double xmin = -0.1796960807612752;
    // double ymax = 0.02655526837492009;
    // double ymin = 0.008455501382753586;

    // driver.init( ss, qq, 1e-4, xmax, xmin, ymax, ymin );
    // // driver.run (  );
    // auto start = high_resolution_clock::now();
    // // printf("time start\n");
    // std::cout << "???" << std::endl;
    // for(int i=0;i<repeat;i++)
    // {
    //     driver.run (  );
    // }
    // driver.p_dev->sync_all_streams(  );

    // auto end = high_resolution_clock::now();
    // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);


    // auto & pool_mag( driver.vp_sol[ 0 ]->pool_mag );
    // pool_mag.cp_d2h( * driver.p_dev );

    // double temptime = double(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;
    // // printf("%.1f ms used for %d repeats and %d n_srcs, n_th = %d\n",temptime*1000, repeat, driver.vp_sol[ 0 ]->n_src * 8, driver.vp_sol[ 0 ]->n_th);
    // printf("%.1f ms used for %d repeats and %d n_srcs, n_th = %d\n",temptime*1000, repeat, driver.vp_sol[ 0 ]->n_src, driver.vp_sol[ 0 ]->n_th);

    // auto & pool_center( driver.vp_sol[ 0 ]->pool_center );
    // pool_center.cp_d2h( * driver.p_dev );

    // auto & pool_extended( driver.vp_sol[ 0 ]->pool_extended );
    // pool_extended.cp_d2h( * driver.p_dev );

    // std::ofstream ofs_mag;     // 创建流对象
    // ofs_mag.precision(16);

    // std::string name = "src_ext";
    // std::string txt = ".txt";
    // std::string suffix = "_sq3500_4090_";
    // std::string path = "./";
    // name = path + name + suffix + txt;
    // {
    //     ofs_mag.open(name, std::ios::out);      // 以 write 模式打开 test.txt
    //     for(int iii=0;iii<driver.vp_sol[ 0 ]->n_src;iii++)
    //     {

    //         ofs_mag<<pool_center.dat_h[iii].loc_centre.re<<",";
    //         ofs_mag<<pool_center.dat_h[iii].loc_centre.im<<",";
    //         ofs_mag<<pool_center.dat_h[iii].rho<<",";
    //         ofs_mag<<pool_mag.dat_h[iii].mag<<",";
    //         ofs_mag<<pool_mag.dat_h[iii].err<<",";
    //         ofs_mag<<pool_extended.dat_h[iii].SolveSucceed<<",";
    //         ofs_mag<<pool_extended.dat_h[iii].Break<<",";
    //         ofs_mag<<std::endl;
            
    //     }
    //     ofs_mag.close();
    //     std::cout<<"print file: "<<name<<std::endl;
    // }

    driver.free();

    return 0;
}
