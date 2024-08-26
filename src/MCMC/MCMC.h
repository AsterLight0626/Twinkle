#pragma once

#include<random>
#include<cmath>
#include"../calculation/init.h"


double LogProbGaussian(double* obsMag, double* mockMag, double* sigma, int len);

void q_Sampling(src_params_t<double>& params_new, src_params_t<double>& params_old, std::default_random_engine generator);

double Alpha(double log_prob_old, double log_prob_new);