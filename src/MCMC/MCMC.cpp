#include "MCMC.h"
#include "../calculation/StructSrc.h"

double LogProbGaussian(double* obsMag, double* mockMag, double* sigma, int len)
{
    double chi2 = 0.;
    for(int i=0;i<len;i++)
    {
        chi2 += ((obsMag[i]-mockMag[i]) * (obsMag[i]-mockMag[i])) / (sigma[i]*sigma[i]);
    }
    return -0.5*chi2;
}


void q_Sampling(src_params_t<double>& params_new, src_params_t<double>& params_old, std::default_random_engine generator)
{
    std::normal_distribution<double> normal01(0.0,1.0);
    double rand01;
    double scale = 0.001;

    rand01 = normal01(generator);
    params_new.t_0 = params_old.t_0 + rand01*scale;

    rand01 = normal01(generator);
    params_new.t_E = params_old.t_E + rand01*scale;
    if(params_new.t_E<0){params_new.t_E = -params_new.t_E;}                  // t_E > 0
    if(params_new.t_E==0){params_new.t_E = params_old.t_E;}


    rand01 = normal01(generator);
    params_new.u_0 = params_old.u_0 + rand01*scale;
    if(params_new.u_0<0){params_new.u_0 = -params_new.u_0;}                  // t_E > 0
    if(params_new.u_0==0){params_new.u_0 = params_old.u_0;}

    rand01 = normal01(generator);
    params_new.s = params_old.s + rand01*scale;
    if(params_new.s<0){params_new.s = -params_new.s;}                        // s > 0
    if(params_new.s==0){params_new.s = params_old.s;}

    rand01 = normal01(generator);
    params_new.q = params_old.q + rand01*scale;
    if(params_new.q<0){params_new.q = -params_new.q;}                        // q > 0
    if(params_new.q==0){params_new.q = params_old.q;}

    rand01 = normal01(generator);
    params_new.alpha = params_old.alpha + rand01*scale;
    params_new.alpha -= floor(params_new.alpha / (2.*PI)) * 2*PI;     // 0 < alpha < 2*PI

    rand01 = normal01(generator);
    params_new.shape.rho = params_old.shape.rho + rand01*scale;
    if(params_new.shape.rho<0){params_new.shape.rho = -params_new.shape.rho;}// shape.rho > 0
    if(params_new.shape.rho==0){params_new.shape.rho = params_old.shape.rho;}


}

double Alpha(double log_prob_old, double log_prob_new)
{
    double alpha = log_prob_new - log_prob_old;
    alpha = exp(alpha);
    alpha = min(1.0,alpha);
    return alpha;
}


