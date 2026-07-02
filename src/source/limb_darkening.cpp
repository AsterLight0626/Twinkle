#include "source_base.h"
#include "../utils/solve.h"

namespace twinkle
{

////////////////////////////////////////////////////////////
//
// Host-side Limb Darkening Integration
//

__host__ bool source_base_t::monotest_single(double val_L, double val_R, double Mag_L, double Mag_R, double phi_self, double& val_self, double& Mag_self)
{
    bool changed = false;
    if (val_self > val_L*10. || val_self < val_R*0.1)         // 如果不满足单调性，给这个点插值
    {
        Mag_self = ( Mag_L + Mag_R )/2.;
        val_self = Mag_self * (1 - phi_self*phi_self);
        changed = true;
    }

        
    return changed;
}

__host__ void source_base_t::monotest(twinkle::ErrorUnit_t<double>& eu)
{
    bool changed_L = monotest_single(eu.val_L(),eu.val_C(),eu.raw_Mag[0],eu.raw_Mag[2], eu.phi[1], eu.val_LC(),eu.raw_Mag[1]);
    bool changed_R = monotest_single(eu.val_C(),eu.val_R(),eu.raw_Mag[2],eu.raw_Mag[4], eu.phi[3], eu.val_RC(),eu.raw_Mag[3]);

    if(astrom)
    {
        if(changed_L)
        {
            eu.astrom_X0[1] = (eu.astrom_X0[0] + eu.astrom_X0[2]) / 2.;
        }
        if(changed_R)
        {
            eu.astrom_X0[3] = (eu.astrom_X0[2] + eu.astrom_X0[4]) / 2.;
        }
    }
    return;
}

__host__ void source_base_t::caustic_cal(double lens_s, double lens_q, c_t* caustic_points)
{

    c_t coef[5];
    image_pt_t<f_t> temp_images[4];

    f_t denom = f_t(1) / (f_t(1) + lens_q);
    f_t m0 = denom;               // 1/(1+q)
    f_t m1 = lens_q * denom;            // q/(1+q)

    f_t A = -lens_s * m1;               // z0
    f_t B =  lens_s * m0;                // z1

    f_t A2 = A * A;
    f_t B2 = B * B;
    f_t AB = A * B;
    f_t A_plus_B = A + B;
    f_t A2B2 = A2 * B2;

    f_t cos_psi = 0.6;
    f_t sin_psi = 0.8;

    coef[4].re = f_t(1);
    coef[4].im = f_t(0);

    coef[3].re = f_t(-2) * A_plus_B;
    coef[3].im = f_t(0);

    f_t real_part_2 = A2 + B2 + f_t(4) * AB;
    f_t sum_m = m0 + m1;
    coef[2].re = real_part_2 - sum_m * cos_psi;
    coef[2].im = sum_m * sin_psi;

    f_t const_term_1 = f_t(-2) * AB * A_plus_B;
    f_t coeff_E_1 = f_t(2) * (A * m1 + B * m0);
    coef[1].re = const_term_1 + coeff_E_1 * cos_psi;
    coef[1].im = - coeff_E_1 * sin_psi;

    f_t const_term_0 = A2B2;
    f_t coeff_E_0 = -(m1 * A2 + m0 * B2);
    coef[0].re = const_term_0 + coeff_E_0 * cos_psi;
    coef[0].im = - coeff_E_0 * sin_psi;


    bool fail = root_finder::cmplx_roots_gen(temp_images, coef, 4,true, false);
    for(int caus_i=0;caus_i<4;caus_i++)
    {
        caustic_points[caus_i] = temp_images[caus_i].position + f_zbar(temp_images[caus_i].position, lens_s );
    }
    if(fail)
    {
        caustic_points = nullptr;        
    }
        
    return;
}

__host__ void source_base_t::hidden_phi(c_t* caustic_pts, const c_t& loc_center, double rho, double* phi_hidden)
{
    if(caustic_pts == nullptr)
    {
        phi_hidden[4] = 0;
    }
    else
    {
        int i_phi=0;
        double rho2 = rho*rho;
        for(int test_i=0;test_i<4;test_i++)
        {
            f_t distance2 = (caustic_pts[test_i] - loc_center).norm2( );
            if(distance2 <= rho2)
            {
                phi_hidden[i_phi] = sqrt(1-distance2/rho2);
                i_phi += 1;
            }
        }
        phi_hidden[4] = i_phi;
    }
    return;
}

__host__ void source_base_t::runLD_beta( device_t & f_dev, double LD_a, int* Nuniform_out, int max_depth )
{
    auto n_bl = ( n_src + n_th - 1 ) / n_th;
    auto s_sh = 0;
    auto lpar = std::make_tuple
        ( dim3( n_bl ), dim3( n_th ), s_sh );
    f_dev.launch( point_approximation, lpar, stream, ( * this ) );
    f_dev.sync_all_streams();
    pool_mag.cp_d2h(f_dev);
    if(astrom)
    {
        pool_astrom_Th.cp_d2h(f_dev);
    }
    n_bl = n_src;
    s_sh = sizeof(int)*n_th + 
    std::max( sizeof(float2_t)*n_point_max, 
              sizeof(float2_t)   * 8*n_th
            + sizeof(shared_info_t< float2_t >) * 1 
            + sizeof(cross_info_t) * n_cross_max * 4
            + sizeof(src_pt_t< float2_t >) * n_th);

    double RelTolS;
    if (LD_a != 0)
    {
        RelTolS = min(1., RelTol / LD_a);
    }
    else{ RelTolS = 1.; }
    double coeff_RelTol = 1./3.;
    int NuniformMax = 512;
    int NuniformMin = 2;
    

    float2_t* LDInte = new float2_t[n_src];
    float2_t* M0 = new float2_t[n_src];
    float2_t* Error_now = new float2_t[n_src];
    float2_t* lowerbound = new float2_t[n_src];
    bool* finished = new bool[n_src];
    bool* ini_finished = new bool[n_src];
    for(int i_src=0;i_src<n_src;i_src++)
    {
        finished[i_src] = false;
        ini_finished[i_src] = false;
    }


    int* Nuniform = new int[n_src];


    complex_t<float2_t>* X_LD = new complex_t<float2_t>[n_src];
    complex_t<float2_t>* X0_max_radius = new complex_t<float2_t>[n_src];



    using ErrorUnit = twinkle::ErrorUnit_t<double>;
    std::priority_queue<ErrorUnit>* heaps = new std::priority_queue<ErrorUnit>[n_src];
    ErrorUnit* EU_buffer = new ErrorUnit[n_src];
    ErrorUnit* EU_buffer_R = new ErrorUnit[n_src];
    ////////////////////////////////////////////////////////////////

    for(int i_src=0; i_src<n_src; i_src++)
    {
        EU_buffer[i_src] = ErrorUnit(0.0,1.0);
        EU_buffer[i_src].depth = 2;
        EU_buffer[i_src].val_R() = 0.0;
        EU_buffer[i_src].Ncross[4] = 0;
        EU_buffer[i_src].raw_Mag[4]=pool_mag.dat_h[i_src].mag;
        EU_buffer[i_src].astrom = astrom;
        if(astrom)
        {
            EU_buffer[i_src].astrom_X0[4] = pool_astrom_Th.dat_h[i_src] * pool_mag.dat_h[i_src].mag;
        }
        Nuniform[i_src] = 0;
    }
    for(int i_src=0;i_src<n_src;i_src++)
    {
        if(pool_mag.dat_h[i_src].err < pool_mag.dat_h[i_src].mag * RelTol
            && 0 >= NuniformMin)
        {
            ini_finished[i_src] = true;
            finished[i_src] = true;
        }
    }

    c_t* caustic_pts = new c_t[5*n_src];
    f_t* phi_hidden = new f_t[5*n_src];
    f_t prev_lens_s, prev_lens_q;
    for(int i_src=0;i_src<n_src;i_src++)
    {
        if(!finished[i_src])
        {
            c_t* caustic_pts_i = caustic_pts + i_src*5;
            if(i_src==0)
            {
                caustic_cal(pool_lens_s.dat_h[ i_src ], lens_q, caustic_pts_i);
                prev_lens_s = pool_lens_s.dat_h[ i_src ];
                prev_lens_q = lens_q;
            }
            else
            {
                if(pool_lens_s.dat_h[ i_src ] == prev_lens_s && lens_q == prev_lens_q)
                {
                    for(int caus_i=0;caus_i<5;caus_i++)
                    {
                        caustic_pts_i[caus_i] = caustic_pts[(i_src-1)*5 + caus_i];
                    }
                }
                else
                {
                    caustic_cal(pool_lens_s.dat_h[ i_src ], lens_q, caustic_pts_i);
                    prev_lens_s = pool_lens_s.dat_h[ i_src ];
                    prev_lens_q = lens_q;
                }
            }
            f_t* phi_hidden_i = phi_hidden + i_src*5;
            hidden_phi(caustic_pts_i, pool_center.dat_h [ i_src ].loc_centre, pool_center.dat_h [ i_src ].rho, phi_hidden_i);
            EU_buffer[i_src].phi_hidden = phi_hidden_i;
        }
    }

    for(int iphi : {0, 2})
    {
        for(int i_src=0; i_src<n_src; i_src++)
        {
            if(!finished[i_src])
            {
                pool_phi.dat_h[i_src] = EU_buffer[i_src].phi[iphi];
            }
            else
            {
                pool_phi.dat_h[i_src] = -1.;
            } 
        }
        f_dev.sync_all_streams();
        pool_phi.cp_h2d(f_dev);
        f_dev.sync_all_streams();
        lpar = std::make_tuple(dim3(n_bl), dim3(n_th), s_sh);
        f_dev.launch(solve_LD_sh, lpar, stream, (*this));
        f_dev.sync_all_streams();
        pool_mag.cp_d2h(f_dev);
        pool_Ncross.cp_d2h(f_dev);
        if(astrom)
        {
            pool_astrom_Th.cp_d2h(f_dev);
        }
        f_dev.sync_all_streams();
        for(int i_src=0; i_src<n_src; i_src++)
        {
            if(!finished[i_src])
            {
                ErrorUnit& eu = EU_buffer[i_src];
                double phi = eu.phi[iphi];
                double mag = pool_mag.dat_h[i_src].mag;
                eu.val[iphi] = mag * (1.0 - phi*phi);
                eu.raw_Mag[iphi] = mag;
                eu.Ncross[iphi] = pool_Ncross.dat_h[i_src];
                Nuniform[i_src] += 1;
                if(astrom)
                {
                    eu.astrom_X0[iphi] = pool_astrom_Th.dat_h[i_src] * mag;
                }
            }
        }
    }
    for(int i_src=0;i_src<n_src;i_src++)
    {
        ErrorUnit& eu = EU_buffer[i_src];
        monotest_single(eu.val_L(),eu.val_R(),eu.raw_Mag[0],eu.raw_Mag[4], eu.phi[2], eu.val_C(),eu.raw_Mag[2]);
    }

    double E1, E2, magL, magC, magR, mag_ini,E_ini;
    double coeff_ini = 0.1;
    for(int i_src=0; i_src<n_src; i_src++)
    {
        if(!finished[i_src])
        {
            ErrorUnit& eu = EU_buffer[i_src];
            magL = eu.raw_Mag[0];
            magC = eu.raw_Mag[2];
            magR = eu.raw_Mag[4];
            mag_ini = 1./6. * (eu.val_L() + 4*eu.val_C() + eu.val_R());
            mag_ini = ((1-LD_a)*magL + LD_a*mag_ini) / (1-LD_a/3);
            E2 = fabs(magL - 2*magC + magR);
            E1 = fabs(magL - magR);
            E_ini = (3*E2<E1) ? E1 : E2;
            E_ini *= coeff_ini;
            f_t* phi_hidden_i = phi_hidden + i_src*5;
            if (E_ini/mag_ini < RelTol * coeff_RelTol
            && 2 >= NuniformMin
            && eu.Ncross[0] == 0
            && eu.Ncross[2] == 0
            && phi_hidden_i[4] ==0)
            {
                finished[i_src] = true;
                ini_finished[i_src] = true;
                pool_mag.dat_h[i_src].mag = mag_ini;
                pool_mag.dat_h[i_src].err = E_ini;
                if(astrom)
                {
                    complex_t<double> Th_ini = 1./6.*(eu.astrom_X0[0]*(1.) + 4*eu.astrom_X0[2]*(0.75) + eu.astrom_X0[4]*(0.));
                    Th_ini = ((1-LD_a)*eu.astrom_X0[0] + LD_a*Th_ini) / (1-LD_a/3);
                    Th_ini /= mag_ini;
                    pool_astrom_Th.dat_h[i_src] = Th_ini;
                }
            }
        }
    }
    
    for(int iphi: {1,3})
    {
        for(int i_src=0; i_src<n_src; i_src++)
        {
            if(!finished[i_src])
            {
                pool_phi.dat_h[i_src] = EU_buffer[i_src].phi[iphi];
            }
            else
            {
                pool_phi.dat_h[i_src] = -1.;
            }            
        }
        f_dev.sync_all_streams();
        pool_phi.cp_h2d(f_dev);
        f_dev.sync_all_streams();
        lpar = std::make_tuple(dim3(n_bl), dim3(n_th), s_sh);
        f_dev.launch(solve_LD_sh, lpar, stream, (*this));
        f_dev.sync_all_streams();
        pool_mag.cp_d2h(f_dev);
        pool_Ncross.cp_d2h(f_dev);
        if(astrom)
        {
            pool_astrom_Th.cp_d2h(f_dev);
        }
        f_dev.sync_all_streams();
        for(int i_src=0; i_src<n_src; i_src++)
        {
            if(!finished[i_src])
            {
                ErrorUnit& eu = EU_buffer[i_src];
                eu.val[iphi] = pool_mag.dat_h[i_src].mag * (1.0 - eu.phi[iphi]*eu.phi[iphi]);
                eu.Ncross[iphi] = pool_Ncross.dat_h[i_src];
                Nuniform[i_src] += 1; 
                if(astrom)
                {
                    eu.astrom_X0[iphi] = pool_astrom_Th.dat_h[i_src] * pool_mag.dat_h[i_src].mag;

                }
            }
        }
    }
    for(int i_src=0;i_src<n_src;i_src++)
    {
        ErrorUnit& eu = EU_buffer[i_src];
        monotest(eu);
    }

    for(int i_src=0; i_src<n_src; i_src++)
    {
        if(!finished[i_src])
        {
            ErrorUnit& eu = EU_buffer[i_src];
            M0[i_src] = eu.val_L();
            eu.SimpsonC = (eu.R() - eu.L()) / 6.0 * (eu.val_L() + 4.0*eu.val_C() + eu.val_R());
            if(astrom)
            {
                X0_max_radius[i_src] = eu.astrom_X0[0];
                eu.SimpsonC_X = (eu.R() - eu.L()) / 6.0 * (eu.astrom_X0[0]*(1-eu.phi[0]*eu.phi[0]) + 4 * eu.astrom_X0[2]*(1-eu.phi[2]*eu.phi[2]) + eu.astrom_X0[4]*(1-eu.phi[4]*eu.phi[4]));
            }
            eu.GetError();
            
            LDInte[i_src] = eu.SimpsonL + eu.SimpsonR;
            if(astrom)
            {
                X_LD[i_src] = eu.SimpsonL_X + eu.SimpsonR_X;
            }
            Error_now[i_src] = eu.Error;
            lowerbound[i_src] = eu.lb;
            heaps[i_src].push(eu);

        }
    }

    bool all_empty = false;
    while(!all_empty)
    {
        all_empty = true;
        for(int i_src=0;i_src<n_src;i_src++)
        {
            if(Error_now[i_src]/lowerbound[i_src] < RelTolS || Nuniform[i_src]>=NuniformMax || ini_finished[i_src])
            {
                EU_buffer[i_src] = ErrorUnit(-2.,-1.);
                EU_buffer_R[i_src] = ErrorUnit(-2.,-1.);
                finished[i_src] = true;
            }
            else
            {
                all_empty = false;
                auto max_error = heaps[i_src].top();
                heaps[i_src].pop();

                auto [euL, euR] = max_error.Split();
                EU_buffer[i_src] = euL;
                EU_buffer_R[i_src] = euR;
                finished[i_src] = false;

                LDInte[i_src] -= (max_error.SimpsonL+max_error.SimpsonR);
                if(astrom)
                {
                    X_LD[i_src] -= (max_error.SimpsonL_X + max_error.SimpsonR_X);
                }
                Error_now[i_src] -= max_error.Error;
                lowerbound[i_src] -= max_error.lb;
            }
        }
        
        for(int iphi: {1,3})
        {
            for(int i_src=0; i_src<n_src; i_src++)
            {
                pool_phi.dat_h[i_src] = EU_buffer[i_src].phi[iphi];
            }
            f_dev.sync_all_streams();
            pool_phi.cp_h2d(f_dev);
            f_dev.sync_all_streams();
            lpar = std::make_tuple(dim3(n_bl), dim3(n_th), s_sh);
            f_dev.launch(solve_LD_sh, lpar, stream, (*this));
            f_dev.sync_all_streams();
            pool_mag.cp_d2h(f_dev);
            pool_Ncross.cp_d2h(f_dev);
            if(astrom)
            {
                pool_astrom_Th.cp_d2h(f_dev);
            }
            f_dev.sync_all_streams();
            for(int i_src=0; i_src<n_src; i_src++)
            {
                if(!finished[i_src])
                {
                    ErrorUnit& euL = EU_buffer[i_src];
                    euL.val[iphi] = pool_mag.dat_h[i_src].mag * (1.0 - euL.phi[iphi]*euL.phi[iphi]);
                    euL.Ncross[iphi] = pool_Ncross.dat_h[i_src];
                    Nuniform[i_src] += 1;
                    if(astrom)
                    {
                        euL.astrom_X0[iphi] = pool_astrom_Th.dat_h[i_src] * pool_mag.dat_h[i_src].mag;
                    }

                }
            }
            for(int i_src=0; i_src<n_src; i_src++)
            {
                pool_phi.dat_h[i_src] = EU_buffer_R[i_src].phi[iphi];
            }
            f_dev.sync_all_streams();
            pool_phi.cp_h2d(f_dev);
            f_dev.sync_all_streams();
            lpar = std::make_tuple(dim3(n_bl), dim3(n_th), s_sh);
            f_dev.launch(solve_LD_sh, lpar, stream, (*this));
            f_dev.sync_all_streams();
            pool_mag.cp_d2h(f_dev);
            pool_Ncross.cp_d2h(f_dev);
            if(astrom)
            {
                pool_astrom_Th.cp_d2h(f_dev);
            }
            f_dev.sync_all_streams();
            for(int i_src=0; i_src<n_src; i_src++)
            {
                if(!finished[i_src])
                {
                    ErrorUnit& euR = EU_buffer_R[i_src];
                    euR.val[iphi] = pool_mag.dat_h[i_src].mag * (1.0 - euR.phi[iphi]*euR.phi[iphi]);
                    euR.Ncross[iphi] = pool_Ncross.dat_h[i_src];
                    Nuniform[i_src] += 1;
                    if(astrom)
                    {
                        euR.astrom_X0[iphi] = pool_astrom_Th.dat_h[i_src] * pool_mag.dat_h[i_src].mag;
                    }
                }
            }
        }
        for(int i_src=0; i_src<n_src; i_src++)
        {
            if(!finished[i_src])
            {
                ErrorUnit& euL = EU_buffer[i_src];
                monotest(euL);
                euL.GetError();
                ErrorUnit& euR = EU_buffer_R[i_src];
                monotest(euR);
                euR.GetError();

                LDInte[i_src] += (euL.SimpsonL+euL.SimpsonR)+(euR.SimpsonL+euR.SimpsonR);
                if(astrom)
                    X_LD[i_src] += (euL.SimpsonL_X+euL.SimpsonR_X)+(euR.SimpsonL_X+euR.SimpsonR_X);
                Error_now[i_src] += euL.Error+euR.Error;
                lowerbound[i_src] += euL.lb+euR.lb;

                heaps[i_src].push(euL);
                heaps[i_src].push(euR);

                if(euL.depth>max_depth)
                {
                    Error_now[i_src] = 0.;
                }
            }
        }

    }


    for(int i_src=0;i_src<n_src;i_src++)
    {
        if(ini_finished[i_src] == false)
        {
            f_t return_mag = ((1-LD_a)*M0[i_src] + LD_a*LDInte[i_src]) / (1-LD_a/3);
            pool_mag.dat_h[ i_src ].mag = return_mag; 
            if(astrom)
            {
                complex_t<float2_t> return_X_LD = ((1-LD_a)*X0_max_radius[i_src] + LD_a*X_LD[i_src]) / (1-LD_a/3); 
                pool_astrom_Th.dat_h[ i_src ] = return_X_LD / return_mag;
            }
        }
    }
    pool_mag.cp_h2d( f_dev );
    if(astrom)
    {
        pool_astrom_Th.cp_h2d( f_dev );        
    }
    f_dev.sync_all_streams();

    if(Nuniform_out)
    {
        for(int i_src=0; i_src<n_src; i_src++)
        {
            Nuniform_out[i_src] = Nuniform[i_src];
        }
    }

    delete[] LDInte;
    delete[] M0;
    delete[] Error_now;
    delete[] lowerbound;
    delete[] finished;
    delete[] ini_finished;
    delete[] EU_buffer;
    delete[] EU_buffer_R;
    delete[] Nuniform;
    delete[] caustic_pts;
    delete[] phi_hidden;
    delete[] heaps;
    delete[] X_LD;
    delete[] X0_max_radius;

    return;
}

}; // namespace twinkle
