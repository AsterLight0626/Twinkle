#include "source_base.h"

namespace twinkle
{
////////////////////////////////////////////////////////////
//

__host__ source_base_t:: source_base_t(  )
{
    // n_th  = 32;
    // n_src = 64*2*2;
    // n_src = 1;
    // n_src = 256;
    n_src = 1;          // 在 driver 中调整

    // Temporarily hard-coded
    n_point_max = 4096;
    n_cross_max =  n_th;
    return;
}

__host__ source_base_t::~source_base_t(  )
{
    return;
}

__host__ void source_base_t::init( device_t & f_dev, const bool pt_only )
{
    if( !pt_only )
    {
        pool_margin.setup_mem( f_dev, n_src * n_point_max,
                                 n_point_max  );
        pool_cross .setup_mem( f_dev, n_src * n_cross_max,
                                 n_cross_max  );
    }
    pool_center.setup_mem( f_dev, n_src );
    pool_mag   .setup_mem( f_dev, n_src );
    pool_extended   .setup_mem( f_dev, n_src );
    pool_lens_s.setup_mem( f_dev, n_src );
    pool_phi.setup_mem( f_dev, n_src );
    pool_Ncross.setup_mem( f_dev, n_src);

    // if(astrom)
    // {
    pool_astrom_Th.setup_mem( f_dev, n_src);
    // }

    stream = f_dev.yield_stream(  );

    // s = ss;
    // q = qq;


    // RelTol    =   1e-4;
    // m1 = 1 / ( q + 1 );
    // m2 = q / ( q + 1 );
    // int Nx = 30;
    // int Ny = 30;

    // for(int y_idx=0;y_idx<Ny;y_idx++)
    // {
    //     for(int x_idx=0;x_idx<Nx;x_idx++)
    //     {
    //         auto & shape  = pool_center.dat_h[ x_idx + y_idx * Nx ];
    //         shape.rho = rho;
    //         shape.loc_centre.re = ((f_t(x_idx)+0.5) / Nx) * (xmax - xmin) + xmin;
    //         shape.loc_centre.im = ((f_t(y_idx)+0.5) / Ny) * (ymax - ymin) + ymin;            
    //     }
    // }


    // pool_center.cp_h2d( f_dev );
    // // std::cerr << "Finishing setting up initial condition"
    // //           << " for tests\n";
    return;
}

// __host__ void source_base_t::set_params_2D( device_t & f_dev, f_t ss, f_t qq, f_t rho, f_t xmax, f_t xmin, f_t ymax, f_t ymin, int Nx, int Ny )
// {
//     s = ss;
//     q = qq;

//     RelTol    =   1e-4;
//     m1 = 1 / ( q + 1 );
//     m2 = q / ( q + 1 );
//     // int Nx = 30;
//     // int Ny = 30;

//     for(int y_idx=0;y_idx<Ny;y_idx++)
//     {
//         for(int x_idx=0;x_idx<Nx;x_idx++)
//         {
//             auto & shape  = pool_center.dat_h[ x_idx + y_idx * Nx ];
//             shape.rho = rho;
//             shape.loc_centre.re = ((f_t(x_idx)+0.5) / f_t(Nx)) * (xmax - xmin) + xmin;
//             shape.loc_centre.im = ((f_t(y_idx)+0.5) / f_t(Ny)) * (ymax - ymin) + ymin;            
//         }
//     }

//     pool_center.cp_h2d( f_dev );
//     return;
// }

// __host__ void source_base_t::set_params_1D( device_t & f_dev, f_t ss, f_t qq, f_t rho, f_t xmax, f_t xmin, f_t ymax, f_t ymin, int Nsrc )
// {
//     s = ss;
//     q = qq;

//     RelTol    =   1e-4;
//     m1 = 1 / ( q + 1 );
//     m2 = q / ( q + 1 );
//     // int Nx = 30;
//     // int Ny = 30;

//     for(int idx=0;idx<Nsrc;idx++)
//     {

//         auto & shape  = pool_center.dat_h[ idx ];
//         shape.rho = rho;
//         shape.loc_centre.re = ((f_t(idx)) / f_t(Nsrc)) * (xmax - xmin) + xmin;
//         shape.loc_centre.im = ((f_t(idx)) / f_t(Nsrc)) * (ymax - ymin) + ymin;            

//     }

//     pool_center.cp_h2d( f_dev );
//     return;
// }


// __host__ void source_base_t::set_same( device_t & f_dev, f_t ss, f_t qq, f_t rho, f_t zeta_x, f_t zeta_y )
// {
//     for(int idx=0;idx<n_src;idx++)
//     {
//         pool_lens_s.dat_h[idx] = ss;
//     }    
//     // lens_s = ss;
//     lens_q = qq;

//     RelTol    =   1e-4;
//     m1 = 1 / ( lens_q + 1 );
//     m2 = lens_q / ( lens_q + 1 );

//     for(int idx=0;idx<n_src;idx++)
//     {
//         auto & shape  = pool_center.dat_h[ idx ];
//         shape.rho = rho;
//         shape.loc_centre.re = zeta_x;
//         shape.loc_centre.im = zeta_y;            
//     }


//     pool_center.cp_h2d( f_dev );
//     pool_lens_s.cp_h2d( f_dev );
//     return;
// }



__host__ void source_base_t::free( device_t & f_dev )
{
    pool_mag   .free( f_dev );    
    pool_cross .free( f_dev );
    pool_center.free( f_dev );
    pool_margin.free( f_dev );
    pool_extended.free( f_dev );
    pool_lens_s.free( f_dev );
    pool_phi.free( f_dev );
    pool_Ncross.free( f_dev );

    // if(astrom)
    // {
    pool_astrom_Th.free( f_dev );
    // }

    return;
}

__host__ void source_base_t::dump_margin( device_t & f_dev )
{
    pool_margin.cp_d2h( f_dev );
    return;
}

__host__ void source_base_t::run( device_t & f_dev )
{
    auto n_bl = ( n_src + n_th - 1 ) / n_th;
    auto s_sh = 0;

    auto lpar = std::make_tuple
        ( dim3( n_bl ), dim3( n_th ), s_sh );

    f_dev.launch( point_approximation, lpar, stream, ( * this ) );

    n_bl = n_src;
    s_sh = sizeof(int)*n_th + 
    std::max( sizeof(float2_t)*n_point_max, 
              sizeof(float2_t)   * 8*n_th
            + sizeof(shared_info_t< float2_t >) * 1 
            + sizeof(cross_info_t) * n_cross_max * 4
            + sizeof(src_pt_t< float2_t >) * n_th);
    // printf("shared memory size: %d\n",s_sh);

    lpar = std::make_tuple
        ( dim3( n_bl ), dim3( n_th ), s_sh );
    f_dev.launch( solve_extended_sh, lpar, stream, ( * this ) );


    return;
}

__host__ void source_base_t::run_pt( device_t & f_dev )
{
    auto n_bl = ( n_src + n_th - 1 ) / n_th;
    auto s_sh = 0;

    auto lpar = std::make_tuple
        ( dim3( n_bl ), dim3( n_th ), s_sh );

    f_dev.launch( point_approximation, lpar, stream, ( * this ) );

    // 这是 host 函数！我忘了，这里是没法访问 pool_extended 的，应该在 Err 上想办法
    // const int i_src = threadIdx.x + n_th * blockIdx.x;
    // if( i_src < n_src )
    //     printf("%d,",pool_extended[i_src].SolveSucceed);

    return;
}


}; // namespace twinkle
