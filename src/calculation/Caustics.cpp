#include"Caustics.h"

// 2 body only
template<isFloating f_T>
_Device void CritCoeff_2b(const complex_t<f_T>* params, const f_T& psi, complex_t<f_T>* res, bool print)
{
    // params[0]: z0
    // params[1]: z1
    // params[2]: m0
    // params[3]: m1
    complex_t<f_T> exp_ipsi;    // exp(-i*psi)
    complex_t<f_T> z0_2[3];     // (z-z0)**2
    complex_t<f_T> z1_2[3];     // (z-z1)**2


    if constexpr (sizeof(f_T)==8)
    {
        exp_ipsi.re = std::cos(psi);
        exp_ipsi.im = -std::sin(psi);
    }
    else
    {
        // exp_ipsi.re = __cosf(psi);
        // exp_ipsi.im = -__sinf(psi);
        exp_ipsi.re = std::cos(psi);
        exp_ipsi.im = -std::sin(psi);        
    }    

    z0_2[1] = params[0];
    z0_2[0] = z0_2[1]*z0_2[1];  // z0_2, 0 order coeff = z0**2
    z0_2[1] = z0_2[1] * f_T(-2.);
    z0_2[2] = complex_t<f_T>(f_T(1.));



    z1_2[1] = params[1];
    z1_2[0] = z1_2[1]*z1_2[1];  // z1_2, 0 order coeff = z1**2
    z1_2[1] = z1_2[1] * f_T(-2.);
    z1_2[2] = complex_t<f_T>(f_T(1.));

    if(print)
    {
        printf("(%.16f+%.16fj),(%.16f+%.16fj),(%.16f+%.16fj),(%.16f+%.16fj)\n",z1_2[0].re,z1_2[0].im,z1_2[1].re,z1_2[1].im,z0_2[0].re,z0_2[0].im,z0_2[1].re,z0_2[1].im);
    }   

    PolyMultiply(z0_2,2,z1_2,2,res);        // res = (z-z0)^2 * (z-z1)^2
    if(print)
    {
        printf("(%.16f+%.16fj),(%.16f+%.16fj),(%.16f+%.16fj),(%.16f+%.16fj),(%.16f+%.16fj)\n",res[0].re,res[0].im,res[1].re,res[1].im,res[2].re,res[2].im,res[3].re,res[3].im,res[4].re,res[4].im);
    }
    PolyMultiply(-params[3]*exp_ipsi,0,z0_2,2,z0_2);    // z0_2 = -m1*exp(-i*psi)*(z-z0)^2
    PolyAdd(z0_2,2,res,4,res);
    if(print)
    {
        printf("(%.16f+%.16fj),(%.16f+%.16fj),(%.16f+%.16fj),(%.16f+%.16fj),(%.16f+%.16fj)\n",res[0].re,res[0].im,res[1].re,res[1].im,res[2].re,res[2].im,res[3].re,res[3].im,res[4].re,res[4].im);
    }
    if(print)
    {
        printf("(%.16f+%.16fj),(%.16f+%.16fj),(%.16f+%.16fj),(%.16f+%.16fj)\n",z1_2[0].re,z1_2[0].im,z1_2[1].re,z1_2[1].im,z0_2[0].re,z0_2[0].im,z0_2[1].re,z0_2[1].im);
    } 
    PolyMultiply(-params[2]*exp_ipsi,0,z1_2,2,z1_2);    // z1_2 = -m0*exp(-i*psi)*(z-z1)^2
    PolyAdd(z1_2,2,res,4,res);              // res = (z-z0)^2 * (z-z1)^2 - m1*exp(-i*psi)*(z-z0)^2 - m0*exp(-i*psi)*(z-z1)^2
    if(print)
    {
        printf("(%.16f+%.16fj),(%.16f+%.16fj),(%.16f+%.16fj),(%.16f+%.16fj)\n",z1_2[0].re,z1_2[0].im,z1_2[1].re,z1_2[1].im,z0_2[0].re,z0_2[0].im,z0_2[1].re,z0_2[1].im);
    } 
    if(print)
    {
        printf("(%.16f+%.16fj),(%.16f+%.16fj),(%.16f+%.16fj),(%.16f+%.16fj),(%.16f+%.16fj)\n",res[0].re,res[0].im,res[1].re,res[1].im,res[2].re,res[2].im,res[3].re,res[3].im,res[4].re,res[4].im);
    }
    // highest coeff of res is 1

}
template _Device void CritCoeff_2b(const complex_t<float>* params, const float& psi, complex_t<float>* res, bool print);
template _Device void CritCoeff_2b(const complex_t<double>* params, const double& psi, complex_t<double>* res, bool print);


// if real solution: true; if additional solution: false
template<isFloating f_T>
_Device bool CritAddRootTst(const complex_t<f_T>& root, const complex_t<f_T>* params)
{
    bool real_root=true;
    for(int i=0;i<NBODY;i++)
    {
        if(norm(root-params[i])<2e-15)
        {
            real_root = false;
            // printf("%.12f,%.12f\n",root.re,root.im);
        }
    }
    return real_root;
}
template _Device bool CritAddRootTst(const complex_t<float>& root, const complex_t<float>* params);
template _Device bool CritAddRootTst(const complex_t<double>& root, const complex_t<double>* params);



// 重载，记录求解失败
template<isFloating f_T>
_Global void SolveCritCaus(const src_params_t<f_T>* src_params,complex_t<f_T>* CritLoc, complex_t<f_T>* CausLoc)
{
    int paramsidx = blockIdx.x;
    int crit_idx = threadIdx.x + blockIdx.z*blockDim.x;
    f_T params_lens[2];
    complex_t<f_T> params[2*NBODY];
    complex_t<f_T> coeff[LENGTH-1];
    f_T& psi = params_lens[0];
    complex_t<f_T> roots[NBODY*NBODY];
    bool fail;
    // bool single_real;

    if(crit_idx < NCRIT)
    {
        // 二体特例
        params_lens[0] = src_params[paramsidx].s;
        params_lens[1] = src_params[paramsidx].q;

        ParamsReal2Complex(params_lens, params);

        psi = (f_T(crit_idx) / f_T(NCRIT)) * 2* PI;

        CritCoeff_2b(params,psi,coeff,false);

        
        fail = cmplx_roots_gen(roots,coeff,LENGTH-2,true,false);
        if(! fail)
        {
            // printf("%d,(%.12f+%.12fj),(%.12f+%.12fj),(%.12f+%.12fj),(%.12f+%.12fj)\n",crit_idx,roots[0].re,roots[0].im,roots[1].re,roots[1].im,roots[2].re,roots[2].im,roots[3].re,roots[3].im);
            for(int j=0;j<NBODY*NBODY;j++)
            {
                CritLoc[crit_idx + j*NCRIT] = roots[j];
                CausLoc[crit_idx + j*NCRIT] = roots[j] + f_zbar(roots[j],params);


            }
        }
    }


}
template _Global void SolveCritCaus(const src_params_t<float>* params,complex_t<float>* CritLoc, complex_t<float>* CausLoc);
template _Global void SolveCritCaus(const src_params_t<double>* params,complex_t<double>* CritLoc, complex_t<double>* CausLoc);


// size  <<<(NWALKERS,1,NCRIT/64),(64)>>>
template<isFloating f_T>
_Global void CNext(complex_t<f_T>* cc, int* next_j)
{
    int idx = threadIdx.x + blockIdx.z*blockDim.x;
    f_T dist2[NBODY*NBODY];
    int rank[NBODY*NBODY];
    if(idx < NCRIT)
    {
        for(int jhere=0;jhere<NBODY*NBODY;jhere++)
        {
            for(int jnext=0;jnext<NBODY*NBODY;jnext++)
            {
                rank[jnext] = jnext;
                dist2[jnext] = norm(cc[idx + jhere*NCRIT] - cc[(idx+1)%NCRIT + jnext*NCRIT] );
            }
            
            BubbleArgSort(rank,dist2,NBODY*NBODY);
            // printf("rank:%d,%d,%d,%d\n",rank[0],rank[1],rank[2],rank[3]);
            next_j[idx + jhere*NCRIT] = rank[0];
            // printf("idx: %d, jhere, %d, jnext, %d\n",idx, jhere, rank[0]);
        }        
    }

    for(int jhere=0;jhere<NBODY*NBODY;jhere++)
    {

        // printf("idx: %d, jhere, %d, jnext, %d\n",idx, jhere, next_j[idx + jhere*NCRIT]);
    }   

    // if(true)
    // {
    //     printf("S");
    // }

}
template _Global void CNext(complex_t<float>* caus, int* next_j);
template _Global void CNext(complex_t<double>* caus, int* next_j);

// 2 * triangle area
template<isFloating f_T>
_Host_Device f_T TriAreaD(f_T x0, f_T y0, f_T x1, f_T y1, f_T x2, f_T y2)
{
    // (x0-x2,y0-y2) cross (x1-x2,y1-y2)
    return (x0-x2)*(y1-y2) - (x1-x2)*(y0-y2);
}
template _Host_Device float TriAreaD(float x0, float y0, float x1, float y1, float x2, float y2);
template _Host_Device double TriAreaD(double x0, double y0, double x1, double y1, double x2, double y2);

// // size<<< (NWALKERS,1,NSRCS),(NCRIT)>>>
// // shared size: int*64 + complex<f_T>*1024*4
// // 暂不考虑x轴对称，考虑的话能减一半比较
// template<isFloating f_T>
// _Global void CollideTest(src_ext_t<f_T>* srcs, const complex_t<f_T>* caus, const int* next_j, int batchidx)
// {
//     int srcidx = blockIdx.z;
//     if(srcidx<NSRCS){
//         if(srcs[srcidx].SolveSucceed == false){
//             extern __shared__ int shared_mem[];
//             int* Ncross = (int*) shared_mem;
//             complex_t<f_T>* causLoc = (complex_t<f_T>*) &shared_mem[BATCH_SIZE];
//             int idx = threadIdx.x;              // 0 ~ 1024
            
//             complex_t<f_T> const * loc0;
//             complex_t<f_T> const * loc1;
//             complex_t<f_T> zeta;
//             int temp;

//             for(int i=0;i<BATCH_SIZE;i++){Ncross[i]=0;}      // 0 ~ 64
//             for(int j=0;j<NBODY*NBODY;j++)
//             {
//                 causLoc[idx+j*NCRIT] = caus[idx+j*NCRIT];       // input to shared, save the reload of caus_loc_1
//                 // printf("(%.12f+%.12fj)\n",caus[idx+j*NCRIT].re,caus[idx+j*NCRIT].im);
//             }

//             // printf("check %d\n",idx);
//             for(int j=0;j<NBODY*NBODY;j++)    // Caustic 有4段
//             {
//                 loc0 = &(causLoc[idx + j*NCRIT]);
//                 loc1 = &(causLoc[(idx+1)%NCRIT + next_j[idx + j*NCRIT]*NCRIT]);
//                 // printf("idx:%d, j:%d, nextj:%d\n",idx,j,next_j[idx + j*NCRIT]);
//                 // printf("%d,%d\n",idx,Ncross[idx]);
//                 for(int i=0;i<BATCH_SIZE;i++)
//                 {
//                     zeta = srcs[srcidx].margin_pts[i + batchidx*BATCH_SIZE].position;
//                     zeta = complex_t<f_T>(zeta.re,abs(zeta.im));                // 所有源点翻到x轴以上
//                     if((zeta.re<(*loc0).re && zeta.re>(*loc1).re)||(zeta.re>(*loc0).re && zeta.re<(*loc1).re)     )// 初步测试，寻找有没有两个点，其x坐标把 i 点夹在中间
//                     // 原则上需要判断==的时候（算作一半），但是概率原则上为0，就算了
//                     {
//                         // 尽量减少计算
//                         // 只看上半截
//                         if((*loc0).im>=0 && (*loc1).im>=0)
//                         {
                            
//                             if(zeta.im<=(*loc0).im && zeta.im<=(*loc1).im)   // 如果比两个点都低，那在下方
//                             {
//                                 temp = atomicAdd(&Ncross[i],1);
//                                 // printf("0 %d,%d,%d,(%.12f+%.12fj),(%.12f+%.12fj),(%.12f+%.12fj),%d\n",i + batchidx*BATCH_SIZE,idx+j*NCRIT,(idx+1)%NCRIT + next_j[idx + j*NCRIT]*NCRIT,zeta.re,zeta.im,(*loc0).re,(*loc0).im,(*loc1).re,(*loc1).im,temp);
//                             }
//                             else{
//                                 if(zeta.im>(*loc0).im || zeta.im>(*loc1).im)    // 在两个点为对角线的矩形之中
//                                 {
                                    
//                                     if(TriAreaD((*loc0).re,(*loc0).im,(*loc1).re,(*loc1).im,zeta.re,zeta.im) >0)    // 012正面积，封闭曲线内部（因为NCRIT是右手的）
//                                     {
//                                         if((*loc0).re>(*loc1).re)       // 0 right, 1 left, area positive ==> 2 below
//                                         {
//                                             temp = atomicAdd(&Ncross[i],1);
//                                             // printf("1 %d,%d,(%.12f+%.12fj),(%.12f+%.12fj),(%.12f+%.12fj),%d\n",i + batchidx*BATCH_SIZE,idx+j*NCRIT,zeta.re,zeta.im,(*loc0).re,(*loc0).im,(*loc1).re,(*loc1).im,temp);
//                                         }
//                                     }
//                                     else
//                                     {
//                                         if((*loc0).re<(*loc1).re)       // 0 left, 1 right, area negative ==> 2 below
//                                         {
//                                             temp = atomicAdd(&Ncross[i],1);
//                                             // printf("2 %d,%d,(%.12f+%.12fj),(%.12f+%.12fj),(%.12f+%.12fj),%d\n",i + batchidx*BATCH_SIZE,idx+j*NCRIT,zeta.re,zeta.im,(*loc0).re,(*loc0).im,(*loc1).re,(*loc1).im,temp);
//                                         }                               
//                                     }
//                                 }
//                                 // else{}    // 比两个点都高，点在散焦线上方（不与向上的射线相交），不做操作
                                
//                             }
//                         }
//                     }
//                 }
//             }
//             // 全检测完了
//             __syncthreads();
//             if(idx<BATCH_SIZE)
//             {
//                 // printf("%d,%d\n",idx + batchidx*BATCH_SIZE,Ncross[idx]);
//                 if((Ncross[idx] & 1)==1)    // mod 2 == 1, inside
//                 {
//                     srcs[srcidx].margin_pts[idx + batchidx*BATCH_SIZE].Nphys = 5;
//                 }
//                 else        // outside
//                 {
//                     srcs[srcidx].margin_pts[idx + batchidx*BATCH_SIZE].Nphys = 3;
//                 }
//                 srcs[srcidx].margin_pts[idx + batchidx*BATCH_SIZE].NcrossCaus = Ncross[idx];
//             }
//         }
//     }
// }
// template _Global void CollideTest(src_ext_t<float>* srcs, const complex_t<float>* caus, const int* next_j, int batchidx);
// template _Global void CollideTest(src_ext_t<double>* srcs, const complex_t<double>* caus, const int* next_j, int batchidx);


// size<<< (1,1,NSRCS),(NCRIT)>>>
// shared size: int*64
// 暂不考虑x轴对称，考虑的话能减一半比较
template<isFloating f_T>
_Global void CollideTest(int* srclist, src_ext_t<f_T>* srcs, const complex_t<f_T>* caus, const int* next_j, int batchidx)
{
    int srcidx = srclist[int(blockIdx.z)];

    // extern __shared__ int shared_mem[];
    __dyn_shared__(int, shared_mem);

    int* Ncross = (int*) shared_mem;
    // complex_t<f_T>* causLoc = (complex_t<f_T>*) &shared_mem[BATCH_SIZE];
    int idx = threadIdx.x;              // 0 ~ 1024
    
    // complex_t<f_T> const * loc0;
    // complex_t<f_T> const * loc1;
    complex_t<f_T> zeta;
    // int temp;
    complex_t<f_T> posi0, posi1;

    for(int i=0;i<BATCH_SIZE;i++){Ncross[i]=0;}      // 0 ~ 64
    // for(int j=0;j<NBODY*NBODY;j++)
    // {
    //     causLoc[idx+j*NCRIT] = caus[idx+j*NCRIT];       // input to shared, save the reload of caus_loc_1
    //     // printf("(%.12f+%.12fj)\n",caus[idx+j*NCRIT].re,caus[idx+j*NCRIT].im);
    // }

    // printf("check %d\n",idx);
    for(int j=0;j<NBODY*NBODY;j++)    // Caustic 有4段
    {
        // loc0 = &(causLoc[idx + j*NCRIT]);
        // loc1 = &(causLoc[(idx+1)%NCRIT + next_j[idx + j*NCRIT]*NCRIT]);
        posi0 = caus[idx+j*NCRIT];
        posi1 = caus[(idx+1)%NCRIT + next_j[idx + j*NCRIT]*NCRIT];


        // printf("idx:%d, j:%d, nextj:%d\n",idx,j,next_j[idx + j*NCRIT]);
        // printf("%d,%d\n",idx,Ncross[idx]);
        for(int i=0;i<BATCH_SIZE;i++)
        {
            zeta = srcs[srcidx].margin_pts[i + batchidx*BATCH_SIZE].position;
            zeta = complex_t<f_T>(zeta.re,fabs(zeta.im));                // 所有源点翻到x轴以上
            if((zeta.re<posi0.re && zeta.re>posi1.re)||(zeta.re>posi0.re && zeta.re<posi1.re)     )// 初步测试，寻找有没有两个点，其x坐标把 i 点夹在中间
            // 原则上需要判断==的时候（算作一半），但是概率原则上为0，就算了
            {
                // 尽量减少计算
                // 只看上半截
                if(posi0.im>=-2e-12 && posi1.im>=-2e-12)        // 需要考虑一定的数值误差
                {
                    
                    if(zeta.im<=posi0.im && zeta.im<=posi1.im)   // 如果比两个点都低，那在下方
                    {
                        // temp = atomicAdd(&Ncross[i],1);
                        atomicAdd(&Ncross[i],1);
                        // printf("0 %d,%d,%d,(%.12f+%.12fj),(%.12f+%.12fj),(%.12f+%.12fj),%d\n",i + batchidx*BATCH_SIZE,idx+j*NCRIT,(idx+1)%NCRIT + next_j[idx + j*NCRIT]*NCRIT,zeta.re,zeta.im,posi0.re,posi0.im,posi1.re,posi1.im,temp);
                    }
                    else{
                        if(zeta.im>posi0.im || zeta.im>posi1.im)    // 在两个点为对角线的矩形之中
                        {
                            
                            if(TriAreaD(posi0.re,posi0.im,posi1.re,posi1.im,zeta.re,zeta.im) >0)    // 012正面积，封闭曲线内部（因为NCRIT是右手的）
                            {
                                if(posi0.re>posi1.re)       // 0 right, 1 left, area positive ==> 2 below
                                {
                                    // temp = atomicAdd(&Ncross[i],1);
                                    atomicAdd(&Ncross[i],1);
                                    // printf("1 %d,%d,(%.12f+%.12fj),(%.12f+%.12fj),(%.12f+%.12fj),%d\n",i + batchidx*BATCH_SIZE,idx+j*NCRIT,zeta.re,zeta.im,posi0.re,posi0.im,posi1.re,posi1.im,temp);
                                }
                            }
                            else
                            {
                                if(posi0.re<posi1.re)       // 0 left, 1 right, area negative ==> 2 below
                                {
                                    // temp = atomicAdd(&Ncross[i],1);
                                    atomicAdd(&Ncross[i],1);
                                    // printf("2 %d,%d,(%.12f+%.12fj),(%.12f+%.12fj),(%.12f+%.12fj),%d\n",i + batchidx*BATCH_SIZE,idx+j*NCRIT,zeta.re,zeta.im,posi0.re,posi0.im,posi1.re,posi1.im,temp);
                                }                               
                            }
                        }
                        // else{}    // 比两个点都高，点在散焦线上方（不与向上的射线相交），不做操作
                        
                    }
                }
            }
        }
    }
    // 全检测完了
    __syncthreads();
    if(idx<BATCH_SIZE)
    {
        // printf("%d,%d\n",idx + batchidx*BATCH_SIZE,Ncross[idx]);
        if((Ncross[idx] & 1)==1)    // mod 2 == 1, inside
        {
            srcs[srcidx].margin_pts[idx + batchidx*BATCH_SIZE].Nphys = 5;
        }
        else        // outside
        {
            srcs[srcidx].margin_pts[idx + batchidx*BATCH_SIZE].Nphys = 3;
        }
        srcs[srcidx].margin_pts[idx + batchidx*BATCH_SIZE].NcrossCaus = Ncross[idx];
    }
}

template _Global void CollideTest(int* srclist, src_ext_t<float>* srcs, const complex_t<float>* caus, const int* next_j, int batchidx);
template _Global void CollideTest(int* srclist, src_ext_t<double>* srcs, const complex_t<double>* caus, const int* next_j, int batchidx);



// size<<< (1,1,NSRCS),(NCRIT)>>>
// shared size: int* 4 * BATCH_SIZE (只看出buried的最大和最小，如果出现多个交替也不管)
// 暂不考虑x轴对称，考虑的话能减一半比较
// 一个线程对应一个焦散线上的区间，考虑任意方向
template<isFloating f_T>
_Global void CollideTest2(int* srclist, src_ext_t<f_T>* srcs, const complex_t<f_T>* caus, const int* next_j, int batchidx)
{
    int srcidx = srclist[int(blockIdx.z)];

    // extern __shared__ int shared_mem[];
    __dyn_shared__(int, shared_mem);
    int* Ncross = (int*) shared_mem;
    int* NcrossNext = (int*) &shared_mem[BATCH_SIZE];
    int* loc_min = (int*) &shared_mem[2*BATCH_SIZE];
    int* loc_max = (int*) &shared_mem[3*BATCH_SIZE];
    // complex_t<f_T>* causLoc = (complex_t<f_T>*) &shared_mem[BATCH_SIZE];
    int idx = threadIdx.x;              // 0 ~ 1024
    
    // complex_t<f_T> const * loc0;
    // complex_t<f_T> const * loc1;
    complex_t<f_T> zeta0,zeta1;
    int temp,temp_idx;
    complex_t<f_T> posi0, posi1;
    // complex_t<float> new0,new1;
    complex_t<f_T> new0,new1;
    bool skip;

    f_T& tempf = zeta1.re;

    f_T right_idx_Q, left_idx_Q;

    if(idx<BATCH_SIZE)
    {
        Ncross[idx] = 0;
        NcrossNext[idx] = 0;
        loc_min[idx] = 100000000;
        loc_max[idx] = 0;
    }
    if(idx==0)
    {
        srcs[srcidx].NBuried = 0;
    }
    __syncthreads();



    // for(int j=0;j<NBODY*NBODY;j++)
    // {
    //     causLoc[idx+j*NCRIT] = caus[idx+j*NCRIT];       // input to shared, save the reload of caus_loc_1
    //     // printf("(%.12f+%.12fj)\n",caus[idx+j*NCRIT].re,caus[idx+j*NCRIT].im);
    // }

    // printf("check %d\n",idx);
    for(int j=0;j<NBODY*NBODY;j++)    // Caustic 有4段
    {
        // loc0 = &(causLoc[idx + j*NCRIT]);
        // loc1 = &(causLoc[(idx+1)%NCRIT + next_j[idx + j*NCRIT]*NCRIT]);
        posi0 = caus[idx+j*NCRIT];
        posi1 = caus[(idx+1)%NCRIT + next_j[idx + j*NCRIT]*NCRIT];


        // printf("idx:%d, j:%d, nextj:%d\n",idx,j,next_j[idx + j*NCRIT]);
        // printf("%d,%d\n",idx,Ncross[idx]);
        for(int i=0;i<BATCH_SIZE;i++)
        {
            zeta0 = srcs[srcidx].margin_pts[i + batchidx*BATCH_SIZE].position;
            temp_idx = srcs[srcidx].margin_pts[i + batchidx*BATCH_SIZE].next_src_idx;
            skip = srcs[srcidx].margin_pts[temp_idx].skip;
            while(skip==true)
            {
                temp_idx = srcs[srcidx].margin_pts[temp_idx].next_src_idx;
                skip = srcs[srcidx].margin_pts[temp_idx].skip;
            }
            zeta1 = srcs[srcidx].margin_pts[temp_idx].position;
            // Delta_zeta = srcs[srcidx].margin_pts[temp].position - srcs[srcidx].margin_pts[i + batchidx*BATCH_SIZE].position;

            // 将zeta0归位原点，zeta1归为（1，0）
            zeta1 = zeta1 - zeta0;

            // new0 = complex_t<float>(posi0-zeta0) / complex_t<float>(zeta1);
            // new1 = complex_t<float>(posi1-zeta0) / complex_t<float>(zeta1);
            new0 = (posi0-zeta0) / (zeta1);         // 归一化的焦散线点
            new1 = (posi1-zeta0) / (zeta1);

            if((new0.im>=0) != (new1.im>0)) // 如果转动参考系之后，锚点的y坐标异号
            {
                // 如果两个参考点都比1大
                if((new0.re > 1.)&&(new1.re > 1.))
                {
                    temp = atomicAdd(&Ncross[i],1);
                    temp = atomicAdd(&NcrossNext[i],1);
                }
                else
                {
                    // 如果参考点既不是特别大也不是特别小
                    if( !( (new0.re < 0.) && (new1.re < 0.) ) )
                    {
                        // tempf is temporary valuable
                        tempf = (new1.im*new0.re - new0.im*new1.re) / (new1.im - new0.im);      // 焦散线连线与x轴点交点
                        if(tempf > 1.)
                        {
                            temp = atomicAdd(&Ncross[i],1);
                            temp = atomicAdd(&NcrossNext[i],1);
                        }
                        else{
                            if(tempf > 0.)
                            {
                                temp = atomicAdd(&Ncross[i],1);
                                temp = atomicMin(&loc_min[i],int(tempf*1e8));
                                temp = atomicMax(&loc_max[i],int(tempf*1e8));
                            }
                            // else{}
                        }
                    }
                    // 如果两个参考点都比0小，什么都不做（+0）
                    // else{}
                    
                }
            }

        }
    }
    // 全检测完了
    __syncthreads();
    if(idx<BATCH_SIZE)
    {
        // printf("%d,%d\n",idx + batchidx*BATCH_SIZE,Ncross[idx]);
        if((Ncross[idx] & 1)==1)    // % 2 == 1， inside caustic
        {
            srcs[srcidx].margin_pts[idx + batchidx*BATCH_SIZE].Nphys = 5;
        }
        else        // outside
        {
            srcs[srcidx].margin_pts[idx + batchidx*BATCH_SIZE].Nphys = 3;
        }
        srcs[srcidx].margin_pts[idx + batchidx*BATCH_SIZE].NcrossCaus = Ncross[idx];



        if(abs(Ncross[idx]-NcrossNext[idx]) > 1)
        {
            temp = atomicAdd(&(srcs[srcidx].NBuried),1);
            srcs[srcidx].idx_buried[temp] = idx + batchidx*BATCH_SIZE;
            temp_idx = srcs[srcidx].margin_pts[idx + batchidx*BATCH_SIZE].next_src_idx;
            skip = srcs[srcidx].margin_pts[temp_idx].skip;
            while(skip==true)
            {
                temp_idx = srcs[srcidx].margin_pts[temp_idx].next_src_idx;
                skip = srcs[srcidx].margin_pts[temp_idx].skip;
            }

            left_idx_Q = srcs[srcidx].margin_pts[idx + batchidx*BATCH_SIZE].Q;
            right_idx_Q = srcs[srcidx].margin_pts[temp_idx].Q;
            if(right_idx_Q < left_idx_Q)
            {
                right_idx_Q += 1.;
            }
            tempf = left_idx_Q + (loc_min[idx]) / f_T(1e8) * (right_idx_Q - left_idx_Q);
            srcs[srcidx].Qmin_buried[temp] = tempf - int(tempf);
            tempf = left_idx_Q + (loc_max[idx]) / f_T(1e8) * (right_idx_Q - left_idx_Q);
            srcs[srcidx].Qmax_buried[temp] = tempf - int(tempf);
            // printf("(%d,%d), buried caustic, Qmin: %.12f, Qmax: %.12f\n",srcidx,idx + batchidx*BATCH_SIZE,srcs[srcidx].Qmin[temp], srcs[srcidx].Qmax[temp]);
        }
    }
}
template _Global void CollideTest2(int* srclist, src_ext_t<float>* srcs, const complex_t<float>* caus, const int* next_j, int batchidx);
template _Global void CollideTest2(int* srclist, src_ext_t<double>* srcs, const complex_t<double>* caus, const int* next_j, int batchidx);


// size  <<<(1,1,NCRIT/32),min(32,NCRIT)>>>
template<isFloating f_T>
_Global void SolveCritCaus(const src_params_t<f_T>* src_params, CC_t<f_T>* CC, int batchidx)
{
    int idx = threadIdx.x + blockIdx.z*blockDim.x + batchidx*BATCH_SIZE;
    f_T params_lens[2];
    complex_t<f_T> params[2*NBODY];
    complex_t<f_T> coeff[LENGTH-1];
    f_T& psi = params_lens[0];
    complex_t<f_T> roots[NBODY*NBODY];
    bool fail;

     if(idx < NCRIT)
    {
        // 二体特例
        params_lens[0] = src_params->s;
        params_lens[1] = src_params->q;

        ParamsReal2Complex(params_lens, params);

        psi = CC[idx].psi;

        CritCoeff_2b(params,psi,coeff,false);



        fail = cmplx_roots_gen(roots,coeff,LENGTH-2,true,false);
        if(! fail)
        {
            // printf("%d,(%.12f+%.12fj),(%.12f+%.12fj),(%.12f+%.12fj),(%.12f+%.12fj)\n",idx,roots[0].re,roots[0].im,roots[1].re,roots[1].im,roots[2].re,roots[2].im,roots[3].re,roots[3].im);
            for(int j=0;j<NBODY*NBODY;j++)
            {
                CC[idx].crit[j] = roots[j];


                CC[idx].caus[j] = roots[j] + f_zbar(roots[j],params);


            }
        }

    }
}
template _Global void SolveCritCaus(const src_params_t<float>* src_params, CC_t<float>* CC, int batchidx);
template _Global void SolveCritCaus(const src_params_t<double>* src_params, CC_t<double>* CC, int batchidx);


template<isFloating f_T>
_Device void ConnectOrderC(complex_t<f_T>* posi_here, complex_t<f_T>* posi_other, int* order, f_T& dist, bool muted=true)
{

    f_T ll[16];			// line length
    f_T Distance;
    f_T min_Distance = 1e99;
    // int min_here_i, min_next_i;
    // int Cantor_num;
    // int rank[4];
    int temp_j[4];
    
    int temp;
    int set;    // 0b1111

    // min_ll = 1e99;
    for(int here_i=0;here_i<4;here_i++)
    {
        for(int next_i=0;next_i<4;next_i++)
        {
            ll[here_i*4 + next_i] = abs_c(posi_here[here_i] - posi_other[next_i]);
        }
    }
    // Cantor Expansion
    for(int Cantor_i=0;Cantor_i<24;Cantor_i++)
    {
        set = 15;
        Distance = 0;
        temp_j[0] = Cantor_i / 6;			// Cantor / 3!
        temp = Cantor_i - temp_j[0]*6;		// Cantor % 3!
        temp_j[1] = temp / 2;				// temp / 2!
        temp -= temp_j[1]*2;					// temp % 2!
        temp_j[2] = temp;					// temp / 1!
                                                // temp % 1!
        temp_j[3] = 0;						// temp / 0!
        
        set -= (1<<temp_j[0]);        // 去除一个元素, << 表示左移，例如元素 2 在集合中表示为 0b0100 = 1<<2 , 去除即为 0b1111 - 0b0100 = 0b1011
        if(temp_j[0] <= temp_j[1])
        {
            temp_j[1] += 1;
        }
        set -= (1<<temp_j[1]);

        for(temp=0;temp<4;temp++)
        {
            // if(!muted){printf("1<<temp:%d, set&(1<<temp)";}
            if((set&(1<<temp))==0)  // 如果 temp 不在集合中
            {
                
                continue;       // 跳过
            }
            else
            {
                if(temp_j[2]==0)        // 用temp_j[2]做标记，正好是0说明到地方了
                {
                    temp_j[2] = temp;
                    break;
                }
                else
                {
                    temp_j[2] -= 1;     // 还没到0 （是1），等下一个符合的点 
                }
            }
            if(temp==3&&(!muted))
            {
                printf("cantor fail\n");
            }
        }

        temp_j[3] = 6 - temp_j[0] - temp_j[1] - temp_j[2];


        if(!muted)
        {
            printf("cantor: %d, sequence: %d%d%d%d\n",Cantor_i,temp_j[0],temp_j[1],temp_j[2],temp_j[3]);
        }

        for(int here_i=0;here_i<4;here_i++)
        {
            Distance += ll[ here_i*4 + temp_j[here_i]];
        }
        if(Distance<min_Distance)
        {
            min_Distance = Distance;
            for(int i=0;i<4;i++)
            {
                order[i] = temp_j[i];
                dist = Distance;
            }
        }
    }

}
template _Device void ConnectOrderC(complex_t<float>* posi_here, complex_t<float>* posi_other, int* order, float& dist, bool muted);
template _Device void ConnectOrderC(complex_t<double>* posi_here, complex_t<double>* posi_other, int* order, double& dist, bool muted);


// size  <<<(1,1,NCRIT/32),min(32,NCRIT)>>>
template<isFloating f_T>
_Global void CNext2(CC_t<f_T>* CC, int batchidx)
{
    int idx = threadIdx.x + blockIdx.z*blockDim.x;
    f_T distance;
    int next_j[NBODY*NBODY];
    int next_idx;
    bool skip;

    if(idx < NCRIT)
    {
        if(!CC[idx].skip)
        {
            next_idx = CC[idx].next_idx;
            skip = CC[next_idx].skip;
            while(skip)
            {
                next_idx = CC[next_idx].next_idx;
                skip = CC[next_idx].skip;
            }
            
            // if(idx==0)
            // {
            //     ConnectOrderC(CC[idx].crit, CC[next_idx].crit,next_j,distance,false);
            // }

            ConnectOrderC(CC[idx].crit, CC[next_idx].crit,next_j,distance);
            if(next_j[0]+next_j[1]+next_j[2]+next_j[3] != 6){printf("error next_j! crit idx: %d\n",idx);};


            for(int j=0;j<NBODY*NBODY;j++)
            {
                CC[idx].next_j[j] = next_j[j];               
            }
            // CC[idx].dist2next = distance;





            // for(int jhere=0;jhere<NBODY*NBODY;jhere++)
            // {
            //     for(int jnext=0;jnext<NBODY*NBODY;jnext++)
            //     {
            //         rank[jnext] = jnext;
            //         // dist2[jnext] = norm(cc[idx + jhere*NCRIT] - cc[(idx+1)%NCRIT + jnext*NCRIT] );
            //         dist2[jnext] = abs_c(CC[idx].crit[jhere] - CC[next_idx].crit[jnext]);
            //     }
                
            //     BubbleArgSort(rank,dist2,NBODY*NBODY);
            //     // printf("rank:%d,%d,%d,%d\n",rank[0],rank[1],rank[2],rank[3]);
            //     // next_j[idx + jhere*NCRIT] = rank[0];
            //     CC[idx].next_j[jhere] = rank[0];
            //     CC[idx].dist2next[jhere] = dist2[0];
            //     // printf("idx: %d, jhere, %d, jnext, %d\n",idx, jhere, rank[0]);
            // }    
        }    
    }


    // if(true)
    // {
    //     printf("S");
    // }

}
template _Global void CNext2(CC_t<float>* CC, int batchidx);
template _Global void CNext2(CC_t<double>* CC, int batchidx);



// size  <<<(1,1,1),max(batchidx,1)*BATCH_SIZE>>>
template<isFloating f_T>
_Global void psiLoc(CC_t<f_T>* CC, int batchidx)
{
    __dyn_shared__(int, shared_data);
    f_T* Weight = (f_T*) shared_data;                  // size = BATCH_SIZE*(batch_idx)
    int* ParentIdx = (int*) &shared_data[BATCH_SIZE*batchidx*sizeof(f_T)/sizeof(int)];


    int idx = threadIdx.x;
    f_T psi;
    int next_idx;
    int piece_length;

    f_T quantile;
    int left_idx;
    int right_idx;
    f_T left_psi;
    f_T right_psi;
    f_T quantile_left;
    f_T normalize;
    bool skip;

    if(batchidx==0)
    {

        psi = (f_T(idx) / f_T(BATCH_SIZE)) * 2* PI;
        next_idx = (idx+1)%BATCH_SIZE;
        CC[idx].psi = psi;

        for(int j=0;j<NBODY*NBODY;j++)
        {
            CC[idx].next_idx = next_idx;
        }              

    }
    else{

        // printf("trigger\n");


        Weight[idx] = 0.;

        if(!CC[idx].skip)
        {
            next_idx = CC[idx].next_idx;
            skip = CC[next_idx].skip;
            while(skip)
            {
                next_idx = CC[next_idx].next_idx;
                skip = CC[next_idx].skip;
            }


            for(int j=0;j<NBODY*NBODY;j++)
            {
                // Weight[idx] += CC[idx].dist2next[j];
                Weight[idx] += abs_c(CC[idx].caus[j]-CC[next_idx].caus[CC[idx].next_j[j]]);
            }            
        }

        // if(!CC[idx].skip)
        // {
        //     Weight[idx] = CC[idx].dist2next;
        // }
        __syncthreads();

        // if(idx==0)
        // {
        //     printf("\n");
        //     for(int iii = 0;iii<BATCH_SIZE*batchidx;iii++)
        //     {
        //         printf("%.12f,",Weight[iii]);
        //     }
        //     printf("\n");
        // }

        // 累加
        piece_length = 2;
        while(piece_length<2*batchidx*BATCH_SIZE)
        {
            if(((idx & (piece_length-1))>=piece_length/2) && (idx<BATCH_SIZE*batchidx))         // 二分的后半截，而且没超界
            {
                // if(idx==31)
                // {
                //     printf("plus idx: %d\n",idx - (idx&(piece_length-1)));
                // }
                Weight[idx] += Weight[idx - (idx&(piece_length-1)) + piece_length/2 -1];      // 加上前半截的尾巴， 和 AdapLoc 一致
            }
            __syncthreads();
            piece_length *= 2;
        }

        // if(idx==0)
        // {
        //     printf("\n");
        //     for(int iii = 0;iii<BATCH_SIZE*batchidx;iii++)
        //     {
        //         printf("%.12f,",Weight[iii]);
        //     }
        //     printf("\n");
        // }

        // 归一化
        normalize = Weight[BATCH_SIZE*batchidx - 1];
        __syncthreads();

        Weight[idx] /= normalize;
        __syncthreads();


        // if(idx==0)
        // {
        //     printf("\n");
        //     printf("norm: %.12f\n",normalize);
        //     for(int iii = 0;iii<BATCH_SIZE*batchidx;iii++)
        //     {
        //         printf("%.12f,",Weight[iii]);
        //     }
        //     printf("\n");
        // }

        // 插值
        if(idx < BATCH_SIZE)
        {
            // if(idx==0)
            // {
            //     printf("triggrt\n");
            // }

            quantile = (f_T(idx)+f_T(0.5)) / f_T(BATCH_SIZE);
            left_idx = BisearchLeft(Weight,quantile,batchidx*BATCH_SIZE);      // 取值范围：[0,len-1]，即parent节点，相邻的left_idx在内存中连续，但是Q不连续
            // right_idx = srcs[srcidx].margin_pts[left_idx].next_src_idx;     // right_idx的Q和left_idx的Q连续，但是内存不连续
            right_idx = CC[left_idx].next_idx;
            skip = CC[right_idx].skip;
            while(skip)
            {
                right_idx = CC[right_idx].next_idx;
                skip = CC[right_idx].skip;
            }

            // printf("idx real: %d, left_idx: %d, right_idx: %d\n",idx+BATCH_SIZE*batchidx,left_idx,right_idx);

            left_psi = CC[left_idx].psi;
            if(left_idx == 0)
            {
                quantile_left = 0;// 左点需要考虑 第0点和第1点之间的情况
            }
            else
            {
                quantile_left = Weight[left_idx-1];         // -1 是因为求和带来的指标错位
            }
            
            right_psi = CC[right_idx].psi;
            if(right_psi == 0)
            {
                right_psi = 2*PI;             // 右点需要将01认同，以满足rightQ总是大于leftQ
            }
            // quantile_right = ErrorWeight[left_idx -1 +1];       // 注意quantile时right和left的差别
            // 不创建变量（塞不下了），直接调用共享内存的数据放在缓存里

            // if(idx==0)
            // {
            //     for(int iii = 0;iii<32;iii++)
            //     {
            //         printf("%.12f,",Weight[iii]);
            //     }
            //     printf("\n");
            // }
            // 线性插值
            psi = left_psi + (right_psi - left_psi) * (quantile-quantile_left)/(Weight[left_idx -1 +1]-quantile_left);         // 不创建变量（塞不下了），quantile_right = ErrorWeight[left_idx -1 +1];

            CC[idx+BATCH_SIZE*batchidx].psi = psi;


            // 连接顺序：
            ParentIdx[idx] = left_idx; 
        }
        __syncthreads();

        if(idx<BATCH_SIZE)
        {
            // 前半截是“谁连自己”
            // 如果前一个点和自己是姊妹节点
            if((idx >= 1) && (ParentIdx[idx-1]==ParentIdx[idx]))      // 0号的前一个要取模，但负数取模有一定问题，而且反正0号线程肯定连parent
            {
                CC[idx-1 + batchidx*BATCH_SIZE].next_idx = idx + batchidx*BATCH_SIZE;     // 前一个点后面连自己
                // CC[idx + batchidx*BATCH_SIZE].prev_src_idx = idx-1 + batchidx*BATCH_SIZE;

            }
            else
            {
                CC[left_idx].next_idx = idx + batchidx*BATCH_SIZE;               // parent节点后面连自己
                // CC[idx + batchidx*BATCH_SIZE].prev_src_idx = left_idx;

            }
            // 后半截是“自己连谁”
            // 如果自己是姊妹节点中的老幺
            if(idx==(BATCH_SIZE-1))     // 如果是本batch中的最后一个
            {
                CC[idx + batchidx*BATCH_SIZE].next_idx = right_idx;                   // 自己连 right idx
                // CC[right_idx].prev_src_idx = idx + batchidx*BATCH_SIZE;

            }
            else{
                if(ParentIdx[idx] != ParentIdx[idx+1])  // 如果是普通情况
                {
                    CC[idx + batchidx*BATCH_SIZE].next_idx = right_idx;                   // 自己连 right idx
                    // CC[right_idx].prev_src_idx = idx + batchidx*BATCH_SIZE;

                }
                // else{}   // 前面连好了已经
            }            
        }


              
    }

}
template _Global void psiLoc(CC_t<float>* CC, int batchidx);
template _Global void psiLoc(CC_t<double>* CC, int batchidx);




// size<<< (1,1,NSRCS),(NCRIT)>>>
// shared size: int*64
// 暂不考虑x轴对称，考虑的话能减一半比较
template<isFloating f_T>
_Global void CollideTest(int* srclist, src_ext_t<f_T>* srcs, const CC_t<f_T>* CC, int batchidx)
{
    int srcidx = srclist[int(blockIdx.z)];

    // extern __shared__ int shared_mem[];
    __dyn_shared__(int, shared_mem);

    int* Ncross = (int*) shared_mem;
    // complex_t<f_T>* causLoc = (complex_t<f_T>*) &shared_mem[BATCH_SIZE];
    int idx = threadIdx.x;              // 0 ~ 1024
    int next_idx;
    bool skip;
    
    // complex_t<f_T> const * loc0;
    // complex_t<f_T> const * loc1;
    complex_t<f_T> zeta;
    // int temp;
    complex_t<f_T> posi0, posi1;


    if(!CC[idx].skip)
    {
        next_idx = CC[idx].next_idx;
        skip = CC[next_idx].skip;
        while(skip)
        {
            next_idx = CC[next_idx].next_idx;
            skip = CC[next_idx].skip;
        }

        for(int i=0;i<BATCH_SIZE;i++){Ncross[i]=0;}      // 0 ~ 64
        // for(int j=0;j<NBODY*NBODY;j++)
        // {
        //     causLoc[idx+j*NCRIT] = caus[idx+j*NCRIT];       // input to shared, save the reload of caus_loc_1
        //     // printf("(%.12f+%.12fj)\n",caus[idx+j*NCRIT].re,caus[idx+j*NCRIT].im);
        // }

        // printf("check %d\n",idx);
        for(int j=0;j<NBODY*NBODY;j++)    // Caustic 有4段
        {
            // loc0 = &(causLoc[idx + j*NCRIT]);
            // loc1 = &(causLoc[(idx+1)%NCRIT + next_j[idx + j*NCRIT]*NCRIT]);
            // posi0 = caus[idx+j*NCRIT];
            // posi1 = caus[(idx+1)%NCRIT + next_j[idx + j*NCRIT]*NCRIT];
            
            posi0 = CC[idx].caus[j];
            posi1 = CC[next_idx].caus[CC[idx].next_j[j]];


            // printf("idx:%d, j:%d, nextj:%d\n",idx,j,next_j[idx + j*NCRIT]);
            // printf("%d,%d\n",idx,Ncross[idx]);
            for(int i=0;i<BATCH_SIZE;i++)
            {
                zeta = srcs[srcidx].margin_pts[i + batchidx*BATCH_SIZE].position;
                zeta = complex_t<f_T>(zeta.re,fabs(zeta.im));                // 所有源点翻到x轴以上
                if((zeta.re<posi0.re && zeta.re>posi1.re)||(zeta.re>posi0.re && zeta.re<posi1.re)     )// 初步测试，寻找有没有两个点，其x坐标把 i 点夹在中间
                // 原则上需要判断==的时候（算作一半），但是概率原则上为0，就算了
                {
                    // 尽量减少计算
                    // 只看上半截
                    if(posi0.im>=-2e-12 && posi1.im>=-2e-12)        // 需要考虑一定的数值误差
                    {
                        
                        if(zeta.im<=posi0.im && zeta.im<=posi1.im)   // 如果比两个点都低，那在下方
                        {
                            // temp = atomicAdd(&Ncross[i],1);
                            atomicAdd(&Ncross[i],1);
                            // printf("0 %d,%d,%d,(%.12f+%.12fj),(%.12f+%.12fj),(%.12f+%.12fj),%d\n",i + batchidx*BATCH_SIZE,idx+j*NCRIT,(idx+1)%NCRIT + next_j[idx + j*NCRIT]*NCRIT,zeta.re,zeta.im,posi0.re,posi0.im,posi1.re,posi1.im,temp);
                        }
                        else{
                            if(zeta.im>posi0.im || zeta.im>posi1.im)    // 在两个点为对角线的矩形之中
                            {
                                
                                if(TriAreaD(posi0.re,posi0.im,posi1.re,posi1.im,zeta.re,zeta.im) >0)    // 012正面积，封闭曲线内部（因为NCRIT是右手的）
                                {
                                    if(posi0.re>posi1.re)       // 0 right, 1 left, area positive ==> 2 below
                                    {
                                        // temp = atomicAdd(&Ncross[i],1);
                                        atomicAdd(&Ncross[i],1);
                                        // printf("1 %d,%d,(%.12f+%.12fj),(%.12f+%.12fj),(%.12f+%.12fj),%d\n",i + batchidx*BATCH_SIZE,idx+j*NCRIT,zeta.re,zeta.im,posi0.re,posi0.im,posi1.re,posi1.im,temp);
                                    }
                                }
                                else
                                {
                                    if(posi0.re<posi1.re)       // 0 left, 1 right, area negative ==> 2 below
                                    {
                                        // temp = atomicAdd(&Ncross[i],1);
                                        atomicAdd(&Ncross[i],1);
                                        // printf("2 %d,%d,(%.12f+%.12fj),(%.12f+%.12fj),(%.12f+%.12fj),%d\n",i + batchidx*BATCH_SIZE,idx+j*NCRIT,zeta.re,zeta.im,posi0.re,posi0.im,posi1.re,posi1.im,temp);
                                    }                               
                                }
                            }
                            // else{}    // 比两个点都高，点在散焦线上方（不与向上的射线相交），不做操作
                            
                        }
                    }
                }
            }
        }
        // 全检测完了
        __syncthreads();
        if(idx<BATCH_SIZE)
        {
            // printf("%d,%d\n",idx + batchidx*BATCH_SIZE,Ncross[idx]);
            if((Ncross[idx] & 1)==1)    // mod 2 == 1, inside
            {
                srcs[srcidx].margin_pts[idx + batchidx*BATCH_SIZE].Nphys = 5;
            }
            else        // outside
            {
                srcs[srcidx].margin_pts[idx + batchidx*BATCH_SIZE].Nphys = 3;
            }
            srcs[srcidx].margin_pts[idx + batchidx*BATCH_SIZE].NcrossCaus = Ncross[idx];
        }        
    }


}

template _Global void CollideTest(int* srclist, src_ext_t<float>* srcs, const CC_t<float>* CC, int batchidx);
template _Global void CollideTest(int* srclist, src_ext_t<double>* srcs, const CC_t<double>* CC, int batchidx);



// size<<< (1,1,NSRCS),(NCRIT)>>>
// shared size: int* 4 * BATCH_SIZE (只看出buried的最大和最小，如果出现多个交替也不管)
// 暂不考虑x轴对称，考虑的话能减一半比较
// 一个线程对应一个焦散线上的区间，考虑任意方向
template<isFloating f_T>
_Global void CollideTest2(int* srclist, src_ext_t<f_T>* srcs, const CC_t<f_T>* CC, int batchidx)
{
    int srcidx = srclist[int(blockIdx.z)];

    // extern __shared__ int shared_mem[];
    __dyn_shared__(int, shared_mem);
    int* Ncross = (int*) shared_mem;
    int* NcrossNext = (int*) &shared_mem[BATCH_SIZE];
    int* loc_min = (int*) &shared_mem[2*BATCH_SIZE];
    int* loc_max = (int*) &shared_mem[3*BATCH_SIZE];
    // complex_t<f_T>* causLoc = (complex_t<f_T>*) &shared_mem[BATCH_SIZE];
    int idx = threadIdx.x;              // 0 ~ 1024
    int next_idx;
    
    // complex_t<f_T> const * loc0;
    // complex_t<f_T> const * loc1;
    complex_t<f_T> zeta0,zeta1;
    int temp,temp_idx;
    complex_t<f_T> posi0, posi1;
    // complex_t<float> new0,new1;
    complex_t<f_T> new0,new1;
    bool skip;

    f_T& tempf = zeta1.re;

    f_T right_idx_Q, left_idx_Q;

    if(idx<BATCH_SIZE)
    {
        Ncross[idx] = 0;
        NcrossNext[idx] = 0;
        loc_min[idx] = 100000000;
        loc_max[idx] = 0;
    }
    if(idx==0)
    {
        srcs[srcidx].NBuried = 0;
    }
    __syncthreads();

    if(!CC[idx].skip)
    {

        next_idx = CC[idx].next_idx;
        skip = CC[next_idx].skip;
        while(skip)
        {
            next_idx = CC[next_idx].next_idx;
            skip = CC[next_idx].skip;
        }

        // for(int j=0;j<NBODY*NBODY;j++)
        // {
        //     causLoc[idx+j*NCRIT] = caus[idx+j*NCRIT];       // input to shared, save the reload of caus_loc_1
        //     // printf("(%.12f+%.12fj)\n",caus[idx+j*NCRIT].re,caus[idx+j*NCRIT].im);
        // }

        // printf("check %d\n",idx);
        for(int j=0;j<NBODY*NBODY;j++)    // Caustic 有4段
        {
            // loc0 = &(causLoc[idx + j*NCRIT]);
            // loc1 = &(causLoc[(idx+1)%NCRIT + next_j[idx + j*NCRIT]*NCRIT]);
            // posi0 = caus[idx+j*NCRIT];
            // posi1 = caus[(idx+1)%NCRIT + next_j[idx + j*NCRIT]*NCRIT];
            posi0 = CC[idx].caus[j];
            posi1 = CC[next_idx].caus[CC[idx].next_j[j]];

            // printf("idx:%d, j:%d, nextj:%d\n",idx,j,next_j[idx + j*NCRIT]);
            // printf("%d,%d\n",idx,Ncross[idx]);
            for(int i=0;i<BATCH_SIZE;i++)
            {
                zeta0 = srcs[srcidx].margin_pts[i + batchidx*BATCH_SIZE].position;
                temp_idx = srcs[srcidx].margin_pts[i + batchidx*BATCH_SIZE].next_src_idx;
                skip = srcs[srcidx].margin_pts[temp_idx].skip;
                while(skip==true)
                {
                    temp_idx = srcs[srcidx].margin_pts[temp_idx].next_src_idx;
                    skip = srcs[srcidx].margin_pts[temp_idx].skip;
                }
                zeta1 = srcs[srcidx].margin_pts[temp_idx].position;
                // Delta_zeta = srcs[srcidx].margin_pts[temp].position - srcs[srcidx].margin_pts[i + batchidx*BATCH_SIZE].position;

                // 将zeta0归位原点，zeta1归为（1，0）
                zeta1 = zeta1 - zeta0;

                // new0 = complex_t<float>(posi0-zeta0) / complex_t<float>(zeta1);
                // new1 = complex_t<float>(posi1-zeta0) / complex_t<float>(zeta1);
                new0 = (posi0-zeta0) / (zeta1);         // 归一化的焦散线点
                new1 = (posi1-zeta0) / (zeta1);

                if((new0.im>=0) != (new1.im>0)) // 如果转动参考系之后，锚点的y坐标异号
                {
                    // 如果两个参考点都比1大
                    if((new0.re > 1.)&&(new1.re > 1.))
                    {
                        temp = atomicAdd(&Ncross[i],1);
                        temp = atomicAdd(&NcrossNext[i],1);
                    }
                    else
                    {
                        // 如果参考点既不是特别大也不是特别小
                        if( !( (new0.re < 0.) && (new1.re < 0.) ) )
                        {
                            // tempf is temporary valuable
                            tempf = (new1.im*new0.re - new0.im*new1.re) / (new1.im - new0.im);      // 焦散线连线与x轴点交点
                            if(tempf > 1.)
                            {
                                temp = atomicAdd(&Ncross[i],1);
                                temp = atomicAdd(&NcrossNext[i],1);
                            }
                            else{
                                if(tempf > 0.)
                                {
                                    temp = atomicAdd(&Ncross[i],1);
                                    temp = atomicMin(&loc_min[i],int(tempf*1e8));
                                    temp = atomicMax(&loc_max[i],int(tempf*1e8));
                                }
                                // else{}
                            }
                        }
                        // 如果两个参考点都比0小，什么都不做（+0）
                        // else{}
                        
                    }
                }

            }
        }
    }
    // 全检测完了
    __syncthreads();

    // 前面idx用于 crit， 后面idx用于新的 points idx
    if(idx<BATCH_SIZE)
    {
        // printf("%d,%d\n",idx + batchidx*BATCH_SIZE,Ncross[idx]);
        if((Ncross[idx] & 1)==1)    // % 2 == 1， inside caustic
        {
            srcs[srcidx].margin_pts[idx + batchidx*BATCH_SIZE].Nphys = 5;
        }
        else        // outside
        {
            srcs[srcidx].margin_pts[idx + batchidx*BATCH_SIZE].Nphys = 3;
        }
        srcs[srcidx].margin_pts[idx + batchidx*BATCH_SIZE].NcrossCaus = Ncross[idx];



        if(abs(Ncross[idx]-NcrossNext[idx]) > 1)
        {
            temp = atomicAdd(&(srcs[srcidx].NBuried),1);
            srcs[srcidx].idx_buried[temp] = idx + batchidx*BATCH_SIZE;
            temp_idx = srcs[srcidx].margin_pts[idx + batchidx*BATCH_SIZE].next_src_idx;
            skip = srcs[srcidx].margin_pts[temp_idx].skip;
            while(skip==true)
            {
                temp_idx = srcs[srcidx].margin_pts[temp_idx].next_src_idx;
                skip = srcs[srcidx].margin_pts[temp_idx].skip;
            }

            left_idx_Q = srcs[srcidx].margin_pts[idx + batchidx*BATCH_SIZE].Q;
            right_idx_Q = srcs[srcidx].margin_pts[temp_idx].Q;
            if(right_idx_Q < left_idx_Q)
            {
                right_idx_Q += 1.;
            }
            tempf = left_idx_Q + (loc_min[idx]) / f_T(1e8) * (right_idx_Q - left_idx_Q);
            srcs[srcidx].Qmin_buried[temp] = tempf - int(tempf);
            tempf = left_idx_Q + (loc_max[idx]) / f_T(1e8) * (right_idx_Q - left_idx_Q);
            srcs[srcidx].Qmax_buried[temp] = tempf - int(tempf);
            // printf("(%d,%d), buried caustic, Qmin: %.12f, Qmax: %.12f\n",srcidx,idx + batchidx*BATCH_SIZE,srcs[srcidx].Qmin[temp], srcs[srcidx].Qmax[temp]);
        }
    }

}
template _Global void CollideTest2(int* srclist, src_ext_t<float>* srcs, const CC_t<float>* CC, int batchidx);
template _Global void CollideTest2(int* srclist, src_ext_t<double>* srcs, const CC_t<double>* CC, int batchidx);
