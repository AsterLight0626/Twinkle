#include "GlobalFuncs_v2.h"

// size: <<<(NSRCS,1,1),BATCH_SIZE>>>
// int:2 f_T:1
// float: 12, double:16
template<isFloating f_T>
_Global void MarginBatch(int const * srclist, src_ext_t<f_T>* srcs)
{
    const int srcidx = srclist[int(blockIdx.x)];
    int idx = threadIdx.x;
    f_T Q;

    Q = (double((idx))) / double(BATCH_SIZE);
    srcs[srcidx].margin_pts[idx].position = Margin_Q(srcs[srcidx].shape,Q);
    srcs[srcidx].margin_pts[idx].next_src_idx = ((idx+1) & (BATCH_SIZE-1));           // BATCH_SIZE = 2^n, then (idx+1)%BATCH_SIZE == (idx+1) & (BATCH_SIZE-1)
    srcs[srcidx].margin_pts[((idx+1) & (BATCH_SIZE-1))].prev_src_idx = idx;
    srcs[srcidx].margin_pts[idx].changed = 1;
    srcs[srcidx].margin_pts[idx].Q = Q;
    srcs[srcidx].margin_pts[idx].special_Err = 0.;
    srcs[srcidx].margin_pts[idx].Nphys = -1;
    for(int j=0;j<LENGTH-1;j++)
    {
        srcs[srcidx].margin_pts[idx].deltaS[j]=0;
        srcs[srcidx].margin_pts[idx].deltaS_Err[j]=0;
        srcs[srcidx].margin_pts[idx].deltaS_new[j]=0;
        srcs[srcidx].margin_pts[idx].Err_new[j]=0;
    }

}
template _Global void MarginBatch(int const * srclist, src_ext_t<float>* srcs);
template _Global void MarginBatch(int const * srclist, src_ext_t<double>* srcs);



// size: <<< <<<(NSRCS,1,1),BATCH_SIZE>>>
// int: 1, f_T: 14, c_f_T: 7, bool:4, img_f_T:5
// float: 170, double 322
template<isFloating f_T>
_Global void SolveBatch(int const * srclist, src_ext_t<f_T>*srcs,const src_params_t<f_T>* src_params, int batchidx, bool muted, bool check)
{
    const int srcidx = srclist[int(blockIdx.x)];
    // int idx = batchidx*BATCH_SIZE + threadIdx.x;          // margin point srcs idx
    src_pt_t<f_T>* pt_here;
    complex_t<f_T> loc_src;    
    f_T m1,m2,s;	
    f_T moving;
    complex_t<f_T> coeffs[LENGTH];
    image_pt_t<f_T> temp_images[LENGTH-1];
    bool fail;

    complex_t<f_T>& temp_z_s = coeffs[4];
    complex_t<f_T>& zc = coeffs[3];
    complex_t<f_T>* items = coeffs;
    f_T temp_jacobi[LENGTH-1];

    complex_t<f_T>& pzeta_pconjz = coeffs[2];
    complex_t<f_T>& ppzeta_ppconjz = coeffs[3];
	complex_t<f_T>& dzeta = items[0];
	complex_t<f_T>& dz = items[1];

    bool fail_phys;
    f_T errors[LENGTH-1];
    bool skip = false;
    bool fail_pp=0;
    bool ambiguous = false;





    pt_here = &(srcs[srcidx].margin_pts[int(batchidx*BATCH_SIZE + threadIdx.x)]);

    loc_src = pt_here->position;

    // 二体特例
    s = src_params->s;
    m2 = src_params->q;       // temprately q
	m1 = 1 / (1 + m2);
	m2 *= m1;

    moving = PolyCoeffTwinkle(s,m1,m2, loc_src, coeffs);


    fail = cmplx_roots_gen(temp_images,coeffs,LENGTH-1,true,false);
    // 先存到局域位置里（反正images不大，float版本只有12Byte
    if(! fail)
    {
        for(int j=0;j<LENGTH-1;j++)
        {
            temp_z_s = temp_images[j].position;
            temp_z_s.re += s;
            temp_z_s.im = -temp_z_s.im;
            zc = conj(temp_images[j].position);

            items[0] = m1 / (temp_z_s * temp_z_s);
            items[1] = m2 / (zc * zc);
            pzeta_pconjz = items[0] + items[1];

            temp_jacobi[j] = 1- norm(pzeta_pconjz);
            if(fabsf(temp_jacobi[j])<1e-5f){ambiguous = true;}
            items[0] = -f_T(2) * items[0] / temp_z_s;
            items[1] = -f_T(2) * items[1] / zc;
            ppzeta_ppconjz = items[0] + items[1];

            dzeta = complex_t<f_T>(0.f,1.f) * (loc_src - srcs[srcidx].shape.loc_centre);
            
            dz = (dzeta - pzeta_pconjz*conj(dzeta)) / temp_jacobi[j];		// equation20 in VBB
            pt_here->dz[j] = dz;

            pt_here->wedge[j] = (srcs[srcidx].shape.rho*srcs[srcidx].shape.rho + imag(dz*dz * dzeta * conj(ppzeta_ppconjz))) / temp_jacobi[j];
            

            temp_images[j].parity = (temp_jacobi[j]>0) ? 0 : 1 ;
            temp_images[j].position.re += moving;
        }
        // 判断实虚
        fail_phys = isPhysicalAll(temp_images,s,m1,m2,loc_src,&(pt_here->Nphys),errors);

        if(fail_phys)
        {
            skip = true;
            for(int j=0;j<LENGTH-1;j++)
            {
                pt_here->images[j] = temp_images[j];
                #if DetailedRecord
                    pt_here->jacobi[j] = temp_jacobi[j];
                    pt_here->phys_err[j] = errors[j];
                #endif
            }
        }
        else
        {
            skip = false;
            // 在这里整合一下checkpp
            if(check)
            {
                if(pt_here->Nphys==5 && ambiguous)
                // if(false)
                {
                    fail_pp = PolishPP2(temp_images,temp_jacobi,loc_src);
                }
                else
                {
                    fail_pp = PolishPP(temp_images,temp_jacobi);
                }
            }    // 如果不检测，fail_pp在开头设置为0


            if(! fail_pp){
                for(int j=0;j<LENGTH-1;j++)
                {
                    pt_here->images[j] = temp_images[j];
                    #if DetailedRecord
                    // {
                        pt_here->jacobi[j] = temp_jacobi[j];
                        pt_here->phys_err[j] = errors[j];
                    // }
                    #endif
                }
            }
            else{
                skip = true;              
            }                     
        }
    }
    else
    {
        skip = true;          
    }

    pt_here->skip = skip;
    if(skip)
    {
        for(int j=0;j<LENGTH-1;j++)
        {
            pt_here->deltaS[j] = 0;
            pt_here->deltaS_Err[j] = 0;
            pt_here->deltaS_new[j] = 0;
            pt_here->Err_new[j] = 0;            
            pt_here->next_idx[j] = -2;
            pt_here->next_j[j] = -2;
            pt_here->wedge[j] = 0;
            pt_here->dz[j] = 0;
        }     
        pt_here->error_interval = 0;
        pt_here->area_interval = 0;
        pt_here->changed = 0;
    }

}
template _Global void SolveBatch(int const * srclist,src_ext_t<float>*srcs,const src_params_t<float>* src_params, int batchidx, bool muted, bool check);
template _Global void SolveBatch(int const * srclist,src_ext_t<double>*srcs,const src_params_t<double>* src_params, int batchidx, bool muted, bool check);


// 修改gridsize，只对新修改过的点计算
// size: <<<(NSRCS,1, 1),BATCH_SIZE>>>
// int: 32, c_f_T: 13, bool: 2
// float: 234, double: 338
template<isFloating f_T>
_Global void FindNextBatch3(int* srclist, src_ext_t<f_T>* srcs, int batchidx, bool muted)
{
    const int srcidx = srclist[int(blockIdx.x)];

    // cross check initialize, will be used later
    if((threadIdx.x==0))
    {
        srcs[srcidx].Ncross=0;
    }            

    
    __syncthreads();


    int prev_idx;    
    int idx;
    int next_idx;
    // bool skip;
    int idx_list[3];

    complex_t<f_T> posi_prev[3];   // only odd parity will connect to previous source points
    int j_prev[3];          // 存储prev的j坐标
    complex_t<f_T> posi_here_odd[3];
    int j_here_odd[3];
    complex_t<f_T> posi_here_even[2];
    int j_here_even[2];
    complex_t<f_T> posi_next[2];   // only even parity will connect to next source points
    int j_next[2];          // 存储next的j坐标

    int Nphys_prev;
    int Nphys_here;
    int Nphys_next;

    int temp0;
    int temp1;

    int j_even_candidate;
    int j_odd_candidate;     // 即将连接自己的像的j坐标,-1表示没有

    int order_e[2];     // positive parity
    int order_o[3];     // negative parity

    complex_t<f_T> posi_next_odd[2];    // Nphys==3
    complex_t<f_T> posi_prev_even[1];   // Nphys==3
    int order_c[3];     // c means cross caustic

    src_pt_t<f_T>* pt_here;



    idx_list[1] = threadIdx.x + batchidx*BATCH_SIZE;
    if(!(srcs[srcidx].margin_pts[idx_list[1]].skip))
    // if(true)
    {

        idx_list[0] = PrevSrcIdx(srcs[srcidx].margin_pts, idx_list[1]);
        idx_list[2] = NextSrcIdx(srcs[srcidx].margin_pts, idx_list[1]);

        for(int loop3=0;loop3<3;loop3++)
        {
            idx = idx_list[loop3];
            pt_here = &(srcs[srcidx].margin_pts[idx]);

            // if(idx<0){temp0=3;}     // 如果idx小于0，说明不用算，temp0=3认同为已经计算过了
            if(loop3!=1 && idx>=batchidx*BATCH_SIZE){temp0=3;}
            else{
                temp0 = atomicExch(&(pt_here->changed),3);      // temp0 可能为0（表示该点信息没更新（不应该存在））， 1（需要计算），3（已经计算过了）
            }

            if(temp0==1)
            {
                if(pt_here->skip==true)
                {
                    temp0 = 3;  // 如果skip了，也算完成计算
                }
            }

            if(temp0==1)   // temp0==1表示该点修改过，需要计算；skip==false说明这一点的信息是真实的
            {
                prev_idx = (loop3<1) ? PrevSrcIdx(srcs[srcidx].margin_pts, idx) : idx_list[loop3-1];
                next_idx = (loop3>1) ? NextSrcIdx(srcs[srcidx].margin_pts, idx) : idx_list[loop3+1];

                // 获取 prev_idx, idx, next_idx, 与之前一样

                // 最后一位有时是未初始化的，保险起见初始化一下
                posi_prev[2] = complex_t<f_T>(0.,0.);
                j_prev[2] = -1;
                posi_here_odd[2] = complex_t<f_T>(0.,0.);
                j_here_odd[2] = -1;
                posi_here_even[1] = complex_t<f_T>(0.,0.);
                j_here_even[1] = -1;
                posi_next[1] = complex_t<f_T>(0.,0.);
                j_next[1] = -1;



                temp0=0;
                temp1=0;

                j_even_candidate = -1;
                j_odd_candidate = -1;     // 即将连接自己的像的j坐标,-1表示没有

                Nphys_prev = srcs[srcidx].margin_pts[prev_idx].Nphys;
                Nphys_here = pt_here->Nphys;
                Nphys_next = srcs[srcidx].margin_pts[next_idx].Nphys;

                for(int j=0;j<5;j++)
                {
                    // imgs_here[j] = pt_here->images[j];
                    if(pt_here->images[j].physical)
                    {
                        if(pt_here->images[j].parity==0)
                        {
                            posi_here_even[temp0] = pt_here->images[j].position;
                            j_here_even[temp0] = j;
                            temp0++;
                        }
                        else{
                            posi_here_odd[temp1] = pt_here->images[j].position;
                            j_here_odd[temp1] = j;
                            temp1++;   
                        }
                    }
                    else    // 当前点是虚的
                    {
                        // output
                        pt_here->next_idx[j] = -1;
                        pt_here->next_j[j] = -1;
                    }
                }
                temp0=0;        // 省一点内存
                temp1=0;
                for(int j=0;j<5;j++)
                {              
                    if((loop3==1)||(prev_idx>=batchidx*BATCH_SIZE))
                    {
                        if((srcs[srcidx].margin_pts[prev_idx].images[j].physical) && (srcs[srcidx].margin_pts[prev_idx].images[j].parity==1))
                        {
                            posi_prev[temp1] = srcs[srcidx].margin_pts[prev_idx].images[j].position;
                            j_prev[temp1] = j;
                            temp1++;
                        }                        
                    }
                    // else: do nothing
                    if((loop3==1)||(next_idx>=batchidx*BATCH_SIZE))
                    {
                        if((srcs[srcidx].margin_pts[next_idx].images[j].physical) && (srcs[srcidx].margin_pts[next_idx].images[j].parity==0))
                        {
                            posi_next[temp0] = srcs[srcidx].margin_pts[next_idx].images[j].position;
                            j_next[temp0] = j;
                            temp0++;
                        }                        
                    }
                    // else: do nothing
                }

                // 如果Nphys是3，imgs最后一位是未初始化的（全都初始化了，231227）

                // input images finished

                // 在1.0版本中，接下来是5对5连接，找最小的总长度，而且还要连两遍（prev 2 here, here 2 next)
                // 现在不需要了
                // 事实上FindNext开销总是很小的，主要是为了健壮性

                // 按手征分类
                // 偶手征只看下，奇手征只看上

                // even parity, 1c1 or 2c2
                if((loop3==1)||(next_idx>=batchidx*BATCH_SIZE))
                {
                    if(Nphys_here==Nphys_next)
                    {
                        if(Nphys_here==3)
                        {
                            // output
                            pt_here->next_idx[j_here_even[0]] = next_idx;
                            pt_here->next_j[j_here_even[0]] = j_next[0];
                        }
                        else    // Nphys_here == 5, 两个比较大小
                        {
                            // int order_e[2];
                            ConnectOrder(posi_here_even,posi_next,2,2,order_e);
                            // output
                            for(int ii=0;ii<2;ii++)
                            {
                                pt_here->next_idx[j_here_even[ii]] = next_idx;
                                pt_here->next_j[j_here_even[ii]] = j_next[order_e[ii]];
                            }
                        }
                    }
                    else    // Nphys_here!=Nphys_next
                    {
                        if(Nphys_here > Nphys_next)     // here: 2, next: 1
                        {
                            if(norm(posi_here_even[0]-posi_next[0]) < norm(posi_here_even[1]-posi_next[0])) // 连近的，现在00是近的
                            {
                                // output
                                pt_here->next_idx[j_here_even[0]] = next_idx;
                                pt_here->next_j[j_here_even[0]] = j_next[0];
                                pt_here->next_idx[j_here_even[1]] = idx;        // connect to here
                                j_even_candidate = j_here_even[1];       // 1 还没连
                            }
                            else    // 10是近的
                            {
                                // output
                                pt_here->next_idx[j_here_even[1]] = next_idx;
                                pt_here->next_j[j_here_even[1]] = j_next[0];
                                pt_here->next_idx[j_here_even[0]] = idx;        // connect to here
                                j_even_candidate = j_here_even[0];       // 0 还没连
                            }
                        }
                        else            // here: 1, next: 2
                        {
                            if(norm(posi_here_even[0]-posi_next[0]) < norm(posi_here_even[0]-posi_next[1])) // 连近的，现在00是近的
                            {
                                // output
                                pt_here->next_idx[j_here_even[0]] = next_idx;
                                pt_here->next_j[j_here_even[0]] = j_next[0];
                            }            
                            else    // 01是近的
                            {
                                // output
                                pt_here->next_idx[j_here_even[0]] = next_idx;
                                pt_here->next_j[j_here_even[0]] = j_next[1];   
                            }            
                        }
                    }
                }


                // odd parity, 2c2 or 3c3
                if((loop3==1)||(prev_idx>=batchidx*BATCH_SIZE))
                {
                    if(Nphys_here==Nphys_prev)
                    {
                        if(Nphys_here==3)
                        {
                            ConnectOrder(posi_here_odd, posi_prev,2,2,order_o);
                            // output
                            for(int ii=0;ii<2;ii++)
                            {
                                pt_here->next_idx[j_here_odd[ii]] = prev_idx;
                                pt_here->next_j[j_here_odd[ii]] = j_prev[order_o[ii]];
                            }
                        }
                        else
                        {
                            ConnectOrder(posi_here_odd, posi_prev,3,3,order_o);
                            for(int ii=0;ii<3;ii++)
                            {
                                pt_here->next_idx[j_here_odd[ii]] = prev_idx;
                                pt_here->next_j[j_here_odd[ii]] = j_prev[order_o[ii]];
                            }                        
                        }
                    }
                    else
                    {
                        if(Nphys_here > Nphys_prev) // here: 3, prev: 2 
                        {
                            ConnectOrder(posi_here_odd,posi_prev,3,2,order_o);      // order is 0,1 or -1 (cross candidate)
                            // output
                            for(int ii=0;ii<3;ii++)
                            {
                                if(order_o[ii]>=0)  // real
                                {
                                    pt_here->next_idx[j_here_odd[ii]] = prev_idx;
                                    pt_here->next_j[j_here_odd[ii]] = j_prev[order_o[ii]];                                
                                }
                                else    // cross candidate
                                {
                                    pt_here->next_idx[j_here_odd[ii]] = idx;    // connect to here
                                    j_odd_candidate = j_here_odd[ii];
                                }
                            }                          
                        }
                        else        // here: 2, prev: 3
                        {

                            ConnectOrder(posi_here_odd,posi_prev,2,3,order_o);

                            // output
                            for(int ii=0;ii<2;ii++)
                            {
                                pt_here->next_idx[j_here_odd[ii]] = prev_idx;
                                pt_here->next_j[j_here_odd[ii]] = j_prev[order_o[ii]];
                            }                        
                        }
                    }
                }

                // 检查跨焦散
                // even parity connect to here
                if(j_even_candidate>=0 && j_even_candidate<LENGTH-1 )
                {
                    // 记录跨焦散位置
                    temp1 = atomicAdd(&srcs[srcidx].Ncross,1);
                    if(temp1>NCROSS_CANDIDATE)
                    {
                        printf("Warning: too many cross happen, NCROSS_CANDIDATE is not adequent. srcidx: %d, batchidx: %d\n",srcidx,batchidx);
                    }
                    srcs[srcidx].idx_cross[temp1] = idx;        // positive means next N=3
                    srcs[srcidx].j_cross[temp1] = j_even_candidate;
                    temp1=0;
                    for(int j=0;j<5;j++)
                    {
                        if((srcs[srcidx].margin_pts[next_idx].images[j].physical)&&(srcs[srcidx].margin_pts[next_idx].images[j].parity==1))
                        {
                            // 由Solve保证，当Nphys==3时，像一定是2个odd，1个even
                            posi_next_odd[temp1] = srcs[srcidx].margin_pts[next_idx].images[j].position;
                            // j_next_odd[temp1] = j;
                            temp1++;
                        }
                    }
                    // posi_here_odd[3] connect to posi_next_odd[2]
                    // int order_c[3];     // c means cross caustic
                    ConnectOrder(posi_here_odd,posi_next_odd,3,2,order_c);
                    for(int ii=0;ii<3;ii++)
                    {
                        if(order_c[ii]==-1)     // 连接会有一个落单的，即离每个nextodd都很远，这个点被even cross连接
                        {
                            // output j
                            pt_here->next_j[j_even_candidate] = j_here_odd[ii];
                        }
                    }


                }

                if(j_odd_candidate>=0 && j_odd_candidate<LENGTH-1)
                {
                    // 记录跨焦散位置
                    temp0 = atomicAdd(&srcs[srcidx].Ncross,1);

                    if(temp0>NCROSS_CANDIDATE)
                    {
                        printf("Warning: too many cross happen, NCROSS_CANDIDATE is not adequent. srcidx: %d, batchidx: %d\n",srcidx,batchidx);
                    }
                    srcs[srcidx].idx_cross[temp0] = -idx-1;       // negative means previous N=3, -1 to avoid case of idx==0
                    srcs[srcidx].j_cross[temp0] = j_odd_candidate;


                    for(int j=0;j<5;j++)
                    {
                        if((srcs[srcidx].margin_pts[prev_idx].images[j].physical)&&(srcs[srcidx].margin_pts[prev_idx].images[j].parity==0))
                        {
                            // 由Solve保证，当Nphys==3时，像一定是2个odd，1个even
                            posi_prev_even[0] = srcs[srcidx].margin_pts[prev_idx].images[j].position;
                        }
                    }
                    // posi_here_even[2] connect to posi_prev_even[1]
                    if(norm(posi_here_even[0]-posi_prev_even[0]) < norm(posi_here_even[1]-posi_prev_even[0]))
                    {
                        // posi_here_even[0] 更近，posi_here_even[1]落单
                        pt_here->next_j[j_odd_candidate] = j_here_even[1];
                    }
                    else    // posi_here_even[0]落单
                    {
                        pt_here->next_j[j_odd_candidate] = j_here_even[0];
                    }


                }


            }
        }
    }
    // else: 当前点skip，无运算
    else
    {
        idx = threadIdx.x + batchidx*BATCH_SIZE;
        prev_idx = PrevSrcIdx(srcs[srcidx].margin_pts, idx);
        next_idx = NextSrcIdx(srcs[srcidx].margin_pts, idx);
        if(prev_idx<batchidx*BATCH_SIZE && next_idx<batchidx*BATCH_SIZE)        // 中间 skip 了，前后都在旧点里
        {
            srcs[srcidx].margin_pts[prev_idx].changed = 0;
            srcs[srcidx].margin_pts[next_idx].changed = 0;
        }
    }            

}
template _Global void FindNextBatch3(int* srclist, src_ext_t<float>* srcs, int batchidx, bool muted);
template _Global void FindNextBatch3(int* srclist, src_ext_t<double>* srcs, int batchidx, bool muted);



// size: <<<(NSRCS,1,1),BATCH_SIZE>>>
// int: 7, f_T: 10, f_T*:1, c_f_T:12
// float:168, double: 304
// 只算普通点，cross不算
template <isFloating f_T>
_Global void AreaErrBatchSrc4(int* srclist, src_ext_t<f_T>* srcs, int batchidx, bool muted)
{
    const int srcidx = srclist[int(blockIdx.x)];

    int idx_list[3];
    int idx;
    complex_t<f_T> zs[2];
    complex_t<f_T> dzs[2];
    f_T ddelta_theta[2];
    f_T delta_theta;

    int next_i;
    int next_j;

    f_T wedge[2];       // x' ^ x''
    f_T E_1;
    f_T E_2;
    f_T E_3;            // 与VBB一致
    f_T E_4;            // difference of deltaS_p1 and deltaS_p2
    f_T deltaS_t;       // Trapezium approximation
    f_T deltaS_p;       // Parabolic correction
    // f_T deltaS_p1,deltaS_p2;
    f_T& deltaS_p1 = E_2;       // save some memory
    f_T& deltaS_p2 = E_3;       // save some memory
    f_T& deter = E_1;           // save some memory
    f_T Qnext,Qhere;
    f_T dtheta2,dtheta3;

    f_T Qnext_src,Qprev_src;
    f_T theta2[2],theta3[2];
    int next_src_idx;
    int prev_src_idx;

    int& temp0 = next_i;
    bool parity;

    // bool skip;

    src_pt_t<f_T>* pt_here;


    idx_list[1] = threadIdx.x + batchidx*BATCH_SIZE;
    if(!(srcs[srcidx].margin_pts[idx_list[1]].skip))
    {
        idx_list[0] = PrevSrcIdx(srcs[srcidx].margin_pts, idx_list[1]);
        idx_list[2] = NextSrcIdx(srcs[srcidx].margin_pts, idx_list[1]);         // 如果前后的点在最新的batch里，那么它们会在idx_list[1]里面被计算    

        for(int loop3=0;loop3<3;loop3++)
        {
            idx = idx_list[loop3];
            pt_here = &srcs[srcidx].margin_pts[idx];

            if(loop3!=1 && idx>=batchidx*BATCH_SIZE){temp0=0;}     // 如果idx小于0，说明不用算，temp0=0认同为已经计算过了
            else{
                temp0 = atomicExch(&(pt_here->changed),4);      // temp0 可能为0（已经算完了）， 1（需要连接（不应该出现）），3（待计算面积）
            }
            if(temp0==1)
            {
                if(!muted)
                {
                    printf("(srcidx: %d, idx: %d, batchidx: %d) Warning: not connect\n",srcidx,idx,batchidx);
                }
                srcs[srcidx].Break=1;
                srcs[srcidx].SolveSucceed = 1;          // 用于跳过后续计算，会在最后一个SumArea设置为fail
                #if DetailedRecord
                    srcs[srcidx].points_used = (batchidx+1) * BATCH_SIZE;  
                #endif
                continue;
            }

            if(temp0==3){

                Qhere = pt_here->Q;
                next_src_idx = (loop3>1) ? NextSrcIdx(srcs[srcidx].margin_pts,idx) : idx_list[loop3+1];
                prev_src_idx = (loop3<1) ? PrevSrcIdx(srcs[srcidx].margin_pts,idx) : idx_list[loop3-1];
                Qnext_src = srcs[srcidx].margin_pts[next_src_idx].Q;
                Qprev_src = srcs[srcidx].margin_pts[prev_src_idx].Q;




                ddelta_theta[0] = (Qnext_src - Qhere) * 2 * PI;      // parity==0, 正手征，连nextsrcidx
                ddelta_theta[1] = (Qprev_src - Qhere) * 2 * PI;      // parity==1, 负手征，连prevsrcidx            
                if(idx==0)      // srcs[srcidx].margin_pts[0].Q  must be equal to 0, and Q = 0 must be idx = 0 
                {
                    ddelta_theta[1] = (Qprev_src - 1.) * 2 * PI;
                }
                if(next_src_idx==0)
                {
                    ddelta_theta[0] = (1. - Qhere) * 2 * PI;                              
                }
                for(int i=0;i<2;i++)
                {
                    theta2[i] = ddelta_theta[i] * ddelta_theta[i];
                    theta3[i] = ddelta_theta[i] * theta2[i];                
                }
            
                // 对于每条链

                for(int j=0;j<LENGTH-1;j++)
                {
                    if(pt_here->images[j].physical)
                    {
                        next_i = pt_here->next_idx[j];
                        if(loop3!=1 && (next_i<batchidx*BATCH_SIZE))
                        {
                            continue;
                        }
                        parity = pt_here->images[j].parity;
                        zs[0] = pt_here->images[j].position;
                        
                        next_j = pt_here->next_j[j];
                        Qnext = srcs[srcidx].margin_pts[next_i].Q;

                        if(!((next_j>=0)&&(next_j<5)))
                        {
                            if(!muted)
                            {
                                printf("invalid next_j in srcidx: %d, idx: %d, next_j: %d\n",srcidx,idx,next_j);
                            }
                            srcs[srcidx].Break=1;
                            srcs[srcidx].SolveSucceed = 1;          // 用于跳过后续计算，会在最后一个SumArea设置为fail
                            #if DetailedRecord
                                srcs[srcidx].points_used = (batchidx+1) * BATCH_SIZE;  
                            #endif
                            continue;
                        }

                        zs[1] = srcs[srcidx].margin_pts[next_i].images[next_j].position;


                        wedge[0] = pt_here->wedge[j];
                        wedge[1] = srcs[srcidx].margin_pts[next_i].wedge[next_j];

                        dzs[0] = pt_here->dz[j];
                        dzs[1] = srcs[srcidx].margin_pts[next_i].dz[next_j];

                        

                        if(idx!=next_i)     // usually
                        {
                            deltaS_t = DeltaS1(zs[0],zs[1]);
                            // deltaS_t = (zs[0].re*zs[1].im - zs[0].im*zs[1].re)/2;
                            deltaS_p1 = (wedge[0] + wedge[1]) * theta3[int(parity)] /24.;
                            deltaS_p2 = wedge_product((zs[1]-zs[0]) , (dzs[1]-dzs[0])) * ddelta_theta[int(parity)] / 12.;
                            // if(srcidx==0 && idx==12)
                            // {
                            //     printf("j = %d, deltaS_p: %.16f\n",j,deltaS_p);
                            //     printf("j = %d, deltaS_p2: %.16f\n",j,deltaS_p2);
                            // }
                            deltaS_p = (deltaS_p1 + deltaS_p2) / 2.;
                            E_4 = fabs((deltaS_p1 - deltaS_p2));
                            deter = theta2[int(parity)] * abs_c(dzs[0]*dzs[1]);
                            if(fabs(deter)>2e-12){E_2 = 1.5 * fabs( deltaS_p *  (norm(zs[0]-zs[1]) / deter -1));}
                            else
                            {
                                E_2=0;
                            }        // usually caused by numericial error
                            E_1 = fabs((wedge[0] - wedge[1]) * theta3[int(parity)]) /48.;
                            // mode = 0;
                            E_3 = 0.1 * fabs(deltaS_p) * theta2[int(parity)];
                        }
                        else                // cross caustic
                        {
                            deltaS_t = 0;
                            deltaS_p = 0;
                            E_1 = 0;
                            E_2 = 0;
                            E_3 = 0;
                            E_4 = 0;
                        }

                        deltaS_t = deltaS_t+deltaS_p;
                        E_1 = E_1+E_2+E_3+E_4;
                        

                        pt_here->deltaS_new[j] = deltaS_t;
                        pt_here->Err_new[j] = E_1;

                        #if DetailedRecord
                        // {
                            pt_here->deltaS_p[j]= deltaS_p;
                        // }
                        #endif

                        if(E_1 < 0)
                        {
                            if(!muted)
                            {
                                printf("(%d,%d,%d),ERROR3: Area error less than 0: Err = %.16f, E1: %.16f, E2: %.16f, E3: %.16f\n",srcidx,idx,j,E_1,E_1-E_2-E_3,E_2,E_3);
                            }
                            srcs[srcidx].Break=1;
                            srcs[srcidx].SolveSucceed = 1;          // 用于跳过后续计算，会在最后一个SumArea设置为fail
                            #if DetailedRecord
                                srcs[srcidx].points_used = (batchidx+1) * BATCH_SIZE;
                            #endif  
                            continue;
                        }

                        #if DetailedRecord
                        // {
                            pt_here->E1[j]= E_1;
                            pt_here->E2[j]= E_2;
                            pt_here->E3[j]= E_3;                            
                        // }
                        #endif

                    }
                    else
                    {
                        pt_here->deltaS_new[j] = 0;
                        pt_here->Err_new[j] = 0;
                        #if DetailedRecord
                        // {
                            pt_here->E1[j]= 0;
                            pt_here->E2[j]= 0;
                            pt_here->E3[j]= 0;
                            pt_here->deltaS_p[j]=0; 
                        // }
                        #endif

                    }
                }
            }       
        }
    }
    // else{} 自己skip，不计算
}
template _Global void AreaErrBatchSrc4(int* srclist, src_ext_t<float>* srcs, int batchidx, bool muted);
template _Global void AreaErrBatchSrc4(int* srclist, src_ext_t<double>* srcs, int batchidx, bool muted);


// size <<<[(NSRCS-1) // 64 +1,1,1],64>>>
template<isFloating f_T>
_Global void AreaErrBatchSrcCross(int* srclist, src_ext_t<f_T>* srcs,const src_params_t<f_T>* src_params, int batchidx, int prev_Ncal, bool muted)
{
    int srcidx = threadIdx.x + blockIdx.x*blockDim.x;
    if(srcidx < prev_Ncal)
    {
        srcidx = srclist[srcidx];
        // printf("(srcidx:%d)(batchidx:%d)",srcidx,batchidx);

        int Ncross = srcs[srcidx].Ncross;
        if(Ncross>0)
        {
            int idx;        // physical
            bool ghost_direction;     // 1 means next, 0 means previous
            int j_test_p[2];
            // int j_test_g[2];

            // int paramsidx = (blockIdx.x)/NSRCS;
            complex_t<f_T> zs[2];
            complex_t<f_T> zeta;
            complex_t<f_T> dzs[2];
            f_T delta_theta;

            f_T params_lens[2];
            f_T rho;
            complex_t<f_T> loc_centre;

            f_T* wedge = &params_lens[0];       // x' ^ x''
            f_T E_1;
            f_T E_2;
            f_T E_3;            // 与VBB一致
            f_T deltaS_t;       // Trapezium approximation
            f_T deltaS_p;       // Parabolic correction
            complex_t<f_T> temp0;
            f_T temp1;

            complex_t<f_T> params[2*NBODY];

            int plus,minus;

            params_lens[0] = src_params->s;
            params_lens[1] = src_params->q;
            rho = src_params->shape.rho;
            loc_centre = srcs[srcidx].shape.loc_centre;

            ParamsReal2Complex(params_lens, params);   



            for(int i=0;i<Ncross;i++)
            {

                idx = srcs[srcidx].idx_cross[i];            


                if(idx<0)
                {
                    idx = -(idx+1);
                    ghost_direction = 0;
                }
                else
                {
                    ghost_direction = 1;
                }
                j_test_p[0] = srcs[srcidx].j_cross[i];
                j_test_p[1] = srcs[srcidx].margin_pts[idx].next_j[j_test_p[0]];


                if(srcs[srcidx].additional_cross[i])
                {
                    if(srcs[srcidx].margin_pts[idx].Nphys==3)
                    {
                        continue;                           
                    }
                    // else{}           do nothing, 这部分会在 AreaSrc里面计算掉
                    else{continue;}
                                     
                }

                if((j_test_p[0]<0)||(j_test_p[0]>4)||(j_test_p[1]<0)||(j_test_p[1]>4))
                {
                    if(!muted)
                    {
                        printf("invalid cross j: i: %d, batchidx: %d, srcidx: %d, idx: %d, j0: %d,j1: %d\n",i,batchidx,srcidx,idx,j_test_p[0],j_test_p[1]);
                    }
                    srcs[srcidx].Break=1;
                    srcs[srcidx].SolveSucceed = 1;          // 用于跳过后续计算，会在最后一个SumArea设置为fail
                    #if DetailedRecord
                        srcs[srcidx].points_used = (batchidx+1) * BATCH_SIZE;  
                    #endif
                    continue;
                }

                // 连接的 i,j拿来了

                if(srcs[srcidx].margin_pts[idx].images[j_test_p[0]].physical)
                {
                    // // gd==1, next is ghost, + connect to -, + is j0, - is j1
                    // // gd==0, prev is ghost, - connect to +, - is j0, + is j1

                    // zs[0] connect to zs[1]
                    zs[0] = srcs[srcidx].margin_pts[idx].images[j_test_p[0]].position;
                    zs[1] = srcs[srcidx].margin_pts[idx].images[j_test_p[1]].position;
                    zeta = srcs[srcidx].margin_pts[idx].position;                

                    wedge[0] = srcs[srcidx].margin_pts[idx].wedge[j_test_p[0]];
                    wedge[1] = srcs[srcidx].margin_pts[idx].wedge[j_test_p[1]];
                    dzs[0] = srcs[srcidx].margin_pts[idx].dz[j_test_p[0]];
                    dzs[1] = srcs[srcidx].margin_pts[idx].dz[j_test_p[1]];                  

                    deltaS_t = DeltaS1(zs[0],zs[1]);

                    delta_theta = abs_c(zs[1]-zs[0]) / sqrt_f_T(abs_c(dzs[0]*dzs[1]));

                    plus = int(!ghost_direction);
                    minus = int(ghost_direction);

                    deltaS_p = (wedge[plus] - wedge[minus]) * delta_theta*delta_theta*delta_theta /24.;
                    E_1 = fabs(wedge[plus] + wedge[minus]) * delta_theta*delta_theta*delta_theta / 48.;

                    temp0 = (zs[plus]-zs[minus])*(dzs[plus]-dzs[minus]);
                    temp1 = 2.*abs_c(zs[plus]-zs[minus])*sqrt_f_T(abs_c(dzs[plus]*dzs[minus]));
                    if(ghost_direction==1)  // next is ghost, destruction
                    {
                        temp0 = temp0 + temp1;
                    }
                    else    // prev is ghost, creation
                    {
                        temp0 = temp0 - temp1;
                    }
                    E_2 = 1.5*abs_c(temp0) * delta_theta;
                    E_3 = 0.1 * fabs(deltaS_p) * delta_theta*delta_theta;

                    deltaS_t = deltaS_t+deltaS_p;
                    E_1 = E_1+E_2+E_3;



                    srcs[srcidx].margin_pts[idx].deltaS_new[j_test_p[0]]= deltaS_t;
                    srcs[srcidx].margin_pts[idx].Err_new[j_test_p[0]]= E_1;                   
                    // srcs[srcidx].margin_pts[idx].deltaS[j_test_p[0]]= deltaS_t+deltaS_p;
                    // srcs[srcidx].margin_pts[idx].deltaS_Err[j_test_p[0]]= E_1+E_2+E_3;
                    #if DetailedRecord
                    // {
                        srcs[srcidx].margin_pts[idx].deltaS_p[j_test_p[0]]= deltaS_p;
                    // }
                    #endif
                    

                    

                    if(E_1 < 0)
                    {
                        if(!muted)
                        {
                            printf("(%d,%d,%d),ERROR: Area error less than 0: Err = %.16f, E1: %.16f, E2: %.16f, E3: %.16f\n",srcidx,idx,j_test_p[0],E_1,E_1-E_2-E_3,E_2,E_3);
                        }
                        srcs[srcidx].Break=1;
                        srcs[srcidx].SolveSucceed = 1;          // 用于跳过后续计算，会在最后一个SumArea设置为fail
                        #if DetailedRecord
                            srcs[srcidx].points_used = (batchidx+1) * BATCH_SIZE;
                        #endif       
                        continue;               
                    }

                    #if DetailedRecord
                    // {
                        srcs[srcidx].margin_pts[idx].E1[j_test_p[0]]= E_1;
                        srcs[srcidx].margin_pts[idx].E2[j_test_p[0]]= E_2;
                        srcs[srcidx].margin_pts[idx].E3[j_test_p[0]]= E_3;                          
                    // }
                    #endif
                  
                }


            }
        }
    }
}
template _Global void AreaErrBatchSrcCross(int* srclist, src_ext_t<float>* srcs,const src_params_t<float>* src_params, int batchidx, int prev_Ncal, bool muted);
template _Global void AreaErrBatchSrcCross(int* srclist, src_ext_t<double>* srcs,const src_params_t<double>* src_params, int batchidx, int prev_Ncal, bool muted);




template <isFloating f_T>
_Global void SumArea0(int* srclist, src_ext_t<f_T>* srcs, int batchidx, bool muted)
{
    int idx = threadIdx.x;    
    const int srcidx = srclist[int(blockIdx.x)];
    int next_src_idx;
    int temp_idx;
    int temp_Nre;
    int sum_length = 2;
    const int sum_size = (batchidx+1)*BATCH_SIZE;

    __dyn_shared__( int, data );
    f_T* deltaS = (f_T*) data;
    f_T* Err = (f_T*) &data[blockDim.x/sizeof(int)*sizeof(f_T)];    // size: blockDim.x = min(1024,sum_size)

    f_T temp_deltaS;
    f_T temp_Err;
    int pointidx;


    if(srcs[srcidx].SolveSucceed == false){

        deltaS[idx] = 0;
        Err[idx] = 0;

        for(int i=0;i<((sum_size-1)/blockDim.x +1);i++)        // ceil(NPOINTS / blockDim.x)
        {
            pointidx = i*blockDim.x + idx;
            if(pointidx < sum_size)
            {
                for(int j=0;j<5;j++)
                {
                    srcs[srcidx].margin_pts[pointidx].deltaS_Err[j] = srcs[srcidx].margin_pts[pointidx].Err_new[j];
                    srcs[srcidx].margin_pts[pointidx].deltaS[j] = srcs[srcidx].margin_pts[pointidx].deltaS_new[j];
                }
            }
        }

        __syncthreads();


        for(int i=0;i<((sum_size-1)/blockDim.x +1);i++)        // ceil(NPOINTS / blockDim.x)
        {
            pointidx = i*blockDim.x + idx;
            if(pointidx < sum_size)
            {
                if(!srcs[srcidx].margin_pts[pointidx].skip)
                {
                    srcs[srcidx].margin_pts[pointidx].changed = 0;

                    temp_deltaS = 0;
                    temp_Err = 0;
                    next_src_idx = NextSrcIdx(srcs[srcidx].margin_pts,pointidx);
                    for(int j=0;j<LENGTH-1;j++)
                    {
                        temp_idx = srcs[srcidx].margin_pts[pointidx].next_idx[j];
                        if(temp_idx == next_src_idx)    // 自己连下一个
                        {
                            temp_deltaS += srcs[srcidx].margin_pts[pointidx].deltaS[j];
                            temp_Err += srcs[srcidx].margin_pts[pointidx].deltaS_Err[j];
                        }
                        if(temp_idx == pointidx)         // 自己连自己，那自己肯定有5个实像
                        {
                            temp_Nre = 0;
                            for(int jj=0;jj<LENGTH-1;jj++)
                            {
                                if(srcs[srcidx].margin_pts[next_src_idx].images[jj].physical){temp_Nre++;}
                            }
                            if(temp_Nre==3)         // 下一个源点只有三个解，说明自己和下一个源点之间有散焦线，误差权重归于自己（即自己和下一个源点之间）
                            {
                                temp_deltaS += srcs[srcidx].margin_pts[pointidx].deltaS[j];
                                temp_Err += srcs[srcidx].margin_pts[pointidx].deltaS_Err[j];
                            }
                            // 如果下个源点是五个解，说明焦散线在自己和上一个源点之间，由上一个源点负责
                        }

                        temp_idx = srcs[srcidx].margin_pts[next_src_idx].next_idx[j];
                        if(temp_idx == pointidx)         // 下一个连自己
                        {
                            temp_deltaS += srcs[srcidx].margin_pts[next_src_idx].deltaS[j];
                            temp_Err += srcs[srcidx].margin_pts[next_src_idx].deltaS_Err[j];
                        }
                        if(temp_idx == next_src_idx)     // 下一个连下一个
                        {
                            temp_Nre = 0;
                            for(int jj=0;jj<LENGTH-1;jj++)
                            {
                                if(srcs[srcidx].margin_pts[pointidx].images[jj].physical){temp_Nre++;}
                            }
                            if(temp_Nre==3)         // 自己源点只有三个解，说明下一个源点和自己之间有散焦线，误差权重归于自己（即下一个源点和自己之间）
                            {
                                temp_deltaS += srcs[srcidx].margin_pts[next_src_idx].deltaS[j];
                                temp_Err += srcs[srcidx].margin_pts[next_src_idx].deltaS_Err[j];
                            }
                            // 如果自己源点是五个解，说明焦散线在下一个和下下一个源点之间，由下一个源点负责                 
                        }                     
                        // }


                        if(srcs[srcidx].margin_pts[pointidx].deltaS_Err[j]<0)
                        {
                            if(!muted)
                            {
                                printf("Negative piece error: srcidx: %d, idx: %d, j: %d, Err= %.16f\n",srcidx,idx,j,srcs[srcidx].margin_pts[idx].deltaS_Err[j]);
                            }
                        }
                        
                    }
                    srcs[srcidx].margin_pts[pointidx].error_interval = temp_Err;
                    srcs[srcidx].margin_pts[pointidx].area_interval = temp_deltaS;
                    deltaS[idx] += temp_deltaS;
                    Err[idx] += temp_Err;
                }
                else
                {
                    srcs[srcidx].margin_pts[pointidx].error_interval = 0;
                    srcs[srcidx].margin_pts[pointidx].area_interval = 0;
                }
            }

        }


        __syncthreads();



        while(sum_length<(blockDim.x*2))
        {
            if(((idx & (sum_length-1))==0) && ((idx + sum_length/2) < blockDim.x))    // sum_length = 2^n, idx%sum_length = idx & (sum_length-1)
            {
                deltaS[idx] += deltaS[idx + sum_length/2];
                // deltaS_p[idx] += deltaS_p[idx + sum_length/2];
                Err[idx] += Err[idx + sum_length/2];

                if(Err[idx]<0)
                {
                    if(!muted)
                    {
                        printf("Negative total error: srcidx: %d, idx: %d, j: %d, Err= %.16f\n",srcidx,idx,sum_length,Err[idx]);
                    }
                }
            }

            sum_length *= 2;
            __syncthreads();
        }
        if(idx == 0)
        {
            srcs[srcidx].Mag = fabs(deltaS[0]) / Area_src(srcs[srcidx].shape);
            srcs[srcidx].Err = Err[0] / Area_src(srcs[srcidx].shape);
            srcs[srcidx].Area = deltaS[0];
            srcs[srcidx].Err_A = Err[0];


            if(Err[0]/fabs(deltaS[0])<=RELTOL+RelErrAllowed)
            {
                srcs[srcidx].SolveSucceed = true;
                // printf("%d: pass, rel err = %.12f\n",srcidx,Err[0]/abs(deltaS[0]));
                #if DetailedRecord
                    srcs[srcidx].points_used = sum_size;
                #endif
            }
            if((sum_size == NPOINTS)&&(Err[0]/fabs(deltaS[0])>=RELTOL+RelErrAllowed))
            {
                if(!muted)
                {
                    printf("srcidx: %d: fail, rel err = %.12f\n",srcidx,Err[0]/fabs(deltaS[0]));
                }
                #if DetailedRecord
                    srcs[srcidx].points_used = sum_size;
                #endif
            }

        }
        
    }      
    if((srcs[srcidx].Break==1))
    {
        if(idx==0)
        {
            srcs[srcidx].SolveSucceed = false;
            if(!muted)
            {
                printf("srcidx: %d: fail, break\n",srcidx);
            }

            srcs[srcidx].Mag = -1;
            srcs[srcidx].Err = -1;
            srcs[srcidx].Area = -0;
            srcs[srcidx].Err_A = -0;
          
        }

    }        

}
template _Global void SumArea0(int* srclist, src_ext_t<float>* srcs, int batchidx, bool muted);
template _Global void SumArea0(int* srclist, src_ext_t<double>* srcs, int batchidx, bool muted);


// 只算 3*BATCHSIZE = 96 个
// 甚至更少一点儿，但位置是96个
// size: <<<NSRCS, BATCH_SIZE, 4*BATCH_SIZE>>>
// 先都放进shared memory，计算好每个点的误差存下来
// 至多2N个区间
template <isFloating f_T>
_Global void SumArea3(int* srclist, src_ext_t<f_T>* srcs, int batchidx, bool muted)
{
    int idx = threadIdx.x;    
    const int srcidx = srclist[int(blockIdx.x)];


    if(srcs[srcidx].SolveSucceed == false){
        int sum_length = 2;
        // const int sum_size = 2*blockDim.x;
        int pointidx;
        int idx_list[2];    
        int temp0;
        int temp_j;
        src_pt_t<f_T>* pt_here;    
        f_T Area,Error;
        f_T tempS;
        int next_src_idx;
        int temp_idx;
        int temp_Nre;
        f_T old_err;
        f_T old_area;

        __dyn_shared__( int, data );
        f_T* deltaS = (f_T*) data;
        f_T* Err = (f_T*) &data[2*blockDim.x/sizeof(int)*sizeof(f_T)];    // size: = 2*blockDim.x;

        for(int loop2=0;loop2<2;loop2++)
        {
            deltaS[idx + loop2*blockDim.x] = 0;
            Err[idx + loop2*blockDim.x] = 0;            
        }

        // 按区间来看，只有 loop0 和 loop1 对应的区间权重被改了
        idx_list[1] = batchidx*BATCH_SIZE + idx;
        idx_list[0] = PrevSrcIdx(srcs[srcidx].margin_pts, idx_list[1]);
        if(idx_list[0]>= batchidx*BATCH_SIZE){idx_list[0]=-1;}

        if(!srcs[srcidx].margin_pts[idx_list[1]].skip)
        {
            for(int loop2=0;loop2<2;loop2++)
            {
                old_err = 0;
                old_area = 0;
                pointidx = idx_list[loop2];
                pt_here = &srcs[srcidx].margin_pts[pointidx];

                // if(pt_here->skip)
                // {
                //     next_src_idx = -1;
                //     continue;
                // }
                if(loop2!=1)
                {
                    old_err = srcs[srcidx].margin_pts[pointidx].error_interval;
                    old_area = srcs[srcidx].margin_pts[pointidx].area_interval;
                }
                if(pointidx<0){temp0 = 0;}
                else
                {
                    temp0 = atomicExch(&(pt_here->changed),0);
                }

                if(temp0==3){printf("strange\n");}
                if(temp0==4)
                {

                    // 如果loop0==-1，不进入这个分支
                    next_src_idx = NextSrcIdx(srcs[srcidx].margin_pts,pointidx);
                    
                    for(int j=0;j<5;j++)
                    {
                        temp_idx = srcs[srcidx].margin_pts[pointidx].next_idx[j];
                        if(temp_idx == next_src_idx)
                        {
                            tempS = srcs[srcidx].margin_pts[pointidx].deltaS_new[j];
                            deltaS[idx + loop2*BATCH_SIZE] += tempS;
                            srcs[srcidx].margin_pts[pointidx].deltaS[j] = tempS;

                            tempS = srcs[srcidx].margin_pts[pointidx].Err_new[j];
                            Err[idx + loop2*BATCH_SIZE] += tempS;
                            srcs[srcidx].margin_pts[pointidx].deltaS_Err[j] = tempS;
                        }
                        if(temp_idx == pointidx)         // 自己连自己，那自己肯定有5个实像
                        {
                            temp_Nre = 0;
                            for(int jj=0;jj<LENGTH-1;jj++)
                            {
                                if(srcs[srcidx].margin_pts[next_src_idx].images[jj].physical){temp_Nre++;}
                            }
                            if(temp_Nre==3)         // 下一个源点只有三个解，说明自己和下一个源点之间有散焦线，误差权重归于自己（即自己和下一个源点之间）
                            {
                                tempS = srcs[srcidx].margin_pts[pointidx].deltaS_new[j];
                                deltaS[idx + loop2*BATCH_SIZE] += tempS;
                                srcs[srcidx].margin_pts[pointidx].deltaS[j] = tempS;

                                tempS = srcs[srcidx].margin_pts[pointidx].Err_new[j];
                                Err[idx + loop2*BATCH_SIZE] += tempS;
                                srcs[srcidx].margin_pts[pointidx].deltaS_Err[j] = tempS;                           
                            }
                            // 如果下个源点是五个解，说明焦散线在自己和上一个源点之间，由上一个源点负责
                        }

                        temp_idx = srcs[srcidx].margin_pts[next_src_idx].next_idx[j];
                        if(temp_idx == pointidx)         // 下一个连自己
                        {
                            tempS = srcs[srcidx].margin_pts[next_src_idx].deltaS_new[j];
                            deltaS[idx + loop2*BATCH_SIZE] += tempS;
                            srcs[srcidx].margin_pts[next_src_idx].deltaS[j] = tempS;

                            tempS = srcs[srcidx].margin_pts[next_src_idx].Err_new[j];
                            Err[idx + loop2*BATCH_SIZE] += tempS;
                            srcs[srcidx].margin_pts[next_src_idx].deltaS_Err[j] = tempS;                    
                        }
                        if(temp_idx == next_src_idx)     // 下一个连下一个
                        {
                            temp_Nre = 0;
                            for(int jj=0;jj<LENGTH-1;jj++)
                            {
                                if(srcs[srcidx].margin_pts[pointidx].images[jj].physical){temp_Nre++;}
                            }
                            if(temp_Nre==3)         // 自己源点只有三个解，说明下一个源点和自己之间有散焦线，误差权重归于自己（即下一个源点和自己之间）
                            {
                                tempS = srcs[srcidx].margin_pts[next_src_idx].deltaS_new[j];
                                deltaS[idx + loop2*BATCH_SIZE] += tempS;
                                srcs[srcidx].margin_pts[next_src_idx].deltaS[j] = tempS;

                                tempS = srcs[srcidx].margin_pts[next_src_idx].Err_new[j];
                                Err[idx + loop2*BATCH_SIZE] += tempS;
                                srcs[srcidx].margin_pts[next_src_idx].deltaS_Err[j] = tempS; 

                            }
                            // 如果自己源点是五个解，说明焦散线在下一个和下下一个源点之间，由下一个源点负责                 
                        }  
                    }
                    srcs[srcidx].margin_pts[pointidx].error_interval = Err[idx + loop2*BATCH_SIZE];
                    srcs[srcidx].margin_pts[pointidx].area_interval = deltaS[idx + loop2*BATCH_SIZE];
                    Err[idx + loop2*BATCH_SIZE] -= old_err;
                    deltaS[idx + loop2*BATCH_SIZE] -= old_area;
                }
                // else{}  // do nothing, error_interval needn't be changed
            }
                
            deltaS[idx] += deltaS[idx + blockDim.x];
            Err[idx] += Err[idx + blockDim.x];
        }
        else
        {
            next_src_idx = -1;
        }
        __syncthreads();
        if(next_src_idx>=0)
        {
            srcs[srcidx].margin_pts[next_src_idx].changed = 0;      // loop3==2, from 4 to 0;
        }
          

        while(sum_length<(2*blockDim.x))        // 只对前 BATCH_SIZE求和，因为后面的已经加过来了
        {
            if(((idx & (sum_length-1))==0) && ((idx + sum_length/2) < blockDim.x))    // sum_length = 2^n, idx%sum_length = idx & (sum_length-1)
            {
                deltaS[idx] += deltaS[idx + sum_length/2];
                // deltaS_p[idx] += deltaS_p[idx + sum_length/2];
                Err[idx] += Err[idx + sum_length/2];

            }
            sum_length *= 2;
            __syncthreads();
        }
        if(idx == 0)
        {
            for(int i=0;i<srcs[srcidx].Ncross;i++)
            {
                temp0 = srcs[srcidx].idx_cross[i];
                if(temp0<0)
                {
                    temp0 = -(temp0+1);
                }
                temp_j = srcs[srcidx].j_cross[i];
                tempS = srcs[srcidx].margin_pts[temp0].deltaS_new[temp_j];
                Area = srcs[srcidx].margin_pts[temp0].deltaS[temp_j];           // 借用一下空间
                if(Area!=tempS)
                {
                    deltaS[0] += tempS - Area;
                    srcs[srcidx].margin_pts[temp0].deltaS[temp_j] = tempS;
                }
                tempS = srcs[srcidx].margin_pts[temp0].Err_new[temp_j];
                Error = srcs[srcidx].margin_pts[temp0].deltaS_Err[temp_j];
                if(Error!=tempS)
                {
                    Err[0] += tempS - Error;
                    srcs[srcidx].margin_pts[temp0].deltaS_Err[temp_j] = tempS;
                }
            }

            Area = deltaS[0] + srcs[srcidx].Area;
            Error = Err[0] + srcs[srcidx].Err_A;
            tempS = srcs[srcidx].src_area;
        

            if(Error<0)
            {
                if(!muted)
                {
                    printf("Negative total error: batchidx: %d, srcidx: %d, idx: %d, j: %d, Err= %.16f, Err0: %.16f\n",batchidx,srcidx,idx,sum_length,Error,Err[0]);
                }
            }            
            srcs[srcidx].Mag = fabs(Area) / tempS;
            srcs[srcidx].Err = Error / tempS;
            srcs[srcidx].Area = Area;
            srcs[srcidx].Err_A = Error;

            if(Error/fabs(Area)<=RELTOL+RelErrAllowed)
            {
                srcs[srcidx].SolveSucceed = true;
                // printf("%d: pass, rel err = %.12f\n",srcidx,Err[0]/abs(deltaS[0]));
                #if DetailedRecord
                    srcs[srcidx].points_used = sum_size;
                #endif
            }
            if(((batchidx+1)*BATCH_SIZE == NPOINTS)&&(Error/fabs(Area)>=RELTOL+RelErrAllowed))          // 终止条件
            {
                if(!muted)
                {
                    printf("srcidx: %d: fail, rel err = %.12f\n",srcidx,Error/fabs(Area));
                }
                #if DetailedRecord
                    srcs[srcidx].points_used = sum_size;
                #endif
            }

        }
        
    }      
    if((srcs[srcidx].Break==1))
    {
        if(idx==0)
        {
            srcs[srcidx].SolveSucceed = false;
            // printf("%d: fail, rel err = %.12f\n",srcidx,Err[0]/abs(deltaS[0]));
            if(!muted)
            {
                printf("srcidx: %d: fail, break\n",srcidx);
            }

            srcs[srcidx].Mag = -1;
            srcs[srcidx].Err = -1;
            srcs[srcidx].Area = -0;
            srcs[srcidx].Err_A = -0;
           
        }

    }        


}
template _Global void SumArea3(int* srclist, src_ext_t<float>* srcs, int batchidx, bool muted);
template _Global void SumArea3(int* srclist, src_ext_t<double>* srcs, int batchidx, bool muted);


// size:<<<(NWALKERS,1,ceil(NSRCS/64)),min(NSRCS,64)>>>
template<isFloating f_T>
_Global void QuadruTst(src_ext_t<f_T>* srcs, const src_params_t<f_T>* src_params, const int Nsrcs)
{
    
    // int srcidx = srclist[int(threadIdx.x + blockDim.x*blockIdx.z)];  // srcs连续排布
    int srcidx = threadIdx.x + blockDim.x*blockIdx.x;

    if(srcidx<NSRCS){
        
        if(srcs[srcidx].SolveSucceed == false){


            complex_t<f_T> loc_src = srcs[srcidx].shape.loc_centre;
            complex_t<f_T> params[ 2*NBODY ];
        
            complex_t<f_T> coeffs[LENGTH];

            f_T temp_jacobi[LENGTH-1];
            complex_t<f_T> lambdas[LENGTH-1];
            f_T params_lens[2];
            bool physical[LENGTH-1];
            bool fail=0;
            int phys_tst = 0;
            bool test=0;
            f_T& Rho = coeffs[0].im;
            
            
            f_T& muQC = coeffs[0].re;
            f_T& muPS = coeffs[1].re;

            complex_t<f_T> z_hat;
            f_T& GT = temp_jacobi[0];   // Ghost Test
            f_T& safedist = temp_jacobi[1]; // Planetary Test
            complex_t<f_T>& zeta_pc = coeffs[2];
            f_T& delta_pc2 = temp_jacobi[2];
            f_T moving;

            f_T physical_err[5];        // 临时输出查看用，本质上不需要
            // 对于每个线程，上述需要
            // float: 4+4+4 + 8*(1+4+6+5) + 4*(2+5+5) + 12*5+1 = 249 Byte, Ampere架构，单个SM里有4个warp，每个warp 65536 Byte 的寄存器。
            // double: 4+4+4 + 16*(16) + 8*(12) + 24*5 +1 = 485 Byte  BATCH_SIZE最多64，建议32

            // 二体特例
            params_lens[0] = src_params->s;
            params_lens[1] = src_params->q;

            ParamsReal2Complex(params_lens, params);
            moving = PolyCoeffTwinkle(params_lens[0],1.f/(1.f+params_lens[1]), params_lens[1]/(1.f+params_lens[1]), loc_src, coeffs);

            fail = cmplx_roots_gen(lambdas,coeffs,LENGTH-1,true,false);

            for(int j=0;j<5;j++)
            {
                lambdas[j].re += moving;
            }

            // 至此，解出中心位置的5个解，存在lambdas里
            if(! fail)
            {
                for(int j=0;j<LENGTH-1;j++)
                {
                    physical[j] = isPhysical(lambdas[j],params,loc_src,physical_err[j]);
                    phys_tst += physical[j];
                    temp_jacobi[j] = Jacobi(lambdas[j],params);
                }

                // 对不同的实像数量分别考虑，由于数值误差，可能有一些点不小心跨过了临界线，导致实像数量不对。但这种情况都是过于接近临界线的极端情况，本来就通过不了Qtst
                if(phys_tst==3 || phys_tst==5)
                {
                    muQC = 0;
                    muPS = 0;
                    Rho = srcs[srcidx].shape.rho;
                    // 先进行四极检验 (VBB2,chap 5.1)
                    for(int j=0;j<LENGTH-1;j++)
                    {
                        if(physical[j]==true)
                        {
                            muQC += fabs(Mu_QC(params,lambdas[j]) * (Rho+1e-3) * (Rho+1e-3));
                            muPS += 1/fabs(temp_jacobi[j]);
                        }
                    }
                    srcs[srcidx].Mag_P = muPS;
                    srcs[srcidx].Err_Mag_P = muQC;
                    if(muQC * C_Q < RELTOL * muPS)
                    {
                        test = true;
                    }

                    // 如果有两个虚像，还要进行虚像检验(VBB2,chap 5.2)
                    if(phys_tst==3 && test==true)
                    {
                        GT = 0;
                        for(int j=0;j<LENGTH-1;j++)
                        {
                            if(physical[j]==false)
                            {
                                z_hat = conj(loc_src) - f_zbar(lambdas[j],params);     // f_zbar: input z, return f(z_conj)
                                GT = max(GT, GhostTst(params,lambdas[j],z_hat));
                            }
                        }
                        test = (GT* (Rho+1e-3)*2*C_G < 1);   // 能到这里说明test已经是true了，接下来test判断只看GT

                    }

                    // 安全距离检测（VBB2, chap 5.3）
                    safedist = 10;
                    if(params_lens[1]<0.01)     // q<0.01
                    {
                        zeta_pc = complex_t<f_T>(-1./params_lens[0],0) + params[1];     // slightly different from VBB2 Eq49, because the frame centre is different. params[1] is z_1
                        delta_pc2 = 9*params_lens[1] / params_lens[0] / params_lens[0];
                        safedist = norm(loc_src-zeta_pc) - C_P*delta_pc2;
                    }
                    if(safedist <= C_P * Rho*Rho)
                    {
                        test = false;
                    }              



                    if(test == true)
                    {
                        srcs[srcidx].SolveSucceed = true;
                        srcs[srcidx].Mag = muPS;
                        srcs[srcidx].Err = muQC;
                        srcs[srcidx].Break = false;
                    }
                }
                else
                {
                    // 点源求解实像数量不对
                    srcs[srcidx].Mag_P = -2;
                }
                for(int iii=0;iii<5;iii++)
                {
                    srcs[srcidx].roots[iii] = lambdas[iii];
                    #if DetailedRecord
                    // {
                        srcs[srcidx].phys_tst[iii] = physical[iii];
                        srcs[srcidx].phys_err[iii] = physical_err[iii];
                    // }
                    #endif
                }                
            }
            else
            {
                // 点源求解失败
                srcs[srcidx].Mag_P = -3;
            }

            if(test==false)
            {
                srcs[srcidx].SolveSucceed = false;
                srcs[srcidx].Break = false;                
            }

        }

    }

    
   return;
}
template _Global void QuadruTst(src_ext_t<float>* srcs, const src_params_t<float>* src_params, const int Nsrcs);
template _Global void QuadruTst(src_ext_t<double>* srcs, const src_params_t<double>* src_params, const int Nsrcs);

// size:<<<(NWALKERS,1,ceil(NSRCS/1024)),min(NSRCS,1024)>>>; 光变曲线点数可能比1024更多
template <isFloating f_T>
_Global void SetPath(src_ext_t<f_T>*srcs, const f_T* time, const src_params_t<f_T>* move_params, int Nsrcs)
{
    int srcidx = threadIdx.x + blockDim.x*blockIdx.z;  // srcs连续排布，一个walker的srcs连续分布，不同walkers的srcs作为整体连续排布
    complex_t<f_T> shape_centre;
    
    if(srcidx < Nsrcs)
    {
        PathLoc(move_params[0],time[srcidx], shape_centre);

        srcs[srcidx].shape.rho = move_params[0].shape.rho;
        srcs[srcidx].shape.loc_centre = shape_centre;
        srcs[srcidx].src_area = Area_src(srcs[srcidx].shape);
    }
}
template _Global void SetPath(src_ext_t<float>*srcs, const float* time, const src_params_t<float>* move_params, int Nsrcs);
template _Global void SetPath(src_ext_t<double>*srcs, const double* time, const src_params_t<double>* move_params, int Nsrcs);

// size:<<<(NWALKERS,1,ceil(NSRCS/1024)),min(NSRCS,1024)>>>; 光变曲线点数可能比1024更多
template <isFloating f_T>
_Global void SetRhoLoc(src_ext_t<f_T>*srcs, const f_T* rhos, const complex_t<f_T>* zetas, int Nsrcs)
{
    int srcidx = threadIdx.x + blockDim.x*blockIdx.z;  // srcs连续排布，一个walker的srcs连续分布，不同walkers的srcs作为整体连续排布
    // complex_t<f_T> shape_centre;
    
    if(srcidx < Nsrcs)
    {

        srcs[srcidx].shape.rho = rhos[srcidx];
        srcs[srcidx].shape.loc_centre = zetas[srcidx];
        srcs[srcidx].src_area = Area_src(srcs[srcidx].shape);
    }
}
template _Global void SetRhoLoc(src_ext_t<float>*srcs, const float* rhos, const complex_t<float>* zetas, int Nsrcs);
template _Global void SetRhoLoc(src_ext_t<double>*srcs, const double* rhos, const complex_t<double>* zetas, int Nsrcs);


// size:<<<NaxisY,NaxisX>>>; 
template <isFloating f_T>
_Global void SetArray(src_ext_t<f_T>*srcs, const src_params_t<f_T>* move_params, const int NaxisX, const int NaxisY, const f_T xmin,const f_T xmax,const f_T ymin,const f_T ymax, int Nsrcs)
{
    int x_idx = threadIdx.x;
    int y_idx = blockIdx.x;
    int srcidx = x_idx + y_idx * NaxisX;    // x连续分布
    complex_t<f_T> shape_centre;

    shape_centre.re = ((f_T(x_idx)+0.5) / NaxisX) * (xmax - xmin) + xmin;
    shape_centre.im = ((f_T(y_idx)+0.5) / NaxisY) * (ymax - ymin) + ymin;
    
    if(srcidx < Nsrcs)
    {
        srcs[srcidx].shape.rho = move_params[0].shape.rho;
        srcs[srcidx].shape.loc_centre = shape_centre;
        srcs[srcidx].src_area = Area_src(srcs[srcidx].shape);
    }


}
template _Global void SetArray(src_ext_t<float>*srcs, const src_params_t<float>* move_params, const int NaxisX, const int NaxisY, const float xmin,const float xmax,const float ymin,const float ymax, int Nsrcs);
template _Global void SetArray(src_ext_t<double>*srcs, const src_params_t<double>* move_params, const int NaxisX, const int NaxisY, const double xmin,const double xmax,const double ymin,const double ymax, int Nsrcs);

// size:<<<1,1>>>; 
template <isFloating f_T>
_Global void SetPoint(src_ext_t<f_T>*srcs, const src_params_t<f_T>* move_params, const f_T x, const f_T y)
{
    int srcidx = 0;    // x连续分布
    complex_t<f_T> shape_centre;

    shape_centre.re = x;
    shape_centre.im = y;
    

    if(srcidx < NSRCS)
    {
        srcs[srcidx].shape.rho = move_params[0].shape.rho;
        srcs[srcidx].shape.loc_centre = shape_centre;
    }


}
template _Global void SetPoint(src_ext_t<float>*srcs, const src_params_t<float>* move_params, const float x, const float y);
template _Global void SetPoint(src_ext_t<double>*srcs, const src_params_t<double>* move_params, const double x, const double y);



// size: <<<(NSRCS,1,1),BATCH_SIZE*(batchidx+1),(BATCH_SIZE)*(batchidx+1)*sizeof(f_T) + (BATCH_SIZE)*sizeof(int)>>>
// 利用前 BATCH_SIZE*(batchidx) 个点，得到 [BATCH_SIZE*(batchidx),BATCH_SIZE*(batchidx+1)) 个点的位置
// int: 7; f_T:5; bool:1
// float: 49 Byte
// double: 69 Byte
// Npts_now: BATCH_SIZE*(batchidx); Npts_next: BATCH_SIZE*(batchidx+1)
// size: <<<(NSRCS,1,1),Npts_next,(BATCH_SIZE)*Npts_now*sizeof(f_T) + (BATCH_SIZE)*sizeof(int)>>>
template <isFloating f_T>
_Global void AdaptiveLocLoop(int* srclist, src_ext_t<f_T>*srcs, int batchidx, bool muted)
{
    const int srcidx = srclist[int(blockIdx.x)];
    // extern __shared__ int shared_data[];                 // 共享内存大小不能用f_T! 可以先用int分配大小（反正传进来的是字节数），再转化指针类型
    __dyn_shared__(int, shared_data);
    f_T* ErrorWeight = (f_T*) shared_data;                  // size = BATCH_SIZE*(batchidx)
    int* ParentIdx = (int*) &shared_data[BATCH_SIZE*batchidx*sizeof(f_T)/sizeof(int)];


    int idx = threadIdx.x;
    int temp_idx;
    int temp_Nre;
    int piece_length;
    f_T normalize;
    f_T& quantile = normalize;
    int& left_idx = temp_idx;
    int& right_idx = temp_Nre;
    f_T left_Q;
    f_T right_Q;
    f_T quantile_left;
    f_T Q;

    f_T right_idx_Q;
    f_T left_idx_Q;

    for(int i=0;i<((batchidx*BATCH_SIZE-1)/blockDim.x +1);i++)        // ceil(NPOINTS / blockDim.x)
    {
        // printf("check %d\n",i);
        if(i*blockDim.x + idx < batchidx*BATCH_SIZE)
        {
            ErrorWeight[i*blockDim.x + idx] = 0;



            // 如果该点求解失败，那就没它的事，权重为0（相当于这个点不存在，将由附近的点负责查找）（而且它也压根没算面积误差）
            // 权重是“当前源点和下个源点之间的所有像面积之和”
            // if((srcs[srcidx].margin_pts[idx].skip == false) || (idx >= blockDim.x - BATCH_SIZE))               // 把最后的点也包进来，这些点还没有skip参数。 
            // {
            if((i*blockDim.x + idx)<(batchidx*BATCH_SIZE) && (srcs[srcidx].margin_pts[(i*blockDim.x + idx)].skip == false))
            {
                ErrorWeight[(i*blockDim.x + idx)] = srcs[srcidx].margin_pts[(i*blockDim.x + idx)].error_interval;
                ErrorWeight[(i*blockDim.x + idx)] += srcs[srcidx].margin_pts[(i*blockDim.x + idx)].special_Err;
                srcs[srcidx].margin_pts[(i*blockDim.x + idx)].special_Err = 0.;  

                ErrorWeight[(i*blockDim.x + idx)] /= fabs(srcs[srcidx].Area);          // 使用相对误差

                // ErrorWeight[idx] = (ErrorWeight[idx] * ErrorWeight[idx] );          // 使用平方根权重，因为误差大小的对delta_theta是超线性相关，误差更敏感，theta距离更迟钝

                ErrorWeight[(i*blockDim.x + idx)] = max(f_T(0.),ErrorWeight[(i*blockDim.x + idx)] - RELTOL/(batchidx*BATCH_SIZE));        // 如果自己的部分误差小于均值了，那算作过关

                // printf("(%d,%d),%.16f\n",srcidx,(i*blockDim.x + idx),ErrorWeight[(i*blockDim.x + idx)]);
                
                if(ErrorWeight[(i*blockDim.x + idx)]<0)
                {
                    if(!muted)
                    {
                        printf("ERROR: srcidx: %d, idx: %d, weight less than zero: %.16f\n",srcidx,(i*blockDim.x + idx),ErrorWeight[(i*blockDim.x + idx)]);
                    }
                }
                
            }
        }
    }

    // input to shared memory finished
    __syncthreads();


    // 多段求和，每段用二分，再累加
    // 单段长度：blockDim.x; 总长度：BATCH_SIZE*(batchidx)
    for(int i=0;i<((batchidx*BATCH_SIZE-1)/blockDim.x +1);i++)        // ceil(batchidx*BATCH_SIZE / blockDim.x)
    {
        piece_length = 2;
        while(piece_length < 2*blockDim.x)
        {
            if((((i*blockDim.x + idx) & (piece_length-1)) >= piece_length/2) && ((i*blockDim.x + idx) < BATCH_SIZE*batchidx))      // piece_length = 2^n, then idx % piece_length = idx & (piece_length-1)
            {
                ErrorWeight[(i*blockDim.x + idx)] += ErrorWeight[(i*blockDim.x + idx) - ((i*blockDim.x + idx) & (piece_length-1)) + piece_length/2 -1];
            }
            __syncthreads();
            piece_length *= 2;
        }
    }
    // 每一段加前一段的尾巴，跳过第0段，从第1段开始，
    for(int i=1;i<((batchidx*BATCH_SIZE-1)/blockDim.x +1);i++)        // ceil(batchidx*BATCH_SIZE / blockDim.x)
    {
        // printf("check %d\n",i);
        if((i*blockDim.x + idx)<(batchidx*BATCH_SIZE))
        {
            ErrorWeight[(i*blockDim.x + idx)] += ErrorWeight[(i*blockDim.x - 1)];
        }
        __syncthreads();
    }

    // 归一化
    normalize = ErrorWeight[BATCH_SIZE*batchidx - 1];
    __syncthreads();

    // 如果noramlize比较大，归一化
    if(normalize>=RelErrAllowed)        // 与idx无关的分支，可以
    {
        for(int i=0;i<((batchidx*BATCH_SIZE-1)/blockDim.x +1);i++)        // ceil(batchidx*BATCH_SIZE / blockDim.x)
        {
            if((i*blockDim.x + idx) < BATCH_SIZE*batchidx)
            {
                ErrorWeight[(i*blockDim.x + idx)] /= normalize;
            }
        }    
    }
    // 如果normalize非常小，即面积误差累加比阈值只超出非常小，为了规避数值误差（乃至NaN），将这种情况计作求解成功
    // 这种情况换成均匀分布！由successcheck判定是否要break
    else
    {
        for(int i=0;i<((batchidx*BATCH_SIZE-1)/blockDim.x +1);i++)        // ceil(batchidx*BATCH_SIZE / blockDim.x)
        {
            if((i*blockDim.x + idx) < BATCH_SIZE*batchidx)
            {
                ErrorWeight[(i*blockDim.x + idx)] = f_T(i*blockDim.x + idx + 1) / f_T(BATCH_SIZE*batchidx);
            }
        }              
    }     

    __syncthreads();


    // ErrorWeight[i]意味着内存地址 i 的点和内存地址 next_i（不是i+1）单点之间的权重
    // 我只要连接left
    // 前64个点插值Q

    if(((idx < (BATCH_SIZE))) && (srcs[srcidx].SolveSucceed==false))
    {

    // 反向插值回去，前64个线程负责，以上述权重找64个新点
    // +0.5 避免总是插值到第0个点上
        quantile = (f_T(idx)+f_T(0.5)) / f_T(BATCH_SIZE);
        left_idx = BisearchLeft(ErrorWeight,quantile,batchidx*BATCH_SIZE);      // 取值范围：[0,len-1]，即parent节点，相邻的left_idx在内存中连续，但是Q不连续
        right_idx = NextSrcIdx(srcs[srcidx].margin_pts, left_idx);

    //    printf("%.6f ",quantile);
        // 这里left_idx right_idx存的是左右points的指标
        // quantile_left的指标在points指标（即Q指标）-1，ErrW[0]对应pts[1]的分度值
        left_Q = srcs[srcidx].margin_pts[left_idx].Q;       
        if(left_idx == 0)
        {
            quantile_left = 0;// 左点需要考虑 第0点和第1点之间的情况
        }
        else
        {
            quantile_left = ErrorWeight[left_idx-1];
        }
        
        right_Q = srcs[srcidx].margin_pts[right_idx].Q;
        if(right_Q == 0)
        {
            right_Q = 1;             // 右点需要将01认同，以满足rightQ总是大于leftQ
        }
        // quantile_right = ErrorWeight[left_idx -1 +1];       // 注意quantile时right和left的差别
        // 不创建变量（塞不下了），直接调用共享内存的数据放在缓存里
        // 线性插值
        Q = left_Q + (right_Q - left_Q) * (quantile-quantile_left)/(ErrorWeight[left_idx -1 +1]-quantile_left);         // 不创建变量（塞不下了），quantile_right = ErrorWeight[left_idx -1 +1]; 


        srcs[srcidx].margin_pts[idx + batchidx*BATCH_SIZE].Q = Q;
        // srcs[srcidx].margin_pts[idx + blockDim.x].quantile = quantile;

        srcs[srcidx].margin_pts[idx + batchidx*BATCH_SIZE].position = Margin_Q(srcs[srcidx].shape,Q);
        srcs[srcidx].margin_pts[idx + batchidx*BATCH_SIZE].changed = 1;
        srcs[srcidx].margin_pts[idx + batchidx*BATCH_SIZE].Nphys = -1;
        for(int j=0;j<LENGTH-1;j++)
        {
            srcs[srcidx].margin_pts[idx + batchidx*BATCH_SIZE].deltaS[j]=0;
            srcs[srcidx].margin_pts[idx + batchidx*BATCH_SIZE].deltaS_Err[j]=0;
            srcs[srcidx].margin_pts[idx + batchidx*BATCH_SIZE].deltaS_new[j]=0;
            srcs[srcidx].margin_pts[idx + batchidx*BATCH_SIZE].Err_new[j]=0;
        }

        // 别急，下面还要连接顺序
        ParentIdx[idx] = left_idx;
    }
    __syncthreads();
    if(((idx < BATCH_SIZE)) && (srcs[srcidx].SolveSucceed==false))
    {
        // 前半截是“谁连自己”
        // 如果前一个点和自己是姊妹节点
        if((idx >= 1) && (ParentIdx[idx-1]==ParentIdx[idx]))      // 0号的前一个要取模，但负数取模有一定问题，而且反正0号线程肯定连parent
        {
            srcs[srcidx].margin_pts[idx-1 + batchidx*BATCH_SIZE].next_src_idx = idx + batchidx*BATCH_SIZE;     // 前一个点后面连自己
            srcs[srcidx].margin_pts[idx + batchidx*BATCH_SIZE].prev_src_idx = idx-1 + batchidx*BATCH_SIZE;
        }
        else
        {
            srcs[srcidx].margin_pts[left_idx].next_src_idx = idx + batchidx*BATCH_SIZE;               // parent节点后面连自己
            srcs[srcidx].margin_pts[idx + batchidx*BATCH_SIZE].prev_src_idx = left_idx;
            srcs[srcidx].margin_pts[left_idx].changed = 1;
        }
        // 后半截是“自己连谁”
        // 如果自己是姊妹节点中的老幺
        if(idx==(BATCH_SIZE-1))     // 如果是本batch中的最后一个
        {
            srcs[srcidx].margin_pts[idx + batchidx*BATCH_SIZE].next_src_idx = right_idx;                   // 自己连 right idx
            srcs[srcidx].margin_pts[right_idx].prev_src_idx = idx + batchidx*BATCH_SIZE; 
            srcs[srcidx].margin_pts[right_idx].changed = 1;
        }
        else{
            if(ParentIdx[idx] != ParentIdx[idx+1])  // 如果是普通情况
            {
                srcs[srcidx].margin_pts[idx + batchidx*BATCH_SIZE].next_src_idx = right_idx;                   // 自己连 right idx
                srcs[srcidx].margin_pts[right_idx].prev_src_idx = idx + batchidx*BATCH_SIZE;
                srcs[srcidx].margin_pts[right_idx].changed = 1;
            }
            // else{}   // 前面连好了已经
        }
    }


}
template _Global void AdaptiveLocLoop(int* srclist, src_ext_t<float>*srcs, int batchidx, bool muted);
template _Global void AdaptiveLocLoop(int* srclist, src_ext_t<double>*srcs, int batchidx, bool muted);



// 不修改g，p往下寻找
// size <<<[(NSRCS-1) // 64 +1,1,1],64>>>
template<isFloating f_T>
_Global void PhysicalTest(int* srclist, src_ext_t<f_T>*srcs, int batchidx, int prev_Ncal, bool muted)
{
    int srcidx = threadIdx.x + blockIdx.x*blockDim.x;
    if(srcidx < prev_Ncal)
    {
        srcidx = srclist[srcidx];
        // printf("(srcidx:%d)(batchidx:%d)",srcidx,batchidx);

        int Ncross = srcs[srcidx].Ncross;
        if(Ncross>0){
            int idx_cross_p;        // physical
            bool ghost_direction;     // 1 means next, 0 means previous
            int idx_cross_g;       // ghost
            complex_t<f_T> posi_p[2]; 
            complex_t<f_T> posi_g[2];
            int n_ghost;
            f_T cos2;
            f_T temp;
            int j_test_p[2];
            int temp_j;
            int iter_i;
            bool unclear;       // if cos2 > 0.1 but < 0.9
            bool test_pass;

            int idx_change[BATCH_SIZE+1];       // ghost images candidates' idx
            int j0[BATCH_SIZE+1];
            int j1[BATCH_SIZE+1];
            f_T norms[BATCH_SIZE+1];
            f_T norm_g;
            f_T& norm_p = norm_g;

            int* j_test_g = j_test_p;

            bool additional_piece;
            bool tempbool;


            for(int i=0;i<Ncross;i++)
            {      
                additional_piece = false;
                test_pass=0;
                unclear = 0;        // set to zero for each cross caustic
                // unclear_start = -1;     // "modify_len" points should be changed   与iter_i的区别是，如果unclear，iter_i++，but modify_len won't change

                idx_cross_p = srcs[srcidx].idx_cross[i];

                if(idx_cross_p<0)
                {
                    idx_cross_p = -(idx_cross_p+1);
                    ghost_direction = 0;
                    idx_cross_g = PrevSrcIdx(srcs[srcidx].margin_pts, idx_cross_p);
                }
                else
                {
                    ghost_direction = 1;
                    idx_cross_g = NextSrcIdx(srcs[srcidx].margin_pts, idx_cross_p);
                }
                if(idx_cross_g<0)
                {
                    if(!muted)
                    {
                        printf("fatal error: idx_g < 0, srcidx: %d, idx_cross_p: %d, cross_i: %d\n",srcidx,idx_cross_p,i);
                    }
                    srcs[srcidx].Break = true;
                    srcs[srcidx].SolveSucceed = 1;          // 用于跳过后续计算，会在最后一个SumArea设置为fail
                    srcs[srcidx].margin_pts[idx_cross_p].Nphys=0;
                }
                j_test_p[0] = srcs[srcidx].j_cross[i];
                j_test_p[1] = srcs[srcidx].margin_pts[idx_cross_p].next_j[j_test_p[0]];
                if(j_test_p[0]<0 || j_test_p[1]<0)
                {
                    if(!muted)
                    {
                        printf("fatal error: fake j: srcidx: %d, idx_cross_p: %d, cross_i: %d, j0: %d, j1: %d\n",srcidx,idx_cross_p,i,j_test_p[0],j_test_p[1]);
                    }
                    srcs[srcidx].Break = true;
                    srcs[srcidx].SolveSucceed = 1;          // 用于跳过后续计算，会在最后一个SumArea设置为fail
                    srcs[srcidx].margin_pts[idx_cross_p].Nphys=0;
                }       
           
                if(srcs[srcidx].margin_pts[idx_cross_p].Nphys==5)       // 这个点判断之初肯定是N5，但是有一种情况（虚像碎片连成环）会在之前的i循环里修改掉
                // if(false)
                {
                    
                    posi_p[0] = srcs[srcidx].margin_pts[idx_cross_p].images[j_test_p[0]].position;
                    posi_p[1] = srcs[srcidx].margin_pts[idx_cross_p].images[j_test_p[1]].position;
                    
                    n_ghost=0;
                    for(int j=0;j<5;j++)
                    {
                        if(!srcs[srcidx].margin_pts[idx_cross_g].images[j].physical)
                        {
                            if(n_ghost<2)
                            {
                                posi_g[n_ghost] = srcs[srcidx].margin_pts[idx_cross_g].images[j].position;
                            }
                            n_ghost++;                                
                        }
                    }
                    if(n_ghost!=2)
                    {
                        if(n_ghost==0)      // 说明这个片段之前被修改过，是由一段虚假的N5错判为N3导致的，此时错判的N3已经被修复，只要去除这个cross点就可以了。当前cross点应该是...
                                            // 所有正常的cross点都是N5，但是之前addition标记的一定是N3，所以现在添加了“addition and N5”这种情况
                        {
                            srcs[srcidx].additional_cross[i] = true;
                            srcs[srcidx].margin_pts[idx_cross_g].Nphys = 5;     // 应该是一定成立的，只是确保一下
                            continue;
                        }
                        else
                        {
                            if(!muted)
                            {
                                printf("fatal error: ghost images number wrong, srcidx: %d, batchidx: %d, idx_cross_g: %d, Nghost: %d\n",srcidx,batchidx,idx_cross_g,n_ghost);
                            }
                            srcs[srcidx].Break = true;
                            srcs[srcidx].SolveSucceed = 1;          // 用于跳过后续计算，会在最后一个SumArea设置为fail
                            #if DetailedRecord
                                srcs[srcidx].points_used = (batchidx+1) * BATCH_SIZE+1;
                            #endif                             
                        }
                       
                    }

                    iter_i = 0;
                    posi_p[0] = posi_p[1] - posi_p[0];      // vector between two points
                    posi_g[0] = posi_g[1] - posi_g[0];
                    norm_g = norm(posi_g[0]);
                    norms[iter_i] = norm(posi_p[0]);

                    if((norm_g>1e-1) ||(norms[iter_i]>1e-1))        // too far away, pass
                    {
                        test_pass = true;
                        // srcs[srcidx].margin_pts[idx_cross_p].special_Err = (norm_g+norms[iter_i])*100;        // 添加权重，强行要求在这里加点
                        // srcs[srcidx].margin_pts[idx_cross_g].special_Err = (norm_g+norms[iter_i])*100;
                        // 目的是在 ghost 和 physical 之间加点
                        if(ghost_direction)     // next is ghost
                        {
                            srcs[srcidx].margin_pts[idx_cross_p].special_Err = (norm_g+norms[iter_i])*100;      // 权重指的是 here 和 next 的区间，对应 physical 和 ghost
                        }
                        else                    // previous is ghost
                        {
                            srcs[srcidx].margin_pts[idx_cross_g].special_Err = (norm_g+norms[iter_i])*100;      // 权重指的是 here 和 next 的区间，对应 ghost 和 physical
                        }

                        continue;
                    }

                    temp = posi_p[0].re * posi_g[0].re + posi_p[0].im * posi_g[0].im;
                    cos2 = (temp*temp) / (norms[iter_i]*norm_g);


                    if(cos2<=0.1){test_pass = true;}
                    else{
                        if(cos2<0.9)
                        {
                            unclear = true;
                            // unclear_start = 0;
                        }
                    }

                    idx_change[iter_i] = idx_cross_p;    // 需要在最后修改的虚点。如果通过检测，iter_i = 0，只会修改指标<0的点，即一个都不改
                    j0[iter_i] = j_test_p[0];
                    j1[iter_i] = j_test_p[1];
                    // cos2s[iter_i] = cos2;                        

                    while((!test_pass) && (iter_i<BATCH_SIZE+1))        // 夹角不够垂直，没有jump，说明当前的physical也是ghost（假设只有N5判断会出错，N3判断严格正确）
                    {
                        iter_i++;
                        if(iter_i == BATCH_SIZE+1){break;}
                        // 是上面注释掉的if-else的整合版
                        j_test_p[1] = srcs[srcidx].margin_pts[idx_cross_p].next_j[j_test_p[1]];
                        if(ghost_direction)
                        {
                            idx_cross_p = PrevSrcIdx(srcs[srcidx].margin_pts, idx_cross_p);
                        }
                        else
                        {       
                            idx_cross_p = NextSrcIdx(srcs[srcidx].margin_pts, idx_cross_p);                  
                        }
                        if(srcs[srcidx].margin_pts[idx_cross_p].Nphys!=5)
                        {
                            // 该片段已走到尽头，也没通过
                            test_pass = true;
                            additional_piece = true;
                            break;
                        }
                        for(int j=0;j<5;j++)        // 连接：因为nextj是单向的，有一个方向需要遍历查找
                        {
                            if(srcs[srcidx].margin_pts[idx_cross_p].images[j].parity != ghost_direction)
                            {
                                temp_j = srcs[srcidx].margin_pts[idx_cross_p].next_j[j];
                                if((temp_j == j_test_p[0]))
                                {
                                    j_test_p[0] = j;
                                    break;
                                }                                    
                            }

                        }
                        posi_p[0] = srcs[srcidx].margin_pts[idx_cross_p].images[j_test_p[0]].position;
                        posi_p[1] = srcs[srcidx].margin_pts[idx_cross_p].images[j_test_p[1]].position;                        

                        posi_p[0] = posi_p[1] - posi_p[0];
                        norms[iter_i] = norm(posi_p[0]);

                        if(norms[iter_i]>1e-1)        // too far away, pass
                        {
                            test_pass = true;
                            // srcs[srcidx].margin_pts[idx_cross_p].special = 100.*norms[iter_i];   // 强行要求在附近加点

                            // 目的是在 ghost 和 physical 之间加点
                            if(ghost_direction)     // next is ghost
                            {
                                srcs[srcidx].margin_pts[idx_cross_p].special_Err = (norms[iter_i])*100;      // 权重指的是 here 和 next 的区间，对应 physical 和 ghost
                            }
                            else                    // previous is ghost
                            {
                                srcs[srcidx].margin_pts[idx_change[iter_i-1]].special_Err = (norms[iter_i])*100;      // 权重指的是 here 和 next 的区间，对应 ghost 和 physical
                            }

                            break;
                        }

                        temp = posi_p[0].re * posi_g[0].re + posi_p[0].im * posi_g[0].im;
                        cos2 = (temp*temp) / (norms[iter_i]*norm_g);

                        idx_change[iter_i] = idx_cross_p;    // 需要在最后修改的虚点
                        j0[iter_i] = j_test_p[0];
                        j1[iter_i] = j_test_p[1];
                        // cos2s[iter_i] = cos2;  


                        if(unclear)
                        {
                           // unclear 区判断标准：距离回升
                            if(norms[iter_i] > norms[iter_i-1])
                            {
                                test_pass=true;
                                // printf("unclear trigger\n");
                            }
                        }
                        else
                        {
                            if(cos2<=0.1){test_pass=true;}
                            else{
                                if(cos2<0.9)
                                {
                                    unclear=true;
                                    // unclear_start = iter_i;       // +1 因为iter_i++放在后面了
                                }
                            }
                            // modify_len++;
                        }
                    }


                    if((iter_i != BATCH_SIZE+1) && (norms[iter_i]<=1e-1))    // 成功了才改，不成功反向检查
                    {
                        if(!additional_piece)
                        {
                            // connect modified points
                            if(!unclear)
                            {
                                srcs[srcidx].margin_pts[idx_cross_p].next_idx[j_test_p[0]] = idx_cross_p;
                                srcs[srcidx].margin_pts[idx_cross_p].next_j[j_test_p[0]] = j_test_p[1];

                                if(ghost_direction)
                                {
                                    srcs[srcidx].idx_cross[i] = idx_cross_p;
                                }
                                else{
                                    srcs[srcidx].idx_cross[i] = -idx_cross_p-1;
                                }                                

                                srcs[srcidx].additional_cross[i] = false;
                                
                                srcs[srcidx].j_cross[i] = j_test_p[0];
                            }
                            else
                            {
                                if(iter_i>1)
                                {                                 
                                    srcs[srcidx].margin_pts[idx_change[iter_i-1]].next_idx[j0[iter_i-1]] = idx_change[iter_i-1];
                                    srcs[srcidx].margin_pts[idx_change[iter_i-1]].next_j[j0[iter_i-1]] = j_test_p[1];
                                    if(ghost_direction)
                                    {
                                        srcs[srcidx].idx_cross[i] = idx_change[iter_i-1];
                                    }
                                    else{
                                        srcs[srcidx].idx_cross[i] = -idx_change[iter_i-1]-1;
                                    }
                                    srcs[srcidx].j_cross[i] = j0[iter_i-1]; 
                                    srcs[srcidx].additional_cross[i] = false;                                   
                                }
                                // else{}       // else, 说明第一个点就unclear，查询下一个点发现没问题，所以谁也不修改
                                // iter_i = 0: 有unclear至少为1
                                // iter_i = 1: 一开始进了unclear，检测一次发现下一个更实，所以0号是实的
                                
                            }                            
                        }
                        else        // additional piece = true
                        {
                            // 需要帮人家再补起来
                            // 就是，进入broken之前，连接出broken之后
                            // 

                            // 找到连接开头的那个点
                            int prev_pt_idx = -1;
                            int prev_pt_j = -1;
                            bool found_prev = false;

                            // 试试是不是自己连自己
                            for(int prevj=0;prevj<5;prevj++)
                            {
                                prev_pt_idx = idx_change[0];
                                if((srcs[srcidx].margin_pts[prev_pt_idx].next_idx[prevj] == prev_pt_idx) && (srcs[srcidx].margin_pts[prev_pt_idx].next_j[prevj] == j0[0]))
                                {
                                    prev_pt_j = prevj;
                                    found_prev = true;
                                }
                            }
                            // 没找到，说明是前后连过来的，取决于 ghost_direction (这玩意好像就是手征啊，我存这么个莫名其妙的东西干啥？)
                            if(!found_prev)
                            {
                                for(int prevj=0;prevj<5;prevj++)
                                {
                                    prev_pt_idx = PrevSrcIdx(srcs[srcidx].margin_pts, idx_change[0]);
                                    if((srcs[srcidx].margin_pts[prev_pt_idx].next_idx[prevj] == prev_pt_idx) && (srcs[srcidx].margin_pts[prev_pt_idx].next_j[prevj] == j0[0]))
                                    {
                                        prev_pt_j = prevj;
                                        found_prev = true;
                                    }
                                }
                            }
                            if(!found_prev)
                            {
                                for(int prevj=0;prevj<5;prevj++)
                                {
                                    prev_pt_idx = NextSrcIdx(srcs[srcidx].margin_pts, idx_change[0]);
                                    if((srcs[srcidx].margin_pts[prev_pt_idx].next_idx[prevj] == prev_pt_idx) && (srcs[srcidx].margin_pts[prev_pt_idx].next_j[prevj] == j0[0]))
                                    {
                                        prev_pt_j = prevj;
                                        found_prev = true;
                                    }
                                }
                            }

                            srcs[srcidx].margin_pts[prev_pt_idx].next_idx[prev_pt_j] = srcs[srcidx].margin_pts[idx_change[0]].next_idx[j1[0]];
                            srcs[srcidx].margin_pts[prev_pt_idx].next_j[prev_pt_j] = srcs[srcidx].margin_pts[idx_change[0]].next_j[j1[0]];


                            srcs[srcidx].additional_cross[i] = true;
                            srcs[srcidx].margin_pts[idx_change[0]].Nphys = 3;       // 实际上下面也有，只是强调一下，增加鲁棒性
                        }



                        // printf("%d,%d,%d,iters: %d, bool: %d, cos2: %.16f, j0: %d, j1: %d\n",idx_cross_p,idx_cross_g,batchidx,iter_i,ghost_direction,cos2,j_test[0],j_test[1]);

                        // 修改错误的信息
                        // 如果除非是unclear段的成功位置，那里旧点也是实像
                        // 其他情况旧点都是虚像
                        for(int iii=0;iii<iter_i - unclear;iii++)       // unclear==0: iter_i 是实像，故修改0～iter_i-1; unclear: 再少一个
                        {
                            // if(srcidx==16 && batchidx==8)
                            // {
                            //     printf("i: %d, idx = %d, j0 = %d, j1 = %d\n",i,idx_change[iii],j0[iii],j1[iii]);
                            // }
                            // change the connection of wrong points
                            // j_test[0] is always the chain to connect another
                            srcs[srcidx].margin_pts[idx_change[iii]].images[j0[iii]].physical = false;
                            srcs[srcidx].margin_pts[idx_change[iii]].images[j1[iii]].physical = false;
                            srcs[srcidx].margin_pts[idx_change[iii]].Nphys = 3;
                            srcs[srcidx].margin_pts[idx_change[iii]].next_idx[j0[iii]] = -2;
                            srcs[srcidx].margin_pts[idx_change[iii]].next_idx[j1[iii]] = -2;

                            srcs[srcidx].margin_pts[idx_change[iii]].deltaS[j0[iii]] = 0;
                            srcs[srcidx].margin_pts[idx_change[iii]].deltaS[j1[iii]] = 0;
                            srcs[srcidx].margin_pts[idx_change[iii]].deltaS_Err[j0[iii]] = 0;
                            srcs[srcidx].margin_pts[idx_change[iii]].deltaS_Err[j1[iii]] = 0;

                            srcs[srcidx].margin_pts[idx_change[iii]].deltaS_new[j0[iii]] = 0;
                            srcs[srcidx].margin_pts[idx_change[iii]].deltaS_new[j1[iii]] = 0;
                            srcs[srcidx].margin_pts[idx_change[iii]].Err_new[j0[iii]] = 0;
                            srcs[srcidx].margin_pts[idx_change[iii]].Err_new[j1[iii]] = 0;
                            // if constexpr (DetailRecord)
                            #if DetailedRecord
                            // {
                                srcs[srcidx].margin_pts[idx_change[iii]].E1[j0[iii]] = 0;
                                srcs[srcidx].margin_pts[idx_change[iii]].E1[j1[iii]] = 0;
                                srcs[srcidx].margin_pts[idx_change[iii]].E2[j0[iii]] = 0;
                                srcs[srcidx].margin_pts[idx_change[iii]].E2[j1[iii]] = 0;
                                srcs[srcidx].margin_pts[idx_change[iii]].E3[j0[iii]] = 0;
                                srcs[srcidx].margin_pts[idx_change[iii]].E3[j1[iii]] = 0;                                
                            // }
                            #endif


                            // if(srcidx==13)
                            // {
                            //     printf("idx to N3: %d\n",idx_change[iii]);
                            // }                                                                                  
                        }
                    }
                    else      // 认为只有虚像N3错判为实N5失败了，实际上有实像N5错判为虚像N3
                    {

                        // if(srcidx==280)
                        // {
                        //     printf("here! idx_cross_p: %d\n",idx_change[0]);
                        // }
                        // 标记方法为：idx_change 系列数组，0号还是srcs[srcidx].idx_cross[i]，从1号起是反向的N3点，注意j0和j1分别连成链
                        test_pass = false;
                        iter_i = 1;
                        idx_cross_p = idx_change[0];        // 前面拿过了
                        posi_p[0] = srcs[srcidx].margin_pts[idx_cross_p].images[j0[0]].position;
                        posi_p[1] = srcs[srcidx].margin_pts[idx_cross_p].images[j1[0]].position;
                                                

                        idx_cross_g = (ghost_direction) ? NextSrcIdx(srcs[srcidx].margin_pts, idx_cross_p) \
                                                        : PrevSrcIdx(srcs[srcidx].margin_pts, idx_cross_p) ;
                        idx_change[iter_i] = idx_cross_g;

                        n_ghost=0;
                        for(int j=0;j<5;j++)
                        {
                            if((!srcs[srcidx].margin_pts[idx_cross_g].images[j].physical) && (n_ghost<2))
                            {
                                posi_g[n_ghost] = srcs[srcidx].margin_pts[idx_cross_g].images[j].position;
                                j_test_g[n_ghost] = j;
                                n_ghost++;                                
                            }
                        }           // 这个不用检查了，实际上前半截已经算过了，只是没存，该报的错都报了

                        posi_p[0] = posi_p[1] - posi_p[0];      // vector between two points
                        posi_g[0] = posi_g[1] - posi_g[0];
                        norm_p = norm(posi_p[0]);
                        norms[iter_i] = norm(posi_g[0]);
                        temp = posi_p[0].re * posi_g[0].re + posi_p[0].im * posi_g[0].im;
                        // if(temp>0)
                        // {
                        //     j0[iter_i] = j_test_g[0]; j1[iter_i] = j_test_g[1];
                        // }
                        // else
                        // {
                        //     j0[iter_i] = j_test_g[1]; j1[iter_i] = j_test_g[0];
                        // }
                        tempbool = (temp>0);
                        j0[iter_i] = j_test_g[int(!tempbool)];
                        j1[iter_i] = j_test_g[int( tempbool)];

                        cos2 = (temp*temp) / (norms[iter_i-1]*norms[iter_i]);
                        if(cos2<=0.1){test_pass = true;}
                        else{
                            if(cos2<0.9)
                            {
                                unclear = true;
                                // unclear_start = 0;
                            }
                        }     
                        while((!test_pass) && (iter_i<BATCH_SIZE+1-1))
                        {
                            // 只改posi_g
                            iter_i++;
                            if(iter_i == BATCH_SIZE+1){break;}
                            idx_cross_p = idx_cross_g;
                            idx_cross_g = (ghost_direction) ? NextSrcIdx(srcs[srcidx].margin_pts, idx_cross_p) \
                                                            : PrevSrcIdx(srcs[srcidx].margin_pts, idx_cross_p) ;
                            idx_change[iter_i] = idx_cross_g;


                            if(srcs[srcidx].margin_pts[idx_cross_g].Nphys!=3)
                            {
                                // 该片段已走到尽头，也没通过
                                test_pass = true;
                                additional_piece = true;

                                break;
                            }



                            n_ghost=0;
                            for(int j=0;j<5;j++)
                            {
                                if((!srcs[srcidx].margin_pts[idx_cross_g].images[j].physical) && (n_ghost<2))
                                {
                                    posi_g[n_ghost] = srcs[srcidx].margin_pts[idx_cross_g].images[j].position;
                                    j_test_g[n_ghost] = j;
                                    n_ghost++;                                
                                }
                            }
                            posi_g[0] = posi_g[1] - posi_g[0];
                            norms[iter_i] = norm(posi_g[0]);
                            if(n_ghost!=2)
                            {
                                // if(n_ghost==0)      // 到头了！
                                // {
                                //     // 实际上应该在前面  if(srcs[srcidx].margin_pts[idx_cross_g].Nphys!=3) 就触发了，不应该
                                // }
                                // else
                                // {
                                if(!muted)
                                {
                                    printf("fatal error: ghost images number wrong, srcidx: %d, batchidx: %d, idx_cross_g: %d, Nghost: %d\n",srcidx,batchidx,idx_cross_g,n_ghost);
                                }
                                srcs[srcidx].Break = true;
                                srcs[srcidx].SolveSucceed = 1;          // 用于跳过后续计算，会在最后一个SumArea设置为fail
                                #if DetailedRecord
                                    srcs[srcidx].points_used = (batchidx+1) * BATCH_SIZE+1;
                                #endif                                    
                                // }
                            }


                            if(norms[iter_i]>1e-1)        // too far away, pass
                            {
                                test_pass = true;
                                // srcs[srcidx].margin_pts[idx_cross_p].special = 100.*norms[iter_i];   // 强行要求在附近加点

                                // 目的是在 ghost 和 physical 之间加点
                                if(ghost_direction)     // next is ghost
                                {
                                    srcs[srcidx].margin_pts[idx_change[iter_i-1]].special_Err = (norms[iter_i])*100;      // 权重指的是 here 和 next 的区间，对应 physical 和 ghost
                                }
                                else                    // previous is ghost
                                {
                                    srcs[srcidx].margin_pts[idx_change[ iter_i ]].special_Err = (norms[iter_i])*100;      // 权重指的是 here 和 next 的区间，对应 ghost 和 physical
                                }
                                break;
                            }

                            temp = posi_p[0].re * posi_g[0].re + posi_p[0].im * posi_g[0].im;
                            cos2 = (temp*temp) / (norms[iter_i]*norm_p);
                            // printf("iter_i: %d, cos2: %.12f\n",iter_i,cos2);

                            tempbool = (temp>0);
                            j0[iter_i] = j_test_g[int(!tempbool)];
                            j1[iter_i] = j_test_g[int( tempbool)];    

                            if(unclear)
                            {
                            // unclear 区判断标准：距离回升
                                if(norms[iter_i] > norms[iter_i-1])
                                {
                                    test_pass=true;
                                    // printf("unclear trigger\n");
                                }
                            }
                            else
                            {
                                if(cos2<=0.1){test_pass=true;}
                                else{
                                    if(cos2<0.9)
                                    {
                                        unclear=true;
                                    }
                                }
                            }

                        }       // end while

                        // printf("iter_i: %d\n",iter_i);
                        if((iter_i != BATCH_SIZE+1) && (norms[iter_i]<=1e-1))    // 成功了才改，还不成功就报错
                        {
                            // if(srcidx==654 && batchidx==12)
                            // {
                            //     printf("i: %d, idx_cross_marked: %d\n",i,srcs[srcidx].idx_cross[i]);
                            // }

                            // 修改错误的信息
                            // 如果除非是unclear段的成功位置，那里旧点也是实像
                            // 其他情况旧点都是虚像
                            tempbool = srcs[srcidx].margin_pts[idx_change[0]].images[j0[0]].parity;
                            srcs[srcidx].margin_pts[idx_change[0]].next_idx[j0[0]] = idx_change[0+1];
                            srcs[srcidx].margin_pts[idx_change[0]].next_j[j0[0]]   = j0[0+1];

                            for(int iii=1;iii<iter_i;iii++)       // unclear==0: iter_i 是实像，故修改0～iter_i-1; unclear: 再少一个
                            {
                                // if(srcidx==16 && batchidx==8)
                                // {
                                //     printf("i: %d, idx = %d, j0 = %d, j1 = %d\n",i,idx_change[iii],j0[iii],j1[iii]);
                                // }
                                // change the connection of wrong points
                                // j_test[0] is always the chain to connect another
                                srcs[srcidx].margin_pts[idx_change[iii]].images[j0[iii]].physical = true;
                                srcs[srcidx].margin_pts[idx_change[iii]].images[j1[iii]].physical = true;
                                srcs[srcidx].margin_pts[idx_change[iii]].Nphys = 5;
                                srcs[srcidx].margin_pts[idx_change[iii]].images[j0[iii]].parity =  tempbool;
                                srcs[srcidx].margin_pts[idx_change[iii]].images[j1[iii]].parity = !tempbool;

                                // ghost = 1: idx_change along next; parity = 0: 0 chain connect to next
                                srcs[srcidx].margin_pts[idx_change[iii]].next_idx[j0[iii]] = idx_change[iii+1];
                                srcs[srcidx].margin_pts[idx_change[iii]].next_j[j0[iii]]   = j0[iii+1];
                                // ghost = 1: idx_change along next; parity = 1: 1 chain connect to next
                                srcs[srcidx].margin_pts[idx_change[iii]].next_idx[j1[iii]] = idx_change[iii-1];
                                srcs[srcidx].margin_pts[idx_change[iii]].next_j[j1[iii]]   = j1[iii-1];

                                                                               
                            }

                            if(!additional_piece)
                            {
                                // connect modified points

                                
                                // srcs[srcidx].margin_pts[idx_cross_p].next_idx[j_test_p[0]] = idx_cross_p;
                                // srcs[srcidx].margin_pts[idx_cross_p].next_j[j_test_p[0]] = j_test_p[1];

                                // if(ghost_direction)
                                // {
                                //     srcs[srcidx].idx_cross[i] = idx_cross_p;
                                // }
                                // else{
                                //     srcs[srcidx].idx_cross[i] = -idx_cross_p-1;
                                // }                                

                                // srcs[srcidx].additional_cross[i] = false;
                                
                                // srcs[srcidx].j_cross[i] = j_test_p[0];

                                // 自己封口

                                srcs[srcidx].margin_pts[idx_change[iter_i-1]].next_idx[j0[iter_i-1]] = idx_change[iter_i-1];
                                srcs[srcidx].margin_pts[idx_change[iter_i-1]].next_j[j0[iter_i-1]] = j1[iter_i-1];

                                srcs[srcidx].idx_cross[i] = (ghost_direction)   ? idx_change[iter_i-1]    \
                                                                                : -idx_change[iter_i-1]-1 ;
                                srcs[srcidx].additional_cross[i] = false;
                                srcs[srcidx].j_cross[i] = j0[iter_i-1];
                        
                            }
                            else        // additional piece = true
                            {
                                // 把两端的封口重新打开，想办法连上
                                // 重新进行一次简化的findnext？但findnext不是以srcpt为单位进行的...

                                // if(srcidx==654&&batchidx==12)
                                // {
                                //     printf("idx_change[iter_i-1]: %d, idx_change[iter_i]: %d\n",idx_change[iter_i-1], idx_change[iter_i]);
                                // }

                                // 示例：
                                // idx:             758, 759, 760, 761, 762
                                // real N:            5,   5,   5,   5,   3
                                // judged N:          5,   3,   5,   5,   3
                                // idx_change[iter_i-1] == 759, idx_change[iter_i] == 760, additional piece

                                for(int jjjj=0;jjjj<5;jjjj++)
                                {
                                    if(srcs[srcidx].margin_pts[idx_change[iter_i]].next_idx[jjjj]==idx_change[iter_i])      // 如果 760 的某个 j 连接自己
                                    {
                                        if(srcs[srcidx].margin_pts[idx_change[iter_i]].images[jjjj].parity == srcs[srcidx].margin_pts[idx_change[iter_i-1]].images[j1[iter_i-1]].parity)       // 如果 760 发出自连的这个 j 手征和 759 被连接的还是一样的
                                        {
                                            // 就是它了（仅限2体）
                                            // 被 jjjj 连的像现在被 759 的 j0 连
                                            srcs[srcidx].margin_pts[idx_change[iter_i-1]].next_idx[j0[iter_i-1]] = idx_change[iter_i];
                                            srcs[srcidx].margin_pts[idx_change[iter_i-1]].next_j[j0[iter_i-1]] = srcs[srcidx].margin_pts[idx_change[iter_i]].next_j[jjjj];
                                            // 760 发起自连的像现在连 759 的 j1
                                            srcs[srcidx].margin_pts[idx_change[iter_i]].next_idx[jjjj] = idx_change[iter_i-1];
                                            srcs[srcidx].margin_pts[idx_change[iter_i]].next_j[jjjj] = j1[iter_i];
                                            break;

                                        }
                                    }
                                    
                                }

                                srcs[srcidx].additional_cross[i] = true;


                            }

                        }

                    }

                // printf("%d,%d\n",idx_cross_p,idx_cross_g);
                }
                else{       // 这个片段在前面就被改掉了，现在是N3，整个小片段都是虚的
                    // if(!muted)
                    // {
                    //     // printf("Waring: ghost cross: srcidx: %d, batchidx: %d, idx_cross_p: %d\n",srcidx,batchidx,idx_cross_p);
                    //     // srcs[srcidx].idx_cross[i] = (1<<31);        // int min, to mark additional case
                    // }
                    // if(srcidx==280&&batchidx==1)
                    // {
                    //     printf("N3,i: %d, iter_i: %d, unclear: %d\n",i,iter_i,unclear);
                    // }
                    srcs[srcidx].additional_cross[i] = true;
                }
            }
        }
    }
}
template _Global void PhysicalTest(int* srclist, src_ext_t<float>*srcs, int batchidx, int prev_Ncal, bool muted);
template _Global void PhysicalTest(int* srclist, src_ext_t<double>*srcs, int batchidx, int prev_Ncal, bool muted);


// size: <<<(Ncal-1) / 64+1 , min(64,Ncal)>>>
template<isFloating f_T>
_Global void SuccessCheck(int* Ncal, int* srclist, src_ext_t<f_T>* srcs, int prev_Ncal, bool muted)
{
    int srcidx = threadIdx.x + blockIdx.x * blockDim.x;
    int temp = -1;
    if(srcidx < prev_Ncal)
    {

        srcidx = srclist[srcidx];
        if((srcs[srcidx].SolveSucceed==false)&&(srcs[srcidx].Break==false))
        {
            temp = atomicAdd(Ncal,1);
        }
    }
    if(temp>=0)
    {
        srclist[temp] = srcidx;
    }
    // printf("%d,%d,%d,%d\n",srcidx,temp,srcs[srcidx].SolveSucceed,srcs[srcidx].Break);
}
template _Global void SuccessCheck(int* Ncal, int* srclist, src_ext_t<float>* srcs, int prev_Ncal, bool muted);
template _Global void SuccessCheck(int* Ncal, int* srclist, src_ext_t<double>* srcs, int prev_Ncal, bool muted);

