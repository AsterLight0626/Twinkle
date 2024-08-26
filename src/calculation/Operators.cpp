#include"Operators.h"

// 使用方法：image_pt_t.Physical = isPhysical(image_pt_t.position,params,src_pt_t.position)
_Device bool isPhysical(const complex_t<double>& lambda, const complex_t<double>* params, const complex_t<double>& loc_src, double& physical_err)
{
	double error;
	complex_t<double> temp;

	temp = f_zbar(lambda,params) + lambda - loc_src;	// temp = z + f(zbar) - zeta
	error = abs_c(temp);
	physical_err = error;
	// double的阈值！
	if(error < THRESHOLD_REAL_DOUBLE){
		return true;
	}
	else{
		return false;
	}
}
_Device bool isPhysical(const complex_t<float>& lambda, const complex_t<float>* params, const complex_t<float>& loc_src, float& physical_err)
{
	float error;
	complex_t<float> temp;

	temp = f_zbar(lambda,params) + lambda - loc_src;	// temp = z + f(zbar) - zeta
	error = abs_c(temp);
	physical_err = error;
	// float的阈值！
	if(error < THRESHOLD_REAL_FLOAT){
		return true;
	}
	else{
		return false;
	}
}


template<isFloating f_T>
_Device bool isPhysicalAll(image_pt_t<f_T>* imgs, const f_T s, const f_T m1, const f_T m2, const complex_t<f_T>& loc_src, int* Nphys, f_T* errors_out)
{
	f_T errors[LENGTH-1];
	complex_t<f_T> temp;
	int rank[LENGTH-1];
	int count;		// used for Nphys == -1
	// int& evencount = count;
	// int oddcount;
	bool fail_test = false;
	// int j_left[2];
	// int j_left_idx=0;
	

	for(int i=0;i<LENGTH-1;i++)
	{
		// temp = f_zbar(imgs[i].position,params) + imgs[i].position - loc_src;	// temp = z + f(zbar) - zeta
		temp = f_zbar(imgs[i].position,s,m1,m2) + imgs[i].position - loc_src;	// temp = z + f(zbar) - zeta
		errors[i] = abs_c(temp);
		rank[i] = i;
	}
	
	BubbleArgSort(rank,errors,LENGTH-1);			// errors also be sorted


	if((*Nphys)==-1)		// 前面没做判断
	{
		for(int j=0;j<5;j++)
		{
			imgs[j].physical = true;		// 先全都是 true，查到false再改
		}

		if(errors[3]*dlmin > (errors[2] + 1e-12))		// same criterion as VBB, "1e-12" is for double
		// if(errors[3] > 1e-4)
		// 在很远离焦散线的地方，由于数值不稳定，有时残差会相差比较大（甚至达到1e-7），在这种安全区域可以用 jump 来判断，效果也比较好
		{
			imgs[rank[4]].physical = false;
			imgs[rank[3]].physical = false;

			// if(rank[4]<0 || rank[4]>4 || rank[3]<0 || rank[3]>4){printf("strange mistake\n");}
			*Nphys = 3;
		}
		else
		{
			count = 0;
			for(int j=4;j>0;j--)
			{
				// if(img_pts[rank[j]].parity!=img_pts[rank[j-1]].parity)
				// {

				if( errors[j]>2e-12 && errors[j-1]>2e-12 && fabs(errors[j] / errors[j-1] - 1)<2e-15 )
				{
					imgs[rank[j]].physical = false;
					imgs[rank[j-1]].physical = false;
					break;
				}
				// }
				count++;
			}
			*Nphys = (count==4) ? 5 : 3;
			// if(count==4)
			// {
			// 	*Nphys = 5;
			// }
			// else
			// {
			// 	*Nphys = 3;
			// }			
		}


	}


	// if((*Nphys)==3)
	// {
	// 	evencount = 1;
	// 	oddcount = 2;
	// 	for(int j=0;j<LENGTH-1;j++)
	// 	{
	// 		if(imgs[rank[j]].parity==0)
	// 		{
	// 			if(evencount>0)
	// 			{
	// 				evencount -= 1;
	// 				imgs[rank[j]].physical = true;
	// 			}
	// 			else{
	// 				imgs[rank[j]].physical = false;
	// 				j_left[j_left_idx] = j;
	// 				j_left_idx++;
	// 			}
	// 		}
	// 		else
	// 		{
	// 			if(oddcount>0)
	// 			{
	// 				oddcount -= 1;
	// 				imgs[rank[j]].physical = true;
	// 			}
	// 			else{
	// 				imgs[rank[j]].physical = false;
	// 				j_left[j_left_idx] = j;
	// 				j_left_idx++;
	// 			}
	// 		}
	// 	}
	// 	if(   abs((errors[j_left[0]]+1e-12) / (errors[j_left[1]]+1e-12) - 1) > 1e-4  )
	// 	{
	// 		fail_test = true;
	// 	}
	// }

	// if((*Nphys)==5)
	// // 5 good roots
	// {
	// 	// if(errors[3] > 1e-4)	// 如果误差比1e-2都大，那不管碰撞检验也要大的删掉，反正后面还有test	（231218, test 在哪呢？）（Kuang's method)
	// 	// if(errors[3]*dlmin > (errors[2] + 1e-12))		// same criterion as VBB, "1e-12" is for double
	// 	if(false)
	// 	{
	// 		// // 修改Nphys
	// 		// *Nphys = 3;
	// 		// // printf("collision fail\n");
	// 		// for(int j=0; j<3;j++)
	// 		// {
	// 		// 	physical[rank[j]] = true;
	// 		// }
	// 		// for(int j=3;j<LENGTH-1;j++)
	// 		// {
	// 		// 	physical[rank[j]] = false;
	// 		// }

	// 	}
	// 	else{
	// 		for(int j=0;j<LENGTH-1;j++)
	// 		{
	// 			imgs[j].physical = true;
	// 		}			
	// 	}

	// }

	// if((*Nphys)!=5 && (*Nphys)!=3 && (*Nphys)!=-1){printf("Nphys wrong! %d\n",*Nphys);}

	for(int i=0;i<5;i++)
	{
		errors_out[rank[i]] = errors[i];
	}


	return fail_test;

}
template _Device bool isPhysicalAll(image_pt_t<float>* imgs,  const float s, const float m1, const float m2, const complex_t<float>& loc_src, int* Nphys, float* errors_out);
template _Device bool isPhysicalAll(image_pt_t<double>* imgs,  const double s, const double m1, const double m2, const complex_t<double>& loc_src, int* Nphys, double* errors_out);


// 从小到大
template<isFloating f_T>
_Device void BubbleArgSort(int* index, f_T* arr, const int len)
// 从小到大
{
	f_T temp;
	int temp_idx;
	int i, j;
	for (i = 0; i < len - 1; i++)
	{
		for (j = 0; j < len - 1 - i; j++)
		{
			if (arr[j] > arr[j + 1])
			{
				// swap(arr[j], arr[j + 1]);
				temp = arr[j];
				arr[j] = arr[j+1];
				arr[j+1] = temp;
				temp_idx = index[j];
				index[j] = index[j+1];
				index[j+1] = temp_idx;
			}		
		}
	}
}
template _Device void BubbleArgSort(int* index, float* arr, const int len);
template _Device void BubbleArgSort(int* index, double* arr, const int len);

template <isFloating f_T>
_Device f_T DeltaS1(const complex_t<f_T>& z_here,const complex_t<f_T>& z_next)
{
	f_T deltaS;
	deltaS = -((z_next.im + z_here.im) * (z_next.re - z_here.re) - (z_next.re + z_here.re) * (z_next.im - z_here.im)) / 4;
	// deltaS = (z_here.re*z_next.im - z_here.im*z_next.re)/2;
	// -((x2_k + x2_k_1) * (x1_k - x1_k_1) - (x1_k + x1_k_1) * (x2_k - x2_k_1)) / 4;
	return deltaS;
}
template _Device float DeltaS1(const complex_t<float>& z_here,const complex_t<float>& z_next);
template _Device double DeltaS1(const complex_t<double>& z_here,const complex_t<double>& z_next);

template<isFloating f_T>
_Device int BisearchLeft(const f_T* values_in_order, const f_T to_insert, const int len)
{
	// values是第一个点（不是第零个点！）到最后一个点的取值
	// 如果插在第一个点左侧，即0到values[0]之间
	if(to_insert < values_in_order[0])
	{
		return 0;
	}
	else
	{
		int left = 0;;
		int right = len-1;			// values[len-1] = 1; to_insert < 1, 一定不会出界
		int compare_now = len/2;	// 中间开始找，反正也用不了几步

		while( (right-left)>1)
		{
			if(to_insert > values_in_order[compare_now])
			{
				left = compare_now;
			}
			else
			{
				right = compare_now;
			}
			compare_now = (left + right)/2;
		}
		// 上面的left和right是value的指标，但是value[0]对应pt[1]的分度值，需要偏移一点儿
		return left + 1;
	}
}
template _Device int BisearchLeft(const float* values_in_order, const float to_insert, const int len);
template _Device int BisearchLeft(const double* values_in_order, const double to_insert, const int len);

// q 是01之间的参数，基本是归一化的弧长参数
_Device complex_t<float> Margin_Q(const src_shape_t<float>& shape, const float& q)
{
    float theta = q*2*PI;
    // float x = __cosf(theta) * shape.rho;
    // float y = __sinf(theta) * shape.rho;
    float x = std::cos(theta) * shape.rho;
    float y = std::sin(theta) * shape.rho;
    return complex_t<float>(x,y)+shape.loc_centre;
}
_Device complex_t<double> Margin_Q(const src_shape_t<double>& shape,const double& q)
{
    double theta = q*2*PI;
    double x = std::cos(theta) * shape.rho;
    double y = std::sin(theta) * shape.rho;
    return complex_t<double>(x,y)+shape.loc_centre;
}


template<isFloating f_T>
_Device f_T Area_src(const src_shape_t<f_T>& shape)
{
    return PI * shape.rho * shape.rho;
}
template _Device float Area_src(const src_shape_t<float>& shape);
template _Device double Area_src(const src_shape_t<double>& shape);


_Device void PathLoc(const src_params_t<float>& move_params, const float time, complex_t<float>& loc)
{
    float distance;
    complex_t<float> PerpendicularFoot;
    // complex_t<float> loc;

    // PerpendicularFoot.re = __cosf(move_params.alpha - PI/2) * move_params.u_0;
    // PerpendicularFoot.im = __sinf(move_params.alpha - PI/2) * move_params.u_0;
    PerpendicularFoot.re = std::cos(move_params.alpha - PI/2) * move_params.u_0;
    PerpendicularFoot.im = std::sin(move_params.alpha - PI/2) * move_params.u_0;


    distance = (time - move_params.t_0)/move_params.t_E;

    loc = PerpendicularFoot + complex_t<float>(std::cos(move_params.alpha),std::sin(move_params.alpha)) * distance;
    // return loc;
}
_Device void PathLoc(const src_params_t<double>& move_params, const double time, complex_t<double>& loc)
{
    double distance;
    complex_t<double> PerpendicularFoot;
    // complex_t<double> loc;

    PerpendicularFoot.re = std::cos(move_params.alpha - PI/2) * move_params.u_0;
    PerpendicularFoot.im = std::sin(move_params.alpha - PI/2) * move_params.u_0;

    distance = (time - move_params.t_0)/move_params.t_E;

    // loc = PerpendicularFoot + complex_t<double>(std::cos(move_params.alpha),std::sin(move_params.alpha)) * distance;
	loc.re = std::cos(move_params.alpha) * distance;
	loc.im = std::sin(move_params.alpha) * distance;
	loc = loc + PerpendicularFoot;
    // return 0;
}


template<isFloating f_T>
_Device bool PolishPP(image_pt_t<f_T>* images, const f_T* jacobi, bool muted)
{
	int neat_parity=0;
	// int ref_neat_parity= -NBODY+1;
	int real_images=0;
	bool fail=false;
	
	for(int j=0;j<LENGTH-1;j++)
	{
		real_images += images[j].physical;
		if(images[j].physical)
		{
			neat_parity += (images[j].parity) ? -1 : 1;
		}
	}
	// 如果实虚像数量不对，直接抛弃
	if(real_images!=3 && real_images!=5)
	{
		if(!muted)
		{
			printf("N real images wrong, %d\n",real_images);
		}
		return true;
	}


	// 不等于就麻烦了！
	if(neat_parity != -NBODY+1)
	{
		int order_p[LENGTH-1];
		f_T temp[LENGTH-1];
		for(int j=0;j<LENGTH-1;j++)
		{
			order_p[j] = j;
			temp[j] = 1/jacobi[j];
		}
		BubbleArgSort(order_p,temp,LENGTH-1);


		// order_p[i]是第i小的Mag（包括正负号）对应的指标


		if(abs(-NBODY+1-neat_parity)>3)		// 0，2，-2是允许的，更多不行，例如4，-4
		{
			// printf("ERROR! Too many parity mistakes\n");
			fail=true;
		}
		else
		{
			// 如果正手征多了，说明有的负手征跨过散焦线变成了正手征，此时放大率最大的嫌疑最大
			if(neat_parity > -NBODY+1)
			{
				for(int jj=LENGTH-1-1;jj>=0;jj--)
				{
					if(images[order_p[jj]].physical)
					{
						images[order_p[jj]].parity = true;
						break;
					}
				}
			}
			// 如果负手征多了，说明有的正手征跨过散焦线变成了负手征，此时放大率最小的嫌疑最大
			else{
				for(int jj=0;jj<=LENGTH-1-1;jj++)
				{
					if(images[order_p[jj]].physical)
					{
						images[order_p[jj]].parity = false;
						break;
					}					
				}
			}			
		}


	}
	return fail;
}
template _Device bool PolishPP(image_pt_t<float>* images, const float* jacobi, bool muted);
template _Device bool PolishPP(image_pt_t<double>* images, const double* jacobi, bool muted);


// 没乘 rho^2
template<isFloating f_T>
_Device f_T Mu_QC(const complex_t<f_T>* params, const complex_t<f_T>& z)
{
	complex_t<f_T> items[NBODY];
	complex_t<f_T> f1,f2,f3;
	f_T J;
	f_T& muQC = items[0].re;

	f1=0;
	for(int i=0;i<NBODY;i++)
	{
		items[i] = params[NBODY + i].re / (z-params[i]) / (z-params[i]);
		f1 = f1 + items[i];
	}
	J = 1-(f1 * conj(f1)).re;

	f2=0;
	for(int i=0;i<NBODY;i++)
	{
		// items[i] /= (z-params[i]);
		// items[i] *= -2;
		items[i] = -f_T(2.) * items[i] / (z-params[i]);
		f2 = f2 + items[i];
	}

	f3=0;
	for(int i=0;i<NBODY;i++)
	{
		// items[i] /= (z-params[i]);
		// items[i] *= -3;
		items[i] = -f_T(3.) * items[i] / (z-params[i]);
		f3 = f3 + items[i];
	}

	f1 = conj(f1);
	muQC = 1./(J*J*J*J*J) * (fabs(-2.*(f_T(3.)* f1*f1*f1 * f2*f2  -  (f_T(3.)-f_T(3.)*J + J*J/f_T(2.))* norm(f2)  +  J*f1*f1 * f3).re) + fabs(6.*(f_T(3.)* f1*f1*f1 * f2*f2).im));

	return muQC;
}
template _Device float Mu_QC(const complex_t<float>* params, const complex_t<float>& z);
template _Device double Mu_QC(const complex_t<double>* params, const complex_t<double>& z);


// 命名逻辑同VBB该模块
template<isFloating f_T>
_Device f_T GhostTst(const complex_t<f_T>* params, const complex_t<f_T>& z_G, const complex_t<f_T>& z_hat)
{
	complex_t<f_T> fz1,fz2,fz_hat1,J_hat;
	f_T J;		// J(z_G)
	complex_t<f_T> items[NBODY];
	complex_t<f_T>& temp = items[0];	// the inverse of LHS of Eq(47) in VBB2

	fz1=0;
	for(int i=0;i<NBODY;i++)
	{
		items[i] = params[NBODY + i].re / (z_G-params[i]) / (z_G-params[i]);
		fz1 = fz1 + items[i];
	}
	J = 1-(fz1 * conj(fz1)).re;

	fz2=0;
	for(int i=0;i<NBODY;i++)
	{
		items[i] = -f_T(2.) * items[i] / (z_G-params[i]);
		fz2 = fz2 + items[i];
	}

	fz_hat1=0;
	for(int i=0;i<NBODY;i++)
	{
		items[i] = params[NBODY + i].re / (z_hat-params[i]) / (z_hat-params[i]);
		fz_hat1 = fz_hat1 + items[i];
	}
	J_hat = f_T(1.) - fz1*fz_hat1;


	temp = J_hat * conj(fz2) * fz1;
	temp = (temp - conj(temp)*fz_hat1) / (J*J_hat*J_hat);

	return abs_c(temp);
}
template _Device float GhostTst(const complex_t<float>* params, const complex_t<float>& z_G, const complex_t<float>& z_hat);
template _Device double GhostTst(const complex_t<double>* params, const complex_t<double>& z_G, const complex_t<double>& z_hat);

template<isFloating f_T>
_Device f_T WedgeDxDDx(const complex_t<f_T>& z, const complex_t<f_T>& zeta, const complex_t<f_T>* params, const f_T& rho, const complex_t<f_T>& loc_centre, complex_t<f_T>& dz_output)
{
	complex_t<f_T> items[NBODY];
	complex_t<f_T> pzeta_pconjz, ppzeta_ppconjz;			// d means derivitive to theta
	f_T J;
	complex_t<f_T>& dzeta = items[0];
	complex_t<f_T>& dz = items[1];
	f_T& res = pzeta_pconjz.re;

	pzeta_pconjz=complex_t<f_T>(0.,0.);
	for(int i=0;i<NBODY;i++)
	{
		items[i] = params[NBODY + i].re / (conj(z)-params[i]) / (conj(z)-params[i]);
		pzeta_pconjz = pzeta_pconjz + items[i];
		// pzeta_pconjz += items[i];
	}
	// J = 1-(pzeta_pconjz * conj(pzeta_pconjz)).re;
	J = 1 - norm(pzeta_pconjz);

	ppzeta_ppconjz=complex_t<f_T>(0.,0.);
	for(int i=0;i<NBODY;i++)
	{
		items[i] = -f_T(2.) * items[i] / (conj(z)-params[i]);
		ppzeta_ppconjz = ppzeta_ppconjz + items[i];
		// ppzeta_ppconjz += items[i];
	}

	dzeta = complex_t<f_T>(0.,1.) * (zeta-loc_centre);		// zeta - loc_centre = rho*exp(i*theta)


	dz = (dzeta - pzeta_pconjz*conj(dzeta)) / J;		// equation20 in VBB

	dz_output = dz;

	res = (rho*rho + imag(dz*dz * dzeta * conj(ppzeta_ppconjz))) / J;

	return res;
}
template _Device float WedgeDxDDx(const complex_t<float>& z, const complex_t<float>& zeta, const complex_t<float>* params, const float& rho, const complex_t<float>& loc_centre, complex_t<float>& dz_output);
template _Device double WedgeDxDDx(const complex_t<double>& z, const complex_t<double>& zeta, const complex_t<double>* params, const double& rho, const complex_t<double>& loc_centre, complex_t<double>& dz_output);


_Device float sqrt_f_T(const float& f)
{
		return __sqrt2(f);
}

_Device double sqrt_f_T(const double& f)
{
	return __sqrt1(f);
}

template<isFloating f_T>
_Device void ConnectOrder(complex_t<f_T>* posi_here, complex_t<f_T>* posi_other, int len_here, int len_other, int* order, bool print)
{
	f_T Distance[6];
	int jbest=0;
	f_T temp;
	complex_t<f_T> temp_c;
	if((len_here==2))
	{
		if(len_other==2)
		{
			// f_T Distance[2];	// 01,10
			// Distance[0] = abs_c(posi_here[0]-posi_other[0]) + abs_c(posi_here[1]-posi_other[1]);
			// Distance[1] = abs_c(posi_here[0]-posi_other[1]) + abs_c(posi_here[1]-posi_other[0]);

			SUB(temp_c,posi_here[0],posi_other[0]);
			temp = abs_c(temp_c);
			Distance[0] = temp;
			SUB(temp_c,posi_here[0],posi_other[1]);
			temp = abs_c(temp_c);
			Distance[1] = temp;
			SUB(temp_c,posi_here[1],posi_other[1]);
			temp = abs_c(temp_c);
			Distance[0] += temp;
			SUB(temp_c,posi_here[1],posi_other[0]);
			temp = abs_c(temp_c);
			Distance[1] += temp;

			if(Distance[0] < Distance[1])
			{
				order[0] = 0;
				order[1] = 1;
			}
			else
			{
				order[0] = 1;
				order[1] = 0;
			}			
		}
		else	// len_other==3
		{
			// f_T Distance[6];	// 01,02,12,10,20,21
			// Distance[0] = abs_c(posi_here[0]-posi_other[0]) + abs_c(posi_here[1]-posi_other[1]);
			// Distance[1] = abs_c(posi_here[0]-posi_other[0]) + abs_c(posi_here[1]-posi_other[2]);
			// Distance[2] = abs_c(posi_here[0]-posi_other[1]) + abs_c(posi_here[1]-posi_other[2]);
			// Distance[3] = abs_c(posi_here[0]-posi_other[1]) + abs_c(posi_here[1]-posi_other[0]);
			// Distance[4] = abs_c(posi_here[0]-posi_other[2]) + abs_c(posi_here[1]-posi_other[0]);
			// Distance[5] = abs_c(posi_here[0]-posi_other[2]) + abs_c(posi_here[1]-posi_other[1]);

			SUB(temp_c,posi_here[0],posi_other[0]);
			temp = abs_c(temp_c);
			Distance[0] = temp;
			Distance[1] = temp;
			SUB(temp_c,posi_here[0],posi_other[1]);
			temp = abs_c(temp_c);
			Distance[2] = temp;
			Distance[3] = temp;	
			SUB(temp_c,posi_here[0],posi_other[2]);
			temp = abs_c(temp_c);
			Distance[4] = temp;
			Distance[5] = temp;

			SUB(temp_c,posi_here[1],posi_other[1]);
			temp = abs_c(temp_c);
			Distance[0] += temp;
			Distance[5] += temp;		
			SUB(temp_c,posi_here[1],posi_other[2]);
			temp = abs_c(temp_c);
			Distance[1] += temp;
			Distance[2] += temp;
			SUB(temp_c,posi_here[1],posi_other[0]);
			temp = abs_c(temp_c);
			Distance[3] += temp;
			Distance[4] += temp;
			// int jbest=0;
			// f_T temp = Distance[0];
			temp = Distance[0];
			for(int j=1;j<6;j++)
			{
				if(Distance[j] < temp)
				{
					temp = Distance[j];
					jbest = j;
				}
			}
			order[0] = jbest/2;		// 0/2=0,1/2=0, 2/2=1.3/2=1, 4/2=2,5/2=2
			order[1] = ((jbest+3)/2) % 3;	// 1,2,2,0,0,1
		}
	}
	else		// len_here==3
	{
		if(len_other==2)
		{
			// f_T Distance[6];	// 01X，0X1, 1X0, 10X, X01, X10
			// Distance[0] = abs_c(posi_here[0]-posi_other[0]) + abs_c(posi_here[1]-posi_other[1]);
			// Distance[1] = abs_c(posi_here[0]-posi_other[0]) + abs_c(posi_here[2]-posi_other[1]);
			// Distance[2] = abs_c(posi_here[0]-posi_other[1]) + abs_c(posi_here[2]-posi_other[0]);
			// Distance[3] = abs_c(posi_here[0]-posi_other[1]) + abs_c(posi_here[1]-posi_other[0]);
			// Distance[4] = abs_c(posi_here[1]-posi_other[0]) + abs_c(posi_here[2]-posi_other[1]);
			// Distance[5] = abs_c(posi_here[1]-posi_other[1]) + abs_c(posi_here[2]-posi_other[0]);

			SUB(temp_c,posi_here[0],posi_other[0]);
			temp = abs_c(temp_c);
			Distance[0] = temp;
			Distance[1] = temp;
			SUB(temp_c,posi_here[0],posi_other[1]);
			temp = abs_c(temp_c);
			Distance[2] = temp;
			Distance[3] = temp;	
			SUB(temp_c,posi_here[1],posi_other[0]);
			temp = abs_c(temp_c);
			Distance[4] = temp;
			Distance[3] += temp;

			SUB(temp_c,posi_here[1],posi_other[1]);
			temp = abs_c(temp_c);
			Distance[0] += temp;
			Distance[5] = temp;					
			SUB(temp_c,posi_here[2],posi_other[1]);
			temp = abs_c(temp_c);
			Distance[1] += temp;
			Distance[4] += temp;
			SUB(temp_c,posi_here[2],posi_other[0]);
			temp = abs_c(temp_c);
			Distance[2] += temp;
			Distance[5] += temp;

			// int jbest=0;
			// f_T temp = Distance[0];
			temp = Distance[0];
			for(int j=1;j<6;j++)
			{
				if(Distance[j] < temp)
				{
					temp = Distance[j];
					jbest = j;
				}
			}
			order[0] = jbest/2;		// 0/2=0,1/2=0, 2/2=1.3/2=1, 4/2=2,5/2=2
			order[1] = ((jbest+3)/2) % 3;	// 1,2,2,0,0,1
			order[2] = 2 - (jbest%3);	// 2,1,0,2,1,0

			for(int j=0;j<6;j++)
			{
				if(order[j]==2){order[j]=-1;}	// 2为虚解，标记为未连接
			}
		}
		else	// len_other==3
		{
			// f_T Distance[6];	// 012, 021, 120, 102, 201, 210
			// Distance[0] = abs_c(posi_here[0]-posi_other[0]) + abs_c(posi_here[1]-posi_other[1]) + abs_c(posi_here[2]-posi_other[2]);
			// Distance[1] = abs_c(posi_here[0]-posi_other[0]) + abs_c(posi_here[1]-posi_other[2]) + abs_c(posi_here[2]-posi_other[1]);
			// Distance[2] = abs_c(posi_here[0]-posi_other[1]) + abs_c(posi_here[1]-posi_other[2]) + abs_c(posi_here[2]-posi_other[0]);
			// Distance[3] = abs_c(posi_here[0]-posi_other[1]) + abs_c(posi_here[1]-posi_other[0]) + abs_c(posi_here[2]-posi_other[2]);
			// Distance[4] = abs_c(posi_here[0]-posi_other[2]) + abs_c(posi_here[1]-posi_other[0]) + abs_c(posi_here[2]-posi_other[1]);
			// Distance[5] = abs_c(posi_here[0]-posi_other[2]) + abs_c(posi_here[1]-posi_other[1]) + abs_c(posi_here[2]-posi_other[0]);

			SUB(temp_c,posi_here[0],posi_other[0]);
			temp = abs_c(temp_c);
			Distance[0] = temp;
			Distance[1] = temp;
			SUB(temp_c,posi_here[0],posi_other[1]);
			temp = abs_c(temp_c);
			Distance[2] = temp;
			Distance[3] = temp;
			SUB(temp_c,posi_here[0],posi_other[2]);
			temp = abs_c(temp_c);
			Distance[4] = temp;
			Distance[5] = temp;

			SUB(temp_c,posi_here[1],posi_other[1]);
			temp = abs_c(temp_c);
			Distance[0] += temp;
			Distance[5] += temp;
			SUB(temp_c,posi_here[1],posi_other[2]);
			temp = abs_c(temp_c);
			Distance[1] += temp;
			Distance[2] += temp;
			SUB(temp_c,posi_here[1],posi_other[0]);
			temp = abs_c(temp_c);
			Distance[3] += temp;
			Distance[4] += temp;

			SUB(temp_c,posi_here[2],posi_other[2]);
			temp = abs_c(temp_c);
			Distance[0] += temp;
			Distance[3] += temp;
			SUB(temp_c,posi_here[2],posi_other[1]);
			temp = abs_c(temp_c);
			Distance[1] += temp;
			Distance[4] += temp;
			SUB(temp_c,posi_here[2],posi_other[0]);
			temp = abs_c(temp_c);
			Distance[2] += temp;
			Distance[5] += temp;

			// int jbest=0;
			// f_T temp = Distance[0];
			temp = Distance[0];
			for(int j=1;j<6;j++)
			{
				if(Distance[j] < temp)
				{
					temp = Distance[j];
					jbest = j;
				}
			}
			order[0] = jbest/2;		// 0/2=0,1/2=0, 2/2=1,3/2=1, 4/2=2,5/2=2
			order[1] = ((jbest+3)/2) % 3;	// 1,2,2,0,0,1
			order[2] = 2 - (jbest%3);	// 2,1,0,2,1,0
		}
	}

}
template _Device void ConnectOrder(complex_t<float>* posi_here, complex_t<float>* posi_other, int len_here, int len_other, int* order, bool print);
template _Device void ConnectOrder(complex_t<double>* posi_here, complex_t<double>* posi_other, int len_here, int len_other, int* order, bool print);


template<isFloating f_T>
_Device int PrevSrcIdx(const src_pt_t<f_T>* margin_pts, const int& idx_here)
{
	bool skip;
	int idx_prev;
	idx_prev = margin_pts[idx_here].prev_src_idx;
	skip = margin_pts[idx_prev].skip;
	while(skip)
	{
		idx_prev = margin_pts[idx_prev].prev_src_idx;
		skip = margin_pts[idx_prev].skip;
	}
	return idx_prev;
}
template _Device int PrevSrcIdx(const src_pt_t<float>* margin_pts, const int& idx_here);
template _Device int PrevSrcIdx(const src_pt_t<double>* margin_pts, const int& idx_here);


template<isFloating f_T>
_Device int NextSrcIdx(const src_pt_t<f_T>* margin_pts, const int& idx_here)
{
	bool skip;
	int idx_next;
	idx_next = margin_pts[idx_here].next_src_idx;
	skip = margin_pts[idx_next].skip;
	while(skip)
	{
		idx_next = margin_pts[idx_next].next_src_idx;
		skip = margin_pts[idx_next].skip;
	}
	return idx_next;
}
template _Device int NextSrcIdx(const src_pt_t<float>* margin_pts, const int& idx_here);
template _Device int NextSrcIdx(const src_pt_t<double>* margin_pts, const int& idx_here);


template<isFloating f_T>
_Host_Device f_T Jacobi_M(const complex_t<f_T>& z, const f_T s, const f_T q)
{
	// light star frame
	f_T jacobi,m1,m2,zr,zi;
	f_T temp1,temp2,temp3;
	m1 = 1./(1.+q);
	m2 = m1*q;
	zr = z.re;
	zi = z.im;

	// 交叉项
	temp1 = (zr+s)*zr + zi*zi;
	temp1 *= temp1;
	temp2 = s*zi;
	temp2 *= temp2;
	temp3 = temp1 + temp2;
	temp3 = (temp1 - temp2)/temp3/temp3;
	temp3 *= 2*m1*m2;

	// 主星项
	temp1 = (zr+s)*(zr+s) + zi*zi;
	temp1 = (m1-temp1)/temp1;
	temp1 = (2+temp1)*temp1;

	// 行星项
	temp2 = m2/norm(z);
	temp2 *= temp2;

	jacobi = - (temp1 + temp2 + temp3);

	return jacobi;
}
template _Host_Device float Jacobi_M(const complex_t<float>& z, const float s, const float q);
template _Host_Device double Jacobi_M(const complex_t<double>& z, const double s, const double q);

template<isFloating f_T>
_Device bool PolishPP2(image_pt_t<f_T>* imgs, f_T* jacobis, const complex_t<f_T>& zeta,int idx)
{
	int LCR[3];
	f_T LCR_re[3];

	int fifth;
	int same_sign = -1;
	f_T& Jmax = LCR_re[0];
	Jmax=0;
	f_T& tempf = LCR_re[1];

	int count=0;

	for(int j=0;j<5;j++)
	{
		if((imgs[j].position.im>0) == (zeta.im>0))
		{
			same_sign = j;			// only one candidate
			continue;
		}
		tempf = fabs(jacobis[j]);
		if(tempf > Jmax)
		{
			Jmax = tempf;
			fifth = j;
		}
	}
	if(same_sign==-1){printf("-1!");return true;}
	for(int j=0;j<5;j++)
	{
		if(j!=same_sign && j!=fifth)
		{
			LCR[count] = j;
			LCR_re[count] = imgs[j].position.re;
			count++;
		}
	}

	BubbleArgSort(LCR,LCR_re,3);

	if(imgs[LCR[0]].parity != 1)
	{
		imgs[LCR[0]].parity = 1;
		jacobis[LCR[0]] = -jacobis[LCR[0]];
		// printf("changeL: idx:%d, j: %d\n",idx,LCR[0]);
	}
	if(imgs[LCR[1]].parity != 0)
	{
		imgs[LCR[1]].parity = 0;
		jacobis[LCR[1]] = -jacobis[LCR[1]];
		// printf("changeC: idx:%d, j: %d\n",idx,LCR[1]);
	}
	if(imgs[LCR[2]].parity != 1)
	{
		imgs[LCR[2]].parity = 1;
		jacobis[LCR[2]] = -jacobis[LCR[2]];
		// printf("changeR: idx:%d, j: %d\n",idx,LCR[2]);
	}


	return false;
}
template _Device bool PolishPP2(image_pt_t<float>* imgs, float* jacobis, const complex_t<float>& zeta,int idx);
template _Device bool PolishPP2(image_pt_t<double>* imgs, double* jacobis, const complex_t<double>& zeta,int idx);