#include"CoeffGenerator.h"

template<isFloating f_T>
_Host_Device void ParamsReal2Complex(const f_T* params_lens, complex_t<f_T>* params)
{
	const f_T& s = params_lens[0];
	const f_T& q = params_lens[1];

	params[2] = complex_t<f_T>(1 / (1 + q), 0);												    	// m0
	params[3] = complex_t<f_T>(q, 0) * params[2];													// m1
	params[0] = complex_t<f_T>(-s, 0) * params[3];													// z0 = - q/(1+q) * s
	params[1] = complex_t<f_T>(s, 0) * params[2];													// z1 = 1/(1+q) * s
}
template _Host_Device void ParamsReal2Complex(const float* params_lens, complex_t<float>* params);
template _Host_Device void ParamsReal2Complex(const double* params_lens, complex_t<double>* params);

// 多项式相乘
// 如要数乘，对一个数取址& 传入，degree = 0
// 仅限数乘时可以把计算结果存回原地
template<typename f_T>
_Host_Device void PolyMultiply(const f_T* coeff1, const int degree1, const f_T* coeff2, const int degree2, f_T* coeff_res)
{
	// 初始化
    // complex_t<f_T> temp;
	for (int i = 0; i < degree1 + degree2 + 1; i++) {
		coeff_res[i] = f_T(0.);
	}
	// 相乘叠加
	for (int i = 0; i < degree1 + 1; i++) {								// 如要检查 i<degree1+1 ，只需记住 i 此时代表 “i次项系数”
		for (int j = 0; j < degree2 + 1; j++) {
            // // MUL(temp,coeff1[i],coeff2[j]);
            // temp = coeff1[i] * coeff2[j];
            // // ADD(coeff_res[i + j],temp,coeff_res[i + j]);
            coeff_res[i + j] = coeff1[i] * coeff2[j] + coeff_res[i + j];
		}
	}
}
template _Host_Device void PolyMultiply(const complex_t<float>* coeff1, const int degree1, const complex_t<float>* coeff2, const int degree2, complex_t<float>* coeff_res);
template _Host_Device void PolyMultiply(const complex_t<double>* coeff1, const int degree1, const complex_t<double>* coeff2, const int degree2, complex_t<double>* coeff_res);
template _Host_Device void PolyMultiply(const float* coeff1, const int degree1, const float* coeff2, const int degree2, float* coeff_res);
template _Host_Device void PolyMultiply(const double* coeff1, const int degree1, const double* coeff2, const int degree2, double* coeff_res);

// 数乘，把数放在前一个输入
template<typename f_T>
_Host_Device void PolyMultiply(const f_T coeff0, const int degree0, const f_T* coeff2, const int degree2, f_T* coeff_res)
{
	for (int i = 0; i < degree2 + 1; i++) {
        // MUL(coeff_res[i], coeff0,coeff2[i]);
        coeff_res[i] = coeff0 * coeff2[i];
	}
}
template _Host_Device void PolyMultiply(const complex_t<float> coeff0, const int degree0, const complex_t<float>* coeff2, const int degree2, complex_t<float>* coeff_res);
template _Host_Device void PolyMultiply(const complex_t<double> coeff0, const int degree0, const complex_t<double>* coeff2, const int degree2, complex_t<double>* coeff_res);
template _Host_Device void PolyMultiply(const float coeff0, const int degree0, const float* coeff2, const int degree2, float* coeff_res);
template _Host_Device void PolyMultiply(const double coeff0, const int degree0, const double* coeff2, const int degree2, double* coeff_res);


// 多项式相加
// 注意：没有初始化
// 可以存在原地
template<typename f_T>
_Host_Device void PolyAdd(const f_T* coeff1,const int degree1, const f_T* coeff2, const int degree2, f_T* coeff_res)
{
	for (int i = 0; i < min(degree1, degree2) + 1; i++) {
        // ADD(coeff_res[i],coeff1[i],coeff2[i]);
        coeff_res[i] = coeff1[i] + coeff2[i];
	}
	for (int i = min(degree1, degree2) + 1; i < max(degree1, degree2) + 1; i++) {
		if (degree1 < degree2) {
			coeff_res[i] = coeff2[i];
		}
		else {
			coeff_res[i] = coeff1[i];
		}
	}
}
template _Host_Device void PolyAdd(const complex_t<float>* coeff1, const int degree1, const complex_t<float>* coeff2, const int degree2, complex_t<float>* coeff_res);
template _Host_Device void PolyAdd(const complex_t<double>* coeff1, const int degree1, const complex_t<double>* coeff2, const int degree2, complex_t<double>* coeff_res);
template _Host_Device void PolyAdd(const float* coeff1, const int degree1, const float* coeff2, const int degree2, float* coeff_res);
template _Host_Device void PolyAdd(const double* coeff1, const int degree1, const double* coeff2, const int degree2, double* coeff_res);

// GPU版和CPU版稍有不同，将loc_src从params中单独拎了出来
template<isFloating f_T>
_Host_Device void PolyCoeff(const complex_t<f_T>* params,const complex_t<f_T>& loc_src, complex_t<f_T>* coeff_res)
{
	// 首先分子分母同乘以 (z-z0)(z-z1)，计算此时分母
	// 定义与初始化
	const int N_body = NBODY;      // 目前只考虑三体
	const int length = N_body * N_body + 2;				// 次数 N^2+1 多一个零次项
	for (int i = 0; i < length; i++) {
		coeff_res[i] = complex_t<f_T>(0, 0);
	}
	complex_t<f_T> denominator[N_body * (N_body + 1)];       // 通分后分母多项式的系数，使用 denominator + i*(N_body+1) 来访问第i个数组 （因为数组名已经是一个指针了，不需要再取址）
	for (int i = 0; i < N_body * (N_body + 1); i++) {
		denominator[i] = complex_t<f_T>(0, 0);
	}
	complex_t<f_T> root01[N_body + 1];						// 临时储存的根式乘积，最终结果(z-z0)(z-z1)
	for (int i = 0; i < N_body + 1; i++) {
		root01[i] = complex_t<f_T>(0, 0);
	}
	complex_t<f_T> temp_res_1[length];						// 根式，存在一起，使用 root + i*2  来访问第i个数组 (z-z_i)（因为数组名已经是一个指针了，不需要再取址）
	for (int i = 0; i < 2 * N_body; i += 2) {			// 因为后面还要用到空间，且这时候root已经功成身退，故它还兼任临时空间的作用
		temp_res_1[i] = -params[i / 2];             	// 零次项 - z_i
		temp_res_1[i + 1] = complex_t<f_T>(1, 0);				// 一次项 1
	}
	complex_t<f_T> temp_poly[N_body + 1];					    // 临时存放处，总是与(z-z0)(z-z1)(z-z2)有关
	complex_t<f_T> temp_res_0[length];						// 临时存放处，总与分母（denominator）的相乘有关
    complex_t<f_T> temp_conj;                                 // 临时存放处，总与共轭计算有关

	// 先在denominator的第一部分储存公用的系数 zeta_bar*(z-z0)(z-z1) + m0(z-z1) + m1(z-z0)
	// 没写循环，仅限二体
	// 先算 m0(z-z1) + m1(z-z0) ，1次 赋给denominator[0:2]
	PolyMultiply(params[2], 0, temp_res_1 + 2, 1, root01);		// root01[0:2] = m0*(z-z1)

	// // for (int i = 0; i < 2; i++) {
	// // 	coeff_res[i] = root01[i];
	// // }
	// // coeff_res[2] = params[2];
	// // coeff_res[3] = *(temp_res_1 + 2);
	// // coeff_res[4] = *(temp_res_1 + 3);

	PolyAdd(denominator, 1, root01, 1, denominator);			// deno[0:2] += m0*(z-z1)

	PolyMultiply(params[3], 0, temp_res_1 + 0, 1, root01);		// root01[0:2] = m1*(z-z0)
	PolyAdd(denominator, 1, root01, 1, denominator);			// deno[0:2] += m1(z-z0)


	// 再算 zeta_bar*(z-z0)(z-z1), 
	PolyMultiply(temp_res_1 + 0, 1, temp_res_1 + 2, 1, root01);	// root012[0:2] = (z-z0)(z-z1)
    PolyMultiply(conj(loc_src), 0, root01, 2, temp_res_1);// temp_res_1 = (z-z0)(z-z1) * zeta_bar      
	PolyAdd(temp_res_1, 2, denominator, 2, denominator);		// deno[0:3] += (z-z0)(z-z1) * zeta_bar



	// 拷贝分母的公共部分
	for (int i = 1; i < N_body; i++) {
		for (int j = 0; j < N_body + 1; j++) {
			denominator[i * (N_body + 1) + j] = denominator[j];
		}
	}

	// for (int i = 0; i < 3; i++) {
	// 	coeff_res[i] = denominator[i];
	// }

	// 至此，denominator 中存了两个公共的多项式 zeta_bar*(z-z0)(z-z1) + m0(z-z1) + m1(z-z0)
	// root01 中存了(z-z0)
	// 真正的分母是公共部分 +  root01 * (-z_i_bar)
	for (int i = 0; i < N_body; i++) {
		PolyMultiply(-conj(params[i]), 0, root01, 2, temp_poly);		// temp_res_1[0:2] = -z_i_bar * (z-z0)(z-z1)	
		PolyAdd(temp_poly, 2, denominator + (i * (N_body + 1)), 2, denominator + (i * (N_body + 1)));	// denominator + (i * (N_body + 1)) 是第i个分母的开头地址
	}
	// 至此，两个分母也得到了
	// 
	// 
	// 
	// 开始通分两个分母
	// 先把两个分式加到结果里去
	// 先不写循环了，先写两体特例
	PolyMultiply(params[2], 0, root01, 2, temp_poly);							// 第 0 个分子，m_0*(z-z0)(z-z1)，存进temp_poly
	PolyMultiply(temp_poly, 2, denominator + 1 * (N_body + 1), 2, temp_res_1);	// nume0 * deno1, 分子0*分母1, 结果放进 temp_res_1（4次）
	PolyAdd(temp_res_1, 4, coeff_res, 4, coeff_res);							// coeff_res = nume0 * deno1

	PolyMultiply(params[3], 0, root01, 2, temp_poly);							// 第 1 个分子，m_1*(z-z0)(z-z1)，存进temp_poly
	PolyMultiply(temp_poly, 2, denominator + 0 * (N_body + 1), 2, temp_res_1);	// nume1 * deno0, 分子1*分母0, 结果放进 temp_res_1（4次）
	PolyAdd(temp_res_1, 4, coeff_res, 4, coeff_res);							// coeff_res = nume0 * deno1 + nume1 * deno0

	// 还有最后一项 -(z - zeta)*dd01
	// 此时其实root01不用了，可优化掉temp_res_0
	temp_poly[1] = complex_t<f_T>(-1, 0);
	temp_poly[0] = loc_src;											            // temp_poly = (-z+zeta)
	PolyMultiply(temp_poly, 1, denominator + 0 * (N_body + 1), 2, temp_res_0);	// temp_res_0 = (-z+zeta)*deno0,3次
	PolyMultiply(temp_res_0, 3, denominator + 1 * (N_body + 1), 2, temp_res_1);	// temp_res_1 = (-z+zeta)*deno0*deno1,5次
	PolyAdd(temp_res_1, 5, coeff_res, 5, coeff_res);							// coeff_res = nd01 + nd10 -(z-zeta)dd01


	// // 由于分母（denominator）中最高次项系数不定，归一化：
	// for (int i = length - 1; i > 0; i--) {
	// 	if (coeff_res[i] != complex_t<f_T>(0, 0)) {
	// 		temp_poly[0] = coeff_res[i];								// 反正是temp，空间借我用用
	// 		break;
	// 	}
	// }
	// for (int i = 0; i < length; i++) {
	// 	coeff_res[i] = coeff_res[i] / temp_poly[0];
	// }
}
template _Host_Device void PolyCoeff(const complex_t<float>* params,const complex_t<float>& loc_src, complex_t<float>* coeff_res);
template _Host_Device void PolyCoeff(const complex_t<double>* params,const complex_t<double>& loc_src, complex_t<double>* coeff_res);

// f(z_bar)，注意输入的是 z 不是 zbar
template<isFloating f_T>
_Host_Device complex_t<f_T> f_zbar(const complex_t<f_T>& z, const complex_t<f_T>* params)
{
	complex_t<f_T> res(0,0);
	for (int i = 0; i < NBODY; i++) {
		res = res - params[NBODY + i] / (conj(z) - conj(params[i]));
	}
	return res;
}
template _Host_Device complex_t<float> f_zbar(const complex_t<float>& z, const complex_t<float>* params);
template _Host_Device complex_t<double> f_zbar(const complex_t<double>& z, const complex_t<double>* params);

// mass centre frame
template<isFloating f_T>
_Host_Device complex_t<f_T> f_zbar(const complex_t<f_T>& z, const f_T s, const f_T m1, const f_T m2)
{
	complex_t<f_T> res(0,0);
	// for (int i = 0; i < NBODY; i++) {
	// 	res = res - params[NBODY + i] / (conj(z) - conj(params[i]));
	// }
	// return res;
	res -= m1 / (conj(z) + s*m2);
	res -= m2 / (conj(z) - s*m1);
	return res;
}
template _Host_Device complex_t<float> f_zbar(const complex_t<float>& z, const float s, const float m1, const float m2);
template _Host_Device complex_t<double> f_zbar(const complex_t<double>& z, const double s, const double m1, const double m2);


// d(f(z_bar)) / d(z_bar)，注意对谁求导
template<isFloating f_T>
_Host_Device complex_t<f_T> D_f_zbar(const complex_t<f_T>& z, const complex_t<f_T>* params)
{
	complex_t<f_T> res(0, 0);
	for (int i = 0; i < NBODY; i++) {
		res = res + params[NBODY + i] / ((conj(z) - conj(params[i])) * (conj(z) - conj(params[i])));
	}
	return res;
}
template _Host_Device complex_t<float> D_f_zbar(const complex_t<float>& z, const complex_t<float>* params);
template _Host_Device complex_t<double> D_f_zbar(const complex_t<double>& z, const complex_t<double>* params);

// Jacobi 放大率的倒数
template <isFloating f_T>
_Host_Device f_T Jacobi(const complex_t<f_T>& z, const complex_t<f_T>* params)
{
	complex_t<f_T> d_f = D_f_zbar(z, params);
	f_T jacobi = (1 - real(d_f * conj(d_f)));
	return jacobi;
}
template _Host_Device float Jacobi(const complex_t<float>& z, const complex_t<float>* params);
template _Host_Device double Jacobi(const complex_t<double>& z, const complex_t<double>* params);

template <typename f_T>
_Host_Device f_T PolyValue(const f_T* coeff, const f_T& z, int degree)
{
	f_T res = coeff[degree];
	for(int i=0;i<degree;i++)
	{
		res = res * z;
		res = res + coeff[degree - i - 1];
	}
	return res;
}
template _Host_Device complex_t<float> PolyValue(const complex_t<float>* coeff, const complex_t<float>& z, int degree);
template _Host_Device complex_t<double> PolyValue(const complex_t<double>* coeff, const complex_t<double>& z, int degree);
template _Host_Device float PolyValue(const float* coeff, const float& z, int degree);
template _Host_Device double PolyValue(const double* coeff, const double& z, int degree);

template <typename f_T>
_Host_Device void PolyDer(const f_T* coeff,f_T* res,const int prev_degree)
{
	for(int i=0;i<prev_degree;i++)
	{
		res[i] = coeff[i+1]* f_T(float(i+1));	// 先转化为浮点泪，再转化为f_T（因为没有整数转complex_t的构造函数）
	}
}
template _Host_Device void PolyDer(const complex_t<float>* coeff,complex_t<float>* res,const int prev_degree);
template _Host_Device void PolyDer(const complex_t<double>* coeff,complex_t<double>* res,const int prev_degree);
template _Host_Device void PolyDer(const float* coeff,float* res,const int prev_degree);
template _Host_Device void PolyDer(const double* coeff,double* res,const int prev_degree);


// 重的在左边
// 定义几何中心和质量中心 geometry系和mass系
// moving 总是非负的，减去就是坐标更小，加上就是坐标更大

template <isFloating f_T>
_Host_Device f_T PolyCoeffGCOld(const f_T s, const f_T q, const complex_t<f_T>& zeta_M, complex_t<f_T>* coeffs)
{
	f_T Dm = (q-1.)/(1.+q) / 2.;	// (m2-m1)/2, always negative since m2 <= m1
	f_T moving = -Dm * s;		// GC坐标系与MC坐标系的距离，moving总是一个正实数，具体怎么用看实际情况
	complex_t<f_T> zeta_G = zeta_M - moving;		// zeta in MC, zeta_G in GC
	f_T z1_G = -s/2;
	// m = (m1+m2)/2 = 0.5

	// f_T z1_G2 = z1_G*z1_G;
	// f_T norm_zeta = norm(zeta);
	// complex_t<f_T> zeta_G_c = conj(zeta_G);
	
	coeffs[5] = z1_G*z1_G - conj(zeta_G)*conj(zeta_G);
	coeffs[4] = -conj(zeta_G) + norm(zeta_G)*conj(zeta_G) - 2.*Dm*z1_G - zeta_G*z1_G*z1_G;
	coeffs[3] = 2.*norm(zeta_G) + 4*Dm*conj(zeta_G)*z1_G + 2.*z1_G*z1_G*conj(zeta_G)*conj(zeta_G) - 2*z1_G*z1_G*z1_G*z1_G;
	coeffs[2] = zeta_G + 2.*Dm*z1_G - 4.*Dm*norm(zeta_G)*z1_G - 2.*norm(zeta_G)*conj(zeta_G)*z1_G*z1_G + 4.*Dm*z1_G*z1_G*z1_G + 2.*z1_G*z1_G*z1_G*zeta_G*z1_G;
	coeffs[1] = -4*Dm*zeta_G*z1_G - 4.*Dm*Dm*z1_G*z1_G - z1_G*z1_G - 2*norm(zeta_G)*z1_G*z1_G - 4.*Dm*conj(zeta_G)*z1_G*z1_G*z1_G - conj(zeta_G)*conj(zeta_G)*z1_G*z1_G*z1_G*z1_G + z1_G*z1_G*z1_G*z1_G*z1_G*z1_G;
	coeffs[0] = z1_G*z1_G* (4.*Dm*Dm*zeta_G + 2.*Dm*z1_G + 4.*Dm*norm(zeta_G)*z1_G + conj(zeta_G)*z1_G*z1_G + norm(zeta_G)*conj(zeta_G)*z1_G*z1_G - 2.*Dm*z1_G*z1_G*z1_G - zeta_G*z1_G*z1_G*z1_G*z1_G);




	// don't normalize, even wrose
	// for(int i=0;i<LENGTH;i++)
	// {
	// 	coeffs[i] = coeffs[i] / coeffs[LENGTH-1];
	// }

	return moving;

}
// template _Host_Device float PolyCoeffGC(const float s, const float q, const complex_t<float>& zeta_M, complex_t<float>* coeffs);
template _Host_Device double PolyCoeffGCOld(const double s, const double q, const complex_t<double>& zeta_M, complex_t<double>* coeffs);

template <isFloating f_T>
_Host_Device f_T PolyCoeffGC(const f_T s, const f_T q, const complex_t<f_T>& zeta_M, complex_t<f_T>* coeffs)
{
	f_T Dm = (q-1.)/(1.+q);	// (m2-m1), always negative since m2 <= m1
	f_T moving = -Dm * s / 2.;		// GC坐标系与MC坐标系的距离，moving总是一个正实数，具体怎么用看实际情况
	complex_t<f_T> zeta_G = zeta_M - moving;		// zeta in MC, zeta_G in GC
	f_T z1_G = -s/2;
	f_T z1_M = -q/(1.+q)*s;


	complex_t<f_T> lenini;
	f_T conini;
	const complex_t<f_T>& pocchama1 = zeta_M;
	complex_t<f_T> pocchama2;

	// complex_t<f_T> temp_c;
	// f_T& temp0 = temp_c.re;
	// f_T& temp1 = temp_c.im;
	f_T temp0, temp1;

	// temp0 = z1_M + zeta_M.re - 2*moving;
	temp0 = zeta_M.re - s/(1.+q);
	temp1 = z1_M - zeta_M.re;

	// lenini = complex_t<f_T>(temp0, -zeta_M.im) * complex_t<f_T>(temp1, zeta_M.im);
	// lenini = complex_t<f_T>(temp0*temp1 + zeta_M.im**2, )
	// conini = temp0*temp1 - zeta_M.im**2;
	lenini.im = (temp0-temp1)*zeta_M.im;
	temp0 = temp0*temp1;
	temp1 = zeta_M.im * zeta_M.im;
	conini = temp0 - temp1;
	lenini.re = temp0 + temp1;
	// printf("lenni: (%.24f,%.24f),conini:%.24f\n",lenini.re,lenini.im,conini);

	temp0 = 1./(1.+q);
	// pocchama1.re = zeta_M.re + Dm*z1_M - temp0;
	// pocchama1.re = zeta_M.re;
	// pocchama1.im = zeta_M.im;
	// pocchama2.re = z1_M + Dm*zeta_M.re - temp0;
	// pocchama2.re = Dm*zeta_M.re - (2*q)/(1+q)/(1+q)*s;
	pocchama2.re = temp0 * ((q-1.)*zeta_M.re + 2.*z1_M);
	pocchama2.im = zeta_M.im * Dm;
	// printf("p1: (%.24f,%.24f),p2:(%.24f,%.24f)\n",pocchama1.re,pocchama1.im,pocchama2.re,pocchama2.im);
	// printf("temp0: %.24f\n",temp0);
	// printf("pocha test: %.16f,%.16f\n",Dm*z1_M,temp0);
	// printf("pocha test: %.16f,%.16f\n",-(2*q)/(1+q)/(1+q)*s,z1_M - temp0);

	temp0 = z1_G*z1_G;

	coeffs[5] = lenini;
	coeffs[4] = -conj(pocchama1) - zeta_G*lenini;
	coeffs[3] = 2.*conj(zeta_G)*pocchama1 - 2.*temp0*lenini;
	coeffs[2] = pocchama1 + 2.*Dm*z1_G*conini + 2.*zeta_G*temp0*lenini;
	coeffs[1] = -Dm*z1_G*pocchama1 - z1_G*pocchama2 - 2.*conj(zeta_G)*temp0*pocchama1 + temp0*temp0*lenini;
	coeffs[0] = temp0*(Dm*pocchama2 + conj(zeta_G)*z1_G*pocchama2 - Dm*z1_G*conini - zeta_G*temp0*lenini);



	// don't normalize, even wrose
	// for(int i=0;i<LENGTH;i++)
	// {
	// 	coeffs[i] = coeffs[i] / coeffs[LENGTH-1];
	// }

	return moving;

}
// template _Host_Device float PolyCoeffGC(const float s, const float q, const complex_t<float>& zeta_M, complex_t<float>* coeffs);
template _Host_Device double PolyCoeffGC(const double s, const double q, const complex_t<double>& zeta_M, complex_t<double>* coeffs);


_Host_Device float PolyCoeffGC(const float s, const float q, const complex_t<float>& zeta_M, complex_t<float>* coeffs)
{
	return 0.f;
}



template <isFloating f_T>
_Host_Device f_T PolyCoeffVB(const f_T s, const f_T q, const complex_t<f_T>& zeta_M, complex_t<f_T>* coeffs_out)
{
	// s>0, 0<q<=1
	f_T m1,m2,a;
	complex_t<f_T> temp[2],y,yc;
	f_T moving;			// moving 总是正数

	m1 = 1.0 / (1.0 + q);
	m2 = q*m1;
	a = -s;
	moving = s*m1;
	y = zeta_M - moving;
	yc = conj(y);
	


	temp[0] = a*a;
	temp[1] = temp[0] * a;
	// temp[2] = m2*m2;
	// temp[3] = temp[0] * m2*m2;
	// temp[4] = a*m2;
	// temp[5] = a*m1;
	// a = a;			// temp[20] in VBBL
	// m1 = m1;			// temp[21] in VBBL
	// m2 = m2;			// temp[22] in VBBL	

	coeffs_out[0] = temp[0] * m2*m2 * y;
	coeffs_out[1] = a*m2 * (a * (m1 + y * (f_T(2) * yc - a)) - f_T(2) * y);
	coeffs_out[2] = y * (f_T(1) - temp[1] * yc) - a * (m1 + f_T(2) * y * yc * (f_T(1) + m2)) + temp[0] * (yc * (m1 - m2) + y * (f_T(1) + m2 + yc * yc));
	coeffs_out[3] = f_T(2) * y * yc + temp[1] * yc + temp[0] * (yc * (f_T(2) * y - yc) - m1) - a * (y + f_T(2) * yc * (yc * y - m2));
	coeffs_out[4] = yc * (f_T(2) * a + y);
	coeffs_out[4] = yc * (coeffs_out[4] - f_T(1)) - a * (coeffs_out[4] - m1);
	coeffs_out[5] = yc * (a - yc);

	// printf("VB\n");
	// // printf("(%.16f+%.16fj),(%.16f+%.16fj),(%.16f+%.16fj),(%.16f+%.16fj),(%.16f+%.16fj),(%.16f+%.16fj)\n",(temp[0]).re,(temp[0]).im,(temp[1]).re,(temp[1]).im,(temp[2]).re,(temp[2]).im,(temp[3]).re,(temp[3]).im,(temp[4]).re,(temp[4]).im,(temp[5]).re,(temp[5]).im);
	// for(int i=0;i<6;i++)
	// {
	// 	printf("%d: (%.30f+%.30fj)\n",i,(coeffs_out[i]).re,(coeffs_out[i]).im);
	// }
	// printf("y: (%.30f+%.30fj)\n",(y).re,(y).im);



	return moving;

}
template _Host_Device float PolyCoeffVB(const float s, const float q, const complex_t<float>& zeta_M, complex_t<float>* coeffs_out);
template _Host_Device double PolyCoeffVB(const double s, const double q, const complex_t<double>& zeta_M, complex_t<double>* coeffs_out);


template <isFloating f_T>
_Host_Device f_T PolyCoeffVB(const f_T s, const f_T m1, const f_T m2, const complex_t<f_T>& zeta_M, complex_t<f_T>* coeffs_out)
{
	// s>0, 0<q<=1
	f_T a;
	complex_t<f_T> temp[2],y,yc;
	f_T moving;			// moving 总是正数

	a = -s;
	moving = s*m1;
	y = zeta_M - moving;
	yc = conj(y);
	


	temp[0] = a*a;
	temp[1] = temp[0] * a;
	// temp[2] = m2*m2;
	// temp[3] = temp[0] * m2*m2;
	// temp[4] = a*m2;
	// temp[5] = a*m1;
	// a = a;			// temp[20] in VBBL
	// m1 = m1;			// temp[21] in VBBL
	// m2 = m2;			// temp[22] in VBBL	

	coeffs_out[0] = temp[0] * m2*m2 * y;
	coeffs_out[1] = a*m2 * (a * (m1 + y * (f_T(2) * yc - a)) - f_T(2) * y);
	coeffs_out[2] = y * (f_T(1) - temp[1] * yc) - a * (m1 + f_T(2) * y * yc * (f_T(1) + m2)) + temp[0] * (yc * (m1 - m2) + y * (f_T(1) + m2 + yc * yc));
	coeffs_out[3] = f_T(2) * y * yc + temp[1] * yc + temp[0] * (yc * (f_T(2) * y - yc) - m1) - a * (y + f_T(2) * yc * (yc * y - m2));
	coeffs_out[4] = yc * (f_T(2) * a + y);
	coeffs_out[4] = yc * (coeffs_out[4] - f_T(1)) - a * (coeffs_out[4] - m1);
	coeffs_out[5] = yc * (a - yc);

	return moving;

}
template _Host_Device float PolyCoeffVB(const float s, const float m1, const float m2, const complex_t<float>& zeta_M, complex_t<float>* coeffs_out);
template _Host_Device double PolyCoeffVB(const double s, const double m1, const double m2, const complex_t<double>& zeta_M, complex_t<double>* coeffs_out);



template <isFloating f_T>
_Host_Device f_T PolyCoeffVB1(const f_T s, const f_T q, const complex_t<f_T>& zeta_M, complex_t<f_T>* coeffs_out)
{
	// s>0, 0<q<=1
	f_T m1,m2,a;
	complex_t<f_T> temp[3],y,yc,zeta_c;
	f_T moving;			// moving 总是正数

	m1 = 1.0 / (1.0 + q);
	m2 = q*m1;
	// a = -s;
	moving = s*m1;
	y = zeta_M - moving;
	yc = conj(y);

	zeta_c = conj(zeta_M);
	


	temp[0] = m2*s;
	temp[1] = zeta_c + temp[0];				// zeta_c + m2*s = zeta_c - z_1		source to m1    ( yc - a ) or ( yc + s )
	temp[2] = (1-s*s) + s*(m2*s + zeta_c);

	coeffs_out[5] = - yc * temp[1];
	coeffs_out[4] = (norm(y) -f_T(1.) - f_T(2.)*s*yc ) * temp[1] + temp[0];
	coeffs_out[3] = (f_T(2.)*y-s) * (  temp[2]*temp[1]  + complex_t<f_T>(f_T(0.),2.*imag(y)) - m2*s ) + complex_t<f_T>(f_T(0.),4.*imag(y)) * (m2*s - y);
	coeffs_out[2] = zeta_M + s * f_T(2.)*real(y*temp[1]) + temp[0] * ( f_T(2.)*norm(y) + y*s - f_T(2.)*s*yc ) + s*s*norm(y)*temp[1];
	coeffs_out[1] = temp[0] * ( (s + f_T(2.)*y) * temp[2] + s * (complex_t<f_T>(f_T(0.),f_T(2.)*imag(y)*s) - m2)  );
	coeffs_out[0] = temp[0] * temp[0] * y;




	// coeffs_out[0] = temp[0] * m2*m2 * y;
	// coeffs_out[1] = a*m2 * (a * (m1 + y * (f_T(2) * yc - a)) - f_T(2) * y);
	// coeffs_out[2] = y * (f_T(1) - temp[1] * yc) - a * (m1 + f_T(2) * y * yc * (f_T(1) + m2)) + temp[0] * (yc * (m1 - m2) + y * (f_T(1) + m2 + yc * yc));
	// coeffs_out[3] = f_T(2) * y * yc + temp[1] * yc + temp[0] * (yc * (f_T(2) * y - yc) - m1) - a * (y + f_T(2) * yc * (yc * y - m2));
	// coeffs_out[4] = yc * (f_T(2) * a + y);
	// coeffs_out[4] = yc * (coeffs_out[4] - f_T(1)) - a * (coeffs_out[4] - m1);
	// coeffs_out[5] = yc * (a - yc);

	// printf("VB1\n");
	// for(int i=0;i<6;i++)
	// {
	// 	printf("%d: (%.30f+%.30fj)\n",i,(coeffs_out[i]).re,(coeffs_out[i]).im);
	// }
	// printf("y: (%.30f+%.30fj)\n",(y).re,(y).im);



	return moving;

}
template _Host_Device float PolyCoeffVB1(const float s, const float q, const complex_t<float>& zeta_M, complex_t<float>* coeffs_out);
template _Host_Device double PolyCoeffVB1(const double s, const double q, const complex_t<double>& zeta_M, complex_t<double>* coeffs_out);


template <isFloating f_T>
_Host_Device f_T PolyCoeffVB1(const f_T s, const f_T m1, const f_T m2, const complex_t<f_T>& zeta_M, complex_t<f_T>* coeffs_out)
{
	// s>0, 0<q<=1
	// f_T m1,m2,a;
	f_T a;
	complex_t<f_T> temp[3],y,yc,zeta_c;
	f_T moving;			// moving 总是正数

	// m1 = 1.0 / (1.0 + q);
	// m2 = q*m1;
	// a = -s;
	moving = s*m1;
	y = zeta_M - moving;
	yc = conj(y);

	zeta_c = conj(zeta_M);
	


	temp[0] = m2*s;
	temp[1] = zeta_c + temp[0];				// zeta_c + m2*s = zeta_c - z_1		source to m1    ( yc - a ) or ( yc + s )
	temp[2] = (1-s*s) + s*(temp[0] + zeta_c);

	coeffs_out[5] = - yc * temp[1];
	coeffs_out[4] = (norm(y) -f_T(1.) - f_T(2.)*s*yc ) * temp[1] + temp[0];
	coeffs_out[3] = (f_T(2.)*y-s) * (  temp[2]*temp[1]  + complex_t<f_T>(f_T(0.),2.*imag(y)) - m2*s ) + complex_t<f_T>(f_T(0.),4.*imag(y)) * (m2*s - y);
	coeffs_out[2] = zeta_M + s * f_T(2.)*real(y*temp[1]) + temp[0] * ( f_T(2.)*norm(y) + y*s - f_T(2.)*s*yc ) + s*s*norm(y)*temp[1];
	coeffs_out[1] = temp[0] * ( (s + f_T(2.)*y) * temp[2] + s * (complex_t<f_T>(f_T(0.),f_T(2.)*imag(y)*s) - m2)  );
	coeffs_out[0] = temp[0] * temp[0] * y;

	return moving;

}
template _Host_Device float PolyCoeffVB1(const float s, const float m1, const float m2, const complex_t<float>& zeta_M, complex_t<float>* coeffs_out);
template _Host_Device double PolyCoeffVB1(const double s, const double m1, const double m2, const complex_t<double>& zeta_M, complex_t<double>* coeffs_out);


template <isFloating f_T>
_Host_Device f_T PolyCoeffTwinkle(const f_T s, const f_T m1, const f_T m2, const complex_t<f_T>& zeta_M, complex_t<f_T>* coeffs_out)
{
	// s>0, 0<q<=1
	// f_T m1,m2,a;
	f_T a;
	complex_t<f_T> temp[3],y,yc,zeta_c;
	f_T moving;			// moving 总是正数

	// m1 = 1.0 / (1.0 + q);
	// m2 = q*m1;
	// a = -s;
	moving = s*m1;
	y = zeta_M - moving;
	yc = conj(y);

	zeta_c = conj(zeta_M);
	


	temp[0] = m2*s;
	temp[1] = zeta_c + temp[0];				// zeta_c + m2*s = zeta_c - z_1		source to m1    ( yc - a ) or ( yc + s )
	temp[2] = (1-s*s) + s*temp[1];

	coeffs_out[5] = - yc * temp[1];
	coeffs_out[4] = (norm(y) -f_T(1) - f_T(2)*s*yc ) * temp[1] + temp[0];
	coeffs_out[3] = (f_T(2)*y-s) * (  temp[2]*temp[1]  + complex_t<f_T>(f_T(0),2.*y.im) - temp[0] ) + complex_t<f_T>(f_T(0),4.*y.im) * (temp[0] - y);
	// coeffs_out[2] = zeta_M + s * f_T(2.)*real(y*temp[1]) + temp[0] * ( f_T(2.)*norm(y) + y*s - f_T(2.)*s*yc ) + s*s*norm(y)*temp[1];
	// coeff[2] = np.real(temp[2])**2 * np.real(temp[1]) + yi*1j*np.real(temp[2])*(1-s*yr) + s*yi*yi*(2 + s*(np.real(temp[1]) - yi*1j)) + temp[0] * ( 2*norm(y) + s * (-yr + 3*yi*1j) - 1)
	// coeffs_out[2] = temp[2].re * temp[2].re * temp[1].re + complex_t<f_T>(0.,temp[2].re*(f_T(1.)-s*y.re)*y.im) + s*y.im*y.im* (f_T(2.) + s * complex_t<f_T>(temp[1].re, -y.im)) + temp[0] * ( 2*norm(y) + s * complex_t<f_T>(-y.re, 3*y.im) - f_T(1.));
	coeffs_out[2] = temp[2].re * (temp[2].re * temp[1].re + complex_t<f_T>(0.,(f_T(1)-s*y.re)*y.im)) + s*y.im*y.im* (f_T(2) + s * temp[1]) + temp[0] * ( 2*norm(y) + s * complex_t<f_T>(-y.re, 3*y.im) - f_T(1));
	coeffs_out[1] = temp[0] * ( (s + f_T(2)*y) * temp[2] + s * (complex_t<f_T>(f_T(0),f_T(2)*y.im*s) - m2)  );
	coeffs_out[0] = temp[0] * temp[0] * y;

	return moving;

}
template _Host_Device float PolyCoeffTwinkle(const float s, const float m1, const float m2, const complex_t<float>& zeta_M, complex_t<float>* coeffs_out);
template _Host_Device double PolyCoeffTwinkle(const double s, const double m1, const double m2, const complex_t<double>& zeta_M, complex_t<double>* coeffs_out);
