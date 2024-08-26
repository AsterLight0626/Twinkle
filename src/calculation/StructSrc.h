#pragma once

#include"StructComplex.h"
#include"MacroVaribles.h"



// 结构体：像
template < typename f_T >
struct image_pt_t {
    complex_t<f_T> position;
    bool physical;
    bool parity;
	_Host_Device image_pt_t();
	_Host_Device image_pt_t(const image_pt_t<f_T>& other);
};

// 结构体：点源
template < typename f_T>
struct src_pt_t {
    // f_T x, y;
    complex_t<f_T> position;
    image_pt_t<f_T> images[NBODY*NBODY+1];
	f_T wedge[NBODY*NBODY+1];
	complex_t<f_T> dz[NBODY*NBODY+1];

	f_T deltaS[NBODY*NBODY+1];
	f_T deltaS_Err[NBODY*NBODY+1];
	f_T deltaS_new[NBODY*NBODY+1];
	f_T Err_new[NBODY*NBODY+1];	

	f_T special_Err;					// for additional weight including physical test
	// if constexpr (d_T)
	// {
	#if DetailedRecord
		f_T deltaS_p[NBODY*NBODY+1];		
		f_T E1[NBODY*NBODY+1];				// for output only
		f_T E2[NBODY*NBODY+1];				// for output only
		f_T E3[NBODY*NBODY+1];				// for output only
		f_T phys_err[NBODY*NBODY+1];		// for output only
		f_T jacobi[NBODY*NBODY+1];			// for output only
	#endif
	// }
	
	int next_idx[NBODY*NBODY+1];
	int next_j[NBODY*NBODY+1];
	int Nphys;
	int NcrossCaus;						// used for collide2


	int next_src_idx;
	int prev_src_idx;
	bool skip;
	f_T Q;

	int changed;

	f_T error_interval;
	f_T area_interval;


};

// 改变形状也总是需要loc_centre，和Path相关函数有依赖，因为在源移动时只有这个参数需要给定
template < typename f_T >
struct src_shape_t{
	f_T rho;
	complex_t<f_T> loc_centre;
	_Host_Device src_shape_t();
	_Host_Device src_shape_t(const f_T Rho,const complex_t<f_T>& centre);
	_Host_Device src_shape_t(const src_shape_t<f_T> & another);
};



template < typename f_T >
struct src_ext_t {
	// complex_t<f_T> loc_centre;
	// f_T rho;
	src_shape_t<f_T> shape;
	f_T src_area;
	f_T Mag;
	// f_T Mag_p;
	f_T Err;
	f_T Area;
	f_T Err_A;
	f_T Mag_P;
	f_T Err_Mag_P;
	bool SolveSucceed;
	bool Break;
	// bool PointApproximation=0;
	complex_t<f_T> roots[LENGTH-1];
	int points_used;	
	#if DetailedRecord
		bool phys_tst[5];
		f_T phys_err[5];
	#endif

	int Ncross;
	int idx_cross[NCROSS_CANDIDATE];
	int j_cross[NCROSS_CANDIDATE];
	bool additional_cross[NCROSS_CANDIDATE];


	// used for collide2
	int NBuried;
	int idx_buried[NCROSS_CANDIDATE];
	f_T Qmin_buried[NCROSS_CANDIDATE];
	f_T Qmax_buried[NCROSS_CANDIDATE];
	


	src_pt_t<f_T> margin_pts[NPOINTS];
	// f_T err_rank[NPOINTS];
	_Host_Device src_ext_t();
};

// alpha 是src运动方向和正方向的夹角
// u_0 是原点到src路径垂足的距离
// 原点到垂足的指向是alpha-PI/2
template < typename f_T >
struct src_params_t {
	f_T t_0;
	f_T t_E;
//	complex_t<f_T> u_0;
	f_T u_0;
	f_T alpha;
//	f_T time;
	src_shape_t<f_T> shape;
	f_T s;
	f_T q;

	// f_T log_prob;

	_Host_Device src_params_t();
	_Host_Device src_params_t(const src_params_t<f_T>& another);


	// _Host_Device f_T& operator[](const int idx);

	_Host_Device void operator=(const src_params_t<f_T>& another);


};



template <typename f_T>
struct params_back_t
{
	complex_t<f_T> zeta;
	f_T s;
	f_T q;
	f_T rho;
	_Host_Device params_back_t();
};


template <typename f_T>
struct CC_t
{
	f_T psi;
	complex_t<f_T> crit[NBODY*NBODY];
	complex_t<f_T> caus[NBODY*NBODY];
	int next_idx;
	int next_j[NBODY*NBODY];
	bool skip;

};