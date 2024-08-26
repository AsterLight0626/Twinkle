#include"StructSrc.h"




template < typename f_T >
_Host_Device image_pt_t<f_T>::image_pt_t()
{
    // position = complex_t<f_T>(0,0);
}
template image_pt_t<float>::image_pt_t();
template image_pt_t<double>::image_pt_t();

template < typename f_T >
_Host_Device image_pt_t<f_T>::image_pt_t(const image_pt_t<f_T>& other)
{
    this->position = other.position;
    this->physical = other.physical;
    this->parity = other.parity;
}
template image_pt_t<float>::image_pt_t(const image_pt_t<float>& other);
template image_pt_t<double>::image_pt_t(const image_pt_t<double>& other);

template < typename f_T >
_Host_Device src_ext_t<f_T>::src_ext_t()
{
    // this->shape = src_shape_t<f_T>();
}
template _Host_Device src_ext_t<float>::src_ext_t();
template _Host_Device src_ext_t<double>::src_ext_t();

template <typename f_T>
_Host_Device src_shape_t<f_T>::src_shape_t()
{
	// this->rho = 0;
	// this->loc_centre = complex_t<f_T>(0.,0.);
};
template _Host_Device src_shape_t<float>::src_shape_t();
template _Host_Device src_shape_t<double>::src_shape_t();

template <typename f_T>
_Host_Device src_shape_t<f_T>::src_shape_t(const f_T Rho,const complex_t<f_T>& centre)
{
	this->rho = Rho;
	this->loc_centre = centre;
};
template _Host_Device src_shape_t<float>::src_shape_t(const float Rho, const complex_t<float>& centre);
template _Host_Device src_shape_t<double>::src_shape_t(const double Rho, const complex_t<double>& centre);

template <typename f_T>
_Host_Device src_shape_t<f_T>::src_shape_t(const src_shape_t<f_T> & another)
{
    this->rho = another.rho;
    this->loc_centre = another.loc_centre;
}
template _Host_Device src_shape_t<float>::src_shape_t(const src_shape_t<float> & another);
template _Host_Device src_shape_t<double>::src_shape_t(const src_shape_t<double> & another);


template < typename f_T >
_Host_Device src_params_t<f_T>::src_params_t()
{
    // this->t_0 = 0;
    // this->t_E = 1;
    // this->u_0 = 1;
    // this->alpha = 0;
    // // this->time = 0;
    // this->shape = src_shape_t<f_T>();
    // this->log_prob=1.;
}
template _Host_Device src_params_t<float>::src_params_t();
template _Host_Device src_params_t<double>::src_params_t();

template < typename f_T >
_Host_Device src_params_t<f_T>::src_params_t(const src_params_t<f_T>& another)
{
    this->t_0 = another.t_0;
    this->t_E = another.t_E;
    this->u_0 = another.u_0;
    this->alpha = another.alpha;
    this->shape = another.shape;
    this->s = another.s;
    this->q = another.q;
}
template _Host_Device src_params_t<float>::src_params_t(const src_params_t<float>& another);
template _Host_Device src_params_t<double>::src_params_t(const src_params_t<double>& another);



template <typename f_T>
_Host_Device void src_params_t<f_T>::operator=(const src_params_t<f_T>& another)
{
    this->t_0 = another.t_0;
    this->t_E = another.t_E;
    this->u_0 = another.u_0;
    this->alpha = another.alpha;
    this->shape = another.shape;
    this->s = another.s;
    this->q = another.q;
}
template _Host_Device void src_params_t<float>::operator=(const src_params_t<float>& another);
template _Host_Device void src_params_t<double>::operator=(const src_params_t<double>& another);



template <typename f_T>
_Host_Device params_back_t<f_T>::params_back_t()
{
    //
}
template _Host_Device params_back_t<float>::params_back_t();
template _Host_Device params_back_t<double>::params_back_t();


