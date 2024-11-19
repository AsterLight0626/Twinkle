#include"StructComplex.h"

// template < isFloating f_T >
// _Host_Device complex_t<f_T>::complex_t(f_T _re, f_T _im) : re(_re), im(_im)   {};
// template _Host_Device complex_t<float>::complex_t(float _re, float _im);
// template _Host_Device complex_t<double>::complex_t(double _re, double _im);

// template < isFloating f_T >
// _Host_Device complex_t<f_T>::complex_t(f_T _re)          : complex_t(_re, 0.) {};
// template _Host_Device complex_t<float>::complex_t(float _re);
// template _Host_Device complex_t<double>::complex_t(double _re);

// template < isFloating f_T >
// _Host_Device complex_t<f_T>::complex_t()                 : complex_t(0. , 0.) {};
// template _Host_Device complex_t<float>::complex_t();
// template _Host_Device complex_t<double>::complex_t();

// template < isFloating f_T >
// template<isFloating f2_T>
// _Host_Device complex_t<f_T>::complex_t(const complex_t<f2_T> & z)
// {
//     this->re = (f_T) z.re;
//     this->im = (f_T) z.im;
// }
// template _Host_Device complex_t<float>::complex_t(const complex_t<double> & z);
// template _Host_Device complex_t<double>::complex_t(const complex_t<float> & z);
// template _Host_Device complex_t<double>::complex_t(const complex_t<double> & z);
// template _Host_Device complex_t<float>::complex_t(const complex_t<float> & z);



template <isFloating f_T>
template <isFloating f2_T>
_Host_Device void complex_t<f_T>::operator=(const complex_t<f2_T> & z)
{
    this->re = (f_T) z.re;      // 转化为f_T类型
    this->im = (f_T) z.im;
}
template _Host_Device void complex_t<float>::operator=<double>(const complex_t<double> & z);
template _Host_Device void complex_t<double>::operator=<float>(const complex_t<float> & z);
template _Host_Device void complex_t<double>::operator=<double>(const complex_t<double> & z);
template _Host_Device void complex_t<float>::operator=<float>(const complex_t<float> & z);

// template <isFloating f_T>
// template <isFloating f2_T>
// _Host_Device void complex_t<f_T>::operator=(const f2_T & z)
// {
//     this->re = (f_T) z;
//     this->im = (f_T) 0;
// }
// template _Host_Device void complex_t<float>::operator=<float>(const float & z);
// template _Host_Device void complex_t<float>::operator=<double>(const double & z);
// template _Host_Device void complex_t<double>::operator=<double>(const double & z);
// template _Host_Device void complex_t<double>::operator=<float>(const float & z);





// 函数声明：复数运算

// Unary

template <isFloating f_T>
_Host_Device f_T norm(const complex_t<f_T> &z)
{
    return z.re*z.re + z.im*z.im;
};
template _Host_Device float norm(const complex_t<float> &z);
template _Host_Device double norm(const complex_t<double> &z);

// real 与 imag 直接通过 .re .im访问
// 好吧，也写一份
template <isFloating f_T>
_Host_Device f_T real(const complex_t<f_T> &z)
{
    return z.re;
}
template _Host_Device float real(const complex_t<float> &z);
template _Host_Device double real(const complex_t<double> &z);

template <isFloating f_T>
_Host_Device f_T imag(const complex_t<f_T> &z)
{
    return z.im;
}
template _Host_Device float imag(const complex_t<float> &z);
template _Host_Device double imag(const complex_t<double> &z);
 

template <isFloating f_T>
_Host_Device complex_t<f_T> conj(const complex_t<f_T> &z)
{
    return complex_t<f_T>(z.re,-z.im);
}
template _Host_Device complex_t<float> conj(const complex_t<float> &z);
template _Host_Device complex_t<double> conj(const complex_t<double> &z);

// 负号
template<isFloating f_T>
_Host_Device complex_t<f_T> operator-(const complex_t<f_T> z)
{
    complex_t<f_T> res;
    res.re = -z.re;
    res.im = -z.im;
    return res;
}
template _Host_Device complex_t<float> operator-<float>(const complex_t<float> z);
template _Host_Device complex_t<double> operator-<double>(const complex_t<double> z);

// 1模数
// template <isFloating f_T>
// _Host_Device f_T abs_c(const complex_t<f_T> &z)
// {
//     return abs(z.re)+abs(z.im);
// };
// template <isFloating f_T>
// _Host f_T abs_c(const complex_t<f_T> &z)
// {
//     return std::sqrt(norm(z));
// }
// template _Host float abs_c(const complex_t<float> &z);
// template _Host double abs_c(const complex_t<double> &z);

// _Host float abs_c(const complex_t<float> &z)
// {
//     return std::sqrt(norm(z));
// }
// _Host double abs_c(const complex_t<double> &z)
// {
//     return std::sqrt(norm(z));
// }

_Device float abs_c(const complex_t<float> &z)
{
    return __sqrt2(norm(z));
}
_Device double abs_c(const complex_t<double> &z)
{
    return __sqrt1(norm(z));
}


// template <isFloating f_T>
// _Host complex_t<f_T> sqrt(const complex_t<f_T> &x)
// {
// 	f_T prev_length = std::sqrt(x.re*x.re + x.im*x.im);             // norm before sqrt
//     f_T after_length = std::sqrt(prev_length);                      // norm after sqrt
//     complex_t<f_T> z;
//     if (prev_length>0) {
//         z.re = after_length * (std::sqrt(abs((1 + x.re / prev_length) / 2))*((x.im>=0) ? 1 : -1));
//         z.im = after_length * std::sqrt(abs((1 - x.re / prev_length) / 2));
//     }
//     return z;
// }
// _Host complex_t<float> sqrt(const complex_t<float> &x)
// {
// 	float prev_length = std::sqrt(x.re*x.re + x.im*x.im);             // norm before sqrt
//     float after_length = std::sqrt(prev_length);                      // norm after sqrt
//     complex_t<float> z;
//     if (prev_length>0) {
//         z.re = after_length * (std::sqrt(abs((1 + x.re / prev_length) / 2))*((x.im>=0) ? 1 : -1));
//         z.im = after_length * std::sqrt(abs((1 - x.re / prev_length) / 2));
//     }
//     return z;
// }
// _Host complex_t<double> sqrt(const complex_t<double> &x)
// {
// 	double prev_length = std::sqrt(x.re*x.re + x.im*x.im);             // norm before sqrt
//     double after_length = std::sqrt(prev_length);                      // norm after sqrt
//     complex_t<double> z;
//     if (prev_length>0) {
//         z.re = after_length * (std::sqrt(abs((1 + x.re / prev_length) / 2))*((x.im>=0) ? 1 : -1));
//         z.im = after_length * std::sqrt(abs((1 - x.re / prev_length) / 2));
//     }
//     return z;
// }
_Device complex_t<float> sqrt(const complex_t<float> &x)
{
	float prev_length = __sqrt2(x.re*x.re + x.im*x.im);             // norm before sqrt
    float after_length = __sqrt2(prev_length);                      // norm after sqrt
    complex_t<float> z;
    if (prev_length>0) {
        z.re = after_length * (__sqrt2(fabs((1 + x.re / prev_length) / 2))*((x.im>=0) ? 1 : -1));
        z.im = after_length * __sqrt2(fabs((1 - x.re / prev_length) / 2));
    }
    return z;
}
_Device complex_t<double> sqrt(const complex_t<double> &x)
{
	double prev_length = __sqrt1(x.re*x.re + x.im*x.im);             // norm before sqrt
    double after_length = __sqrt1(prev_length);                      // norm after sqrt
    complex_t<double> z;
    if (prev_length>0) {
        z.re = after_length * (__sqrt1(fabs((1 + x.re / prev_length) / 2))*((x.im>=0) ? 1 : -1));
        z.im = after_length * __sqrt1(fabs((1 - x.re / prev_length) / 2));
    }
    return z;
}

// 支持cout输出
template <isFloating f_T>
std::ostream & operator<<(std::ostream &out, const complex_t<f_T> &A)
{
    out <<"("<< A.re <<","<< A.im <<")";
    return out;
}
template std::ostream & operator<<<float>(std::ostream &out, const complex_t<float> &A);
template std::ostream & operator<<<double>(std::ostream &out, const complex_t<double> &A);



template <isFloating f_T>
_Host_Device bool operator==(const complex_t<f_T> &z1, const complex_t<f_T> &z2)
{
    if(z1.re==z2.re && z1.im==z2.im){return true;}
    else{return false;}
}
template _Host_Device bool operator==<float>(const complex_t<float> &z1, const complex_t<float> &z2);
template _Host_Device bool operator==<double>(const complex_t<double> &z1, const complex_t<double> &z2);


template <isFloating f_T>
_Host_Device bool operator!=(const complex_t<f_T> &z1, const complex_t<f_T> &z2)
{
    if(z1.re!=z2.re || z1.im!=z2.im){return true;}
    else{return false;}
}
template _Host_Device bool operator!=<float>(const complex_t<float> &z1, const complex_t<float> &z2);
template _Host_Device bool operator!=<double>(const complex_t<double> &z1, const complex_t<double> &z2);

// same type
template <isFloating f_T>
_Host_Device complex_t<f_T> operator+(const complex_t<f_T> &z1,const complex_t<f_T> &z2)
{
    complex_t<f_T> res;
    res.re = z1.re + z2.re;
    res.im = z1.im + z2.im;
    return res;
}
template _Host_Device complex_t<float> operator+<float>(const complex_t<float> &z1,const complex_t<float> &z2);
template _Host_Device complex_t<double> operator+<double>(const complex_t<double> &z1,const complex_t<double> &z2);

template <isFloating f_T>
_Host_Device complex_t<f_T> operator-(const complex_t<f_T> &z1,const complex_t<f_T> &z2)
{
    complex_t<f_T> res;
    res.re = z1.re - z2.re;
    res.im = z1.im - z2.im;
    return res;
}
template _Host_Device complex_t<float> operator-<float>(const complex_t<float> &z1,const complex_t<float> &z2);
template _Host_Device complex_t<double> operator-<double>(const complex_t<double> &z1,const complex_t<double> &z2);

template <isFloating f_T>
_Host_Device complex_t<f_T> operator*(const complex_t<f_T> &z1,const complex_t<f_T> &z2)
{
    complex_t<f_T> res;
    res.re = z1.re * z2.re - z1.im * z2.im;
    res.im = z1.im * z2.re + z1.re * z2.im;
    return res;
}
template _Host_Device complex_t<float> operator*<float>(const complex_t<float> &z1,const complex_t<float> &z2);
template _Host_Device complex_t<double> operator*<double>(const complex_t<double> &z1,const complex_t<double> &z2);

// _Host_Device complex_t<double> operator*(const complex_t<double> &z1,const complex_t<double> &z2)
// {
//     complex_t<double> res;
//     double temp;
//     // res.re = z1.re * z2.re - z1.im * z2.im;
//     // res.im = z1.im * z2.re + z1.re * z2.im;
//     res.re = __dmul_rn(z1.re,z2.re);
//     temp = __dmul_rn(z1.im,z2.im);
//     res.re = __dsub_rn(res.re,temp);

//     res.im = __dmul_rn(z1.im,z2.re);
//     temp = __dmul_rn(z1.re,z2.im);
//     res.im = __dadd_rn(res.im,temp);
    
//     return res;
// }


template <isFloating f_T>
_Host_Device complex_t<f_T> operator/(const complex_t<f_T> &z1,const complex_t<f_T> &z2)
{
    f_T norm2;
    complex_t<f_T> res;
    norm2 = norm(z2);
    res.re = (z1.re * z2.re + z1.im * z2.im)/norm2;
    res.im = (z1.im * z2.re - z1.re * z2.im)/norm2;
    return res;
}
template _Host_Device complex_t<float> operator/<float>(const complex_t<float> &z1,const complex_t<float> &z2);
template _Host_Device complex_t<double> operator/<double>(const complex_t<double> &z1,const complex_t<double> &z2);

template <isFloating f_T>
_Host_Device void operator+=(complex_t<f_T> &z1,const complex_t<f_T> &z2)
{
    z1.re = z1.re + z2.re;
    z1.im = z1.im + z2.im;
}
template _Host_Device void operator+=<float>(complex_t<float> &z1,const complex_t<float> &z2);
template _Host_Device void operator+=<double>(complex_t<double> &z1,const complex_t<double> &z2);


// Real & Complex 只有同精度的才可以直接四则运算
template <isFloating f_T>
_Host_Device complex_t<f_T> operator+(const f_T &f1,const complex_t<f_T> &z2)
{
    complex_t<f_T> res=z2;
    res.re = res.re + f1;
    return res;
}
template _Host_Device complex_t<float> operator+<float>(const float &f1,const complex_t<float> &z2);
template _Host_Device complex_t<double> operator+<double>(const double &f1,const complex_t<double> &z2);

template <isFloating f_T>
_Host_Device complex_t<f_T> operator+(const complex_t<f_T> &z1,const f_T &f2)
{
    complex_t<f_T> res=z1;
    res.re = res.re + f2;
    return res;
}
template _Host_Device complex_t<float> operator+<float>(const complex_t<float> &z1,const float &f2);
template _Host_Device complex_t<double> operator+<double>(const complex_t<double> &z1,const double &f2);


template <isFloating f_T>
_Host_Device complex_t<f_T> operator-(const f_T &f1,const complex_t<f_T> &z2)
{
    complex_t<f_T> res=z2;
    res.re = f1 - res.re;
    res.im = -res.im;
    return res;
}
template _Host_Device complex_t<float> operator-<float>(const float &f1,const complex_t<float> &z2);
template _Host_Device complex_t<double> operator-<double>(const double &f1,const complex_t<double> &z2);

template <isFloating f_T>
_Host_Device complex_t<f_T> operator-(const complex_t<f_T> &z1,const f_T &f2)
{
    complex_t<f_T> res=z1;
    res.re = res.re - f2;
    return res;
}
template _Host_Device complex_t<float> operator-<float>(const complex_t<float> &z1,const float &f2);
template _Host_Device complex_t<double> operator-<double>(const complex_t<double> &z1,const double &f2);

template <isFloating f_T>
_Host_Device complex_t<f_T> operator*(const f_T &f1,const complex_t<f_T> &z2)
{
    complex_t<f_T> res=z2;
    res.re *= f1;
    res.im *= f1;
    return res;
}
template _Host_Device complex_t<float> operator*<float>(const float &f1,const complex_t<float> &z2);
template _Host_Device complex_t<double> operator*<double>(const double &f1,const complex_t<double> &z2);

template <isFloating f_T>
_Host_Device complex_t<f_T> operator*(const complex_t<f_T> &z1,const f_T &f2)
{
    complex_t<f_T> res=z1;
    res.re *= f2;
    res.im *= f2;
    return res;
}
template _Host_Device complex_t<float> operator*<float>(const complex_t<float> &z1,const float &f2);
template _Host_Device complex_t<double> operator*<double>(const complex_t<double> &z1,const double &f2);

// _Host_Device complex_t<double> operator*(const complex_t<double> &z1,const double &f2)
// {
//     complex_t<double> res=z1;
//     res.re = __dmul_rn(res.re,f2);
//     res.im = __dmul_rn(res.im,f2);
//     return res;
// }


template <isFloating f_T>
_Host_Device complex_t<f_T> operator/(const f_T &f1,const complex_t<f_T> &z2)
{
    complex_t<f_T> res; //=conj(z2);
    f_T coeff = f1 / norm(z2);
    // f_T coeff = f1 / (z2.re*z2.re + z2.im*z2.im);
    res.re = z2.re*coeff;
    res.im = -z2.im*coeff;
    return res;
}
template _Host_Device complex_t<float> operator/<float>(const float &f1,const complex_t<float> &z2);
template _Host_Device complex_t<double> operator/<double>(const double &f1,const complex_t<double> &z2);

template <isFloating f_T>
_Host_Device complex_t<f_T> operator/(const complex_t<f_T> &z1,const f_T &f2)
{
    complex_t<f_T> res=z1;
    res.re /= f2;
    res.im /= f2;
    return res;
}
template _Host_Device complex_t<float> operator/<float>(const complex_t<float> &z1,const float &f2);
template _Host_Device complex_t<double> operator/<double>(const complex_t<double> &z1,const double &f2);

// template <typename f_T>
// _Host_Device f_T abs(const f_T &f)
// {
//     if(f<0){return -f;}
//     else{return f;}
// }
// template _Host_Device float abs(const float &f);
// template _Host_Device double abs(const double &f);
// // template _Host_Device int abs(const int &f);

template <typename f_T>
_Host_Device f_T min(const f_T &f1, const f_T &f2)
{
    if(f1<f2){return f1;}
    else{return f2;}
}
template _Host_Device float min(const float &f1, const float &f2);
template _Host_Device double min(const double &f1, const double &f2);
template _Host_Device int min(const int &f1, const int &f2);

template <typename f_T>
_Host_Device f_T max(const f_T &f1, const f_T &f2)
{
    if(f1>f2){return f1;}
    else{return f2;}
}
template _Host_Device float max(const float &f1, const float &f2);
template _Host_Device double max(const double &f1, const double &f2);
template _Host_Device int max(const int &f1, const int &f2);



template<typename f_T>
_Host_Device void ADD(complex_t<f_T> &z, const complex_t<f_T> &x, const complex_t<f_T> &y)
{
    z.re = x.re + y.re;
    z.im = x.im + y.im;
}
template _Host_Device void ADD(complex_t<float> &z, const complex_t<float> &x, const complex_t<float> &y);
template _Host_Device void ADD(complex_t<double> &z, const complex_t<double> &x, const complex_t<double> &y);

template<typename f_T>
_Host_Device void SUB(complex_t<f_T> &z, const complex_t<f_T> &x, const complex_t<f_T> &y)
{
    z.re = x.re - y.re;
    z.im = x.im - y.im;
}
template _Host_Device void SUB(complex_t<float> &z, const complex_t<float> &x, const complex_t<float> &y);
template _Host_Device void SUB(complex_t<double> &z, const complex_t<double> &x, const complex_t<double> &y);



template<typename f_T>
_Host_Device void MUL(complex_t<f_T> &z, const complex_t<f_T> &x, const complex_t<f_T> &y)
{
    // f_T temp; // in case of MUL(a,a,a)
    // f_T temp1;



    // // temp = __dmul_rn(x.re,y.re);
    // // temp1 = __dmul_rn(x.im,y.im);
    // // temp = __dsub_rn(temp,temp1);
    // temp = x.re * y.re;
    // temp1 = x.im * y.im;
    // temp = temp - temp1;

    // // z.im = __dmul_rn(x.im,y.re);
    // // temp1 = __dmul_rn(x.re,y.im);
    // // z.im = __dadd_rn(z.im,temp1);
    // z.im = x.im*y.re;
    // temp1 = x.re * y.im;
    // z.im = z.im + temp1;

    // // temp = x.re*y.re - x.im*y.im;
    // // z.im = x.re*y.im + x.im*y.re;
    // z.re = temp;
    f_T temp; // in case of MUL(a,a,a)
    temp = x.re*y.re - x.im*y.im;
    z.im = x.re*y.im + x.im*y.re;
    z.re = temp;

}
template _Host_Device void MUL(complex_t<float> &z, const complex_t<float> &x, const complex_t<float> &y);
template _Host_Device void MUL(complex_t<double> &z, const complex_t<double> &x, const complex_t<double> &y);

template<typename f_T>
_Host_Device void DIV(complex_t<f_T> &z, const complex_t<f_T> &x, const complex_t<f_T> &y)
{
    f_T md;// = __dadd_rn(__dmul_rn(y.re,y.re) , __dmul_rn(y.im,y.im));
    md = y.re * y.re + y.im * y.im;
    f_T temp; // in case of MUL(a,a,a)
    // temp = __ddiv_rn(__dadd_rn(__dmul_rn(x.re,y.re) , __dmul_rn(x.im,y.im)) , md);
    // z.im = __ddiv_rn(__dsub_rn( __dmul_rn(x.im,y.re) , __dmul_rn(x.re,y.im)) , md);
    temp = (x.re*y.re + x.im * y.im) / md;
    z.im = (x.im*y.re - x.re * y.im) / md;
    z.re = temp;


}
template _Host_Device void DIV(complex_t<float> &z, const complex_t<float> &x, const complex_t<float> &y);
template _Host_Device void DIV(complex_t<double> &z, const complex_t<double> &x, const complex_t<double> &y);



_Host_Device complex_t<double> INV(const complex_t<double>& z)
{
    complex_t<double> res; //=conj(z2);
    double coeff = norm(z);
    res.re = z.re / coeff;
    res.im = -z.im / coeff;
    return res;
}

template<typename f_T>
_Host_Device void SUB(complex_t<f_T> &z, const f_T &x, const complex_t<f_T> &y)
{
    z.re = x - y.re;
    z.im = - y.im;
}
template _Host_Device void SUB(complex_t<float> &z, const float &x, const complex_t<float> &y);
template _Host_Device void SUB(complex_t<double> &z, const double &x, const complex_t<double> &y);

template<typename f_T>
_Host_Device void conj(complex_t<f_T> &z, const complex_t<f_T> &x)
{
    z.re = x.re;
    z.im = -x.im;
}
template _Host_Device void conj(complex_t<float> &z, const complex_t<float> &x);
template _Host_Device void conj(complex_t<double> &z, const complex_t<double> &x);


template<typename f_T>
_Host_Device void MUL(complex_t<f_T> &z, const f_T &x, const complex_t<f_T> &y)
{
    z.re = x * y.re;
    z.im = x * y.im;
}
template _Host_Device void MUL(complex_t<float> &z, const float &x, const complex_t<float> &y);
template _Host_Device void MUL(complex_t<double> &z, const double &x, const complex_t<double> &y);


_Device void sqrt(complex_t<float> &z, const complex_t<float> &x)
{
	float prev_length = __sqrt2(x.re*x.re + x.im*x.im);             // norm before sqrt
    float after_length = __sqrt2(prev_length);                      // norm after sqrt

    if (prev_length>0) {
        z.re = after_length * (__sqrt2(fabs((1 + x.re / prev_length) / 2))*((x.im>=0) ? 1 : -1));
        z.im = after_length * __sqrt2(fabs((1 - x.re / prev_length) / 2));
    }
}

_Device void sqrt(complex_t<double> &z, const complex_t<double> &x)
{
	double prev_length = __sqrt1(x.re*x.re + x.im*x.im);             // norm before sqrt
    double after_length = __sqrt1(prev_length);                      // norm after sqrt

    if (prev_length>0) {
        z.re = after_length * (__sqrt1(fabs((1 + x.re / prev_length) / 2))*((x.im>=0) ? 1 : -1));
        z.im = after_length * __sqrt1(fabs((1 - x.re / prev_length) / 2));
    }  
}

_Host_Device void SUB(complex_t<double> &z, const complex_t<double> &zD, const complex_t<float> &zF)
{
	z.re = zD.re - zF.re; 
    z.im = zD.im - zF.im;    
}

_Host_Device void SUB(complex_t<double> &z, const complex_t<float> &zF, const complex_t<double> &zD)
{
	z.re = zF.re - zD.re;
    z.im = zF.im - zD.im;
}


_Host_Device void SUB(complex_t<float> &z, const complex_t<double> &zD1, const complex_t<double> &zD2)
{
	z.re = float(zD1.re) - float(zD2.re);
    z.im = float(zD1.im) - float(zD2.im);
}


_Host_Device void MUL(complex_t<double> &z, const complex_t<float> &x, const complex_t<double> &y)
{
    double temp; // in case of MUL(a,a,a)
    temp = x.re*y.re - x.im*y.im;
    z.im = x.re*y.im + x.im*y.re;
    z.re = temp;
}

template <isFloating f_T>
_Host_Device void operator-=(complex_t<f_T> &z1,const complex_t<f_T> &z2)
{
    z1.re = z1.re - z2.re;
    z1.im = z1.im - z2.im;
}
template _Host_Device void operator-=<float>(complex_t<float> &z1,const complex_t<float> &z2);
template _Host_Device void operator-=<double>(complex_t<double> &z1,const complex_t<double> &z2);

template<typename f_T>
_Host_Device f_T wedge_product(const complex_t<f_T> &z1, const complex_t<f_T> &z2)
{
    return z1.re*z2.im - z1.im*z2.re;
}
template _Host_Device float wedge_product(const complex_t<float> &z1, const complex_t<float> &z2);
template _Host_Device double wedge_product(const complex_t<double> &z1, const complex_t<double> &z2);