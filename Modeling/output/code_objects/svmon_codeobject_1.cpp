#include "code_objects/svmon_codeobject_1.h"
#include "objects.h"
#include "brianlib/common_math.h"
#include "brianlib/stdint_compat.h"
#include<cmath>
#include<ctime>
#include<iostream>
#include<fstream>
#include<climits>
#include "gsl_const_mksa.h"
#include "math.h"

////// SUPPORT CODE ///////
namespace {
        
    #include <math.h>
    double Hill(double x,double K,double n)
    {
        return pow(x,n)/(pow(x,n)+pow(K,n));
    };
    #include <gsl/gsl_const_mksa.h>
    const double R = GSL_CONST_MKSA_MOLAR_GAS; 
    const double F = GSL_CONST_MKSA_FARADAY;
    double ThermalVoltage(const double T)
    {
        return R*(T+273.15)/F;
    };
    double NernstPotential(double x_e,double x_i,double z,const double T)
    {
        return ThermalVoltage(T)*log(x_e/x_i)/z;
    };
    template <typename T>
    static inline T _clip(const T value, const double a_min, const double a_max)
    {
        if (value < a_min)
            return a_min;
        if (value > a_max)
            return a_max;
        return value;
    }
    template < typename T1, typename T2 > struct _higher_type;
    template < > struct _higher_type<int32_t,int32_t> { typedef int32_t type; };
    template < > struct _higher_type<int32_t,int64_t> { typedef int64_t type; };
    template < > struct _higher_type<int32_t,float> { typedef float type; };
    template < > struct _higher_type<int32_t,double> { typedef double type; };
    template < > struct _higher_type<int32_t,long double> { typedef long double type; };
    template < > struct _higher_type<int64_t,int32_t> { typedef int64_t type; };
    template < > struct _higher_type<int64_t,int64_t> { typedef int64_t type; };
    template < > struct _higher_type<int64_t,float> { typedef float type; };
    template < > struct _higher_type<int64_t,double> { typedef double type; };
    template < > struct _higher_type<int64_t,long double> { typedef long double type; };
    template < > struct _higher_type<float,int32_t> { typedef float type; };
    template < > struct _higher_type<float,int64_t> { typedef float type; };
    template < > struct _higher_type<float,float> { typedef float type; };
    template < > struct _higher_type<float,double> { typedef double type; };
    template < > struct _higher_type<float,long double> { typedef long double type; };
    template < > struct _higher_type<double,int32_t> { typedef double type; };
    template < > struct _higher_type<double,int64_t> { typedef double type; };
    template < > struct _higher_type<double,float> { typedef double type; };
    template < > struct _higher_type<double,double> { typedef double type; };
    template < > struct _higher_type<double,long double> { typedef long double type; };
    template < > struct _higher_type<long double,int32_t> { typedef long double type; };
    template < > struct _higher_type<long double,int64_t> { typedef long double type; };
    template < > struct _higher_type<long double,float> { typedef long double type; };
    template < > struct _higher_type<long double,double> { typedef long double type; };
    template < > struct _higher_type<long double,long double> { typedef long double type; };
    // General template, used for floating point types
    template < typename T1, typename T2 >
    static inline typename _higher_type<T1,T2>::type
    _brian_mod(T1 x, T2 y)
    {
        return x-y*floor(1.0*x/y);
    }
    // Specific implementations for integer types
    // (from Cython, see LICENSE file)
    template <>
    inline int32_t _brian_mod(int32_t x, int32_t y)
    {
        int32_t r = x % y;
        r += ((r != 0) & ((r ^ y) < 0)) * y;
        return r;
    }
    template <>
    inline int64_t _brian_mod(int32_t x, int64_t y)
    {
        int64_t r = x % y;
        r += ((r != 0) & ((r ^ y) < 0)) * y;
        return r;
    }
    template <>
    inline int64_t _brian_mod(int64_t x, int32_t y)
    {
        int64_t r = x % y;
        r += ((r != 0) & ((r ^ y) < 0)) * y;
        return r;
    }
    template <>
    inline int64_t _brian_mod(int64_t x, int64_t y)
    {
        int64_t r = x % y;
        r += ((r != 0) & ((r ^ y) < 0)) * y;
        return r;
    }
    // General implementation, used for floating point types
    template < typename T1, typename T2 >
    static inline typename _higher_type<T1,T2>::type
    _brian_floordiv(T1 x, T2 y)
    {{
        return floor(1.0*x/y);
    }}
    // Specific implementations for integer types
    // (from Cython, see LICENSE file)
    template <>
    inline int32_t _brian_floordiv<int32_t, int32_t>(int32_t a, int32_t b) {
        int32_t q = a / b;
        int32_t r = a - q*b;
        q -= ((r != 0) & ((r ^ b) < 0));
        return q;
    }
    template <>
    inline int64_t _brian_floordiv<int32_t, int64_t>(int32_t a, int64_t b) {
        int64_t q = a / b;
        int64_t r = a - q*b;
        q -= ((r != 0) & ((r ^ b) < 0));
        return q;
    }
    template <>
    inline int64_t _brian_floordiv<int64_t, int>(int64_t a, int32_t b) {
        int64_t q = a / b;
        int64_t r = a - q*b;
        q -= ((r != 0) & ((r ^ b) < 0));
        return q;
    }
    template <>
    inline int64_t _brian_floordiv<int64_t, int64_t>(int64_t a, int64_t b) {
        int64_t q = a / b;
        int64_t r = a - q*b;
        q -= ((r != 0) & ((r ^ b) < 0));
        return q;
    }
    #ifdef _MSC_VER
    #define _brian_pow(x, y) (pow((double)(x), (y)))
    #else
    #define _brian_pow(x, y) (pow((x), (y)))
    #endif

}

////// HASH DEFINES ///////



void _run_svmon_codeobject_1()
{
    using namespace brian;


    ///// CONSTANTS ///////////
    const double C0_Cl_E = 130.0;
const double C0_Cl_N = 5.0;
const double C0_K_E = 3.0;
const double C0_K_N = 130.0;
const double C0_Na_E = 145.0;
const double C0_Na_N = 10.0;
const double I_NKP_max = 0.07071461769188153;
const size_t _numN = 1;
const int64_t T_exp = 37;
const double V_E = 5e-16;
const double V_N = 5e-16;
const double V_T = 0.003188413306993716;
const size_t _num_____source_I_NKP_hh_1_f_NaK_hh_1_sigma_hh_1_C_Na_E_hh_1_n_Na_E = 1;
const size_t _num____source_I_Cl_L_hh_1_E_Cl_hh_1_C_Cl_E_hh_1_n_Cl_E = 1;
const size_t _num____source_I_Cl_L_hh_1_E_Cl_hh_1_C_Cl_N_hh_1_n_Cl_N = 1;
const size_t _num____source_I_KCC_hh_1_E_Cl_hh_1_C_Cl_E_hh_1_n_Cl_E = 1;
const size_t _num____source_I_KCC_hh_1_E_Cl_hh_1_C_Cl_N_hh_1_n_Cl_N = 1;
const size_t _num____source_I_KCC_hh_1_E_K_hh_1_C_K_E_hh_1_n_K_E = 1;
const size_t _num____source_I_KCC_hh_1_E_K_hh_1_C_K_N_hh_1_n_K_N = 1;
const size_t _num____source_I_K_hh_1_E_K_hh_1_C_K_E_hh_1_n_K_E = 1;
const size_t _num____source_I_K_hh_1_E_K_hh_1_C_K_N_hh_1_n_K_N = 1;
const size_t _num____source_I_Na_L_hh_1_E_Na_hh_1_C_Na_E_hh_1_n_Na_E = 1;
const size_t _num____source_I_Na_L_hh_1_E_Na_hh_1_C_Na_N_hh_1_n_Na_N = 1;
const size_t _num____source_I_Na_hh_1_E_Na_hh_1_C_Na_E_hh_1_n_Na_E = 1;
const size_t _num____source_I_Na_hh_1_E_Na_hh_1_C_Na_N_hh_1_n_Na_N = 1;
const size_t _num___source_E_Cl_hh_1_C_Cl_E_hh_1_n_Cl_E = 1;
const size_t _num___source_E_Cl_hh_1_C_Cl_N_hh_1_n_Cl_N = 1;
const size_t _num___source_E_K_hh_1_C_K_E_hh_1_n_K_E = 1;
const size_t _num___source_E_K_hh_1_C_K_N_hh_1_n_K_N = 1;
const size_t _num___source_E_Na_hh_1_C_Na_E_hh_1_n_Na_E = 1;
const size_t _num___source_E_Na_hh_1_C_Na_N_hh_1_n_Na_N = 1;
const size_t _num___source_I_NKP_hh_1_C_K_E_hh_1_n_K_E = 1;
const size_t _num___source_I_NKP_hh_1_C_Na_N_hh_1_n_Na_N = 1;
const size_t _num___source_I_NKP_hh_1_f_NaK_hh_1_v = 1;
const size_t _num__source_C_Cl_N_hh_1_n_Cl_N = 1;
const size_t _num__source_C_K_N_hh_1_n_K_N = 1;
const size_t _num__source_C_Na_N_hh_1_n_Na_N = 1;
const size_t _num__source_I_Cl_L_hh_1_v = 1;
const size_t _num__source_I_K_hh_1_n = 1;
const size_t _num__source_I_K_hh_1_v = 1;
const size_t _num__source_I_Na_L_hh_1_v = 1;
const size_t _num__source_I_Na_hh_1_h = 1;
const size_t _num__source_I_Na_hh_1_m = 1;
const size_t _num__source_I_Na_hh_1_v = 1;
const size_t _num_clock_t = 1;
const size_t _num_indices = 1;
const size_t _num_source_v = 1;
const double g_Cl_L = 0.49999999999999994;
const double g_K = 6929.999999999999;
const double g_KCC = - 0.0454665638343149;
const double g_Na = 20400.0;
const double g_Na_L = 0.3270276859728401;
const double mM = true;
double* const _array_svmon_t = _dynamic_array_svmon_t.empty()? 0 : &_dynamic_array_svmon_t[0];
const size_t _numt = _dynamic_array_svmon_t.size();
const double volt = true;
const double zeta_K = 0.2;
const double zeta_Na = 13.0;
    ///// POINTERS ////////////
        
    int32_t*   _ptr_array_svmon_N = _array_svmon_N;
    double*   _ptr_array_hh_1_n_Na_E = _array_hh_1_n_Na_E;
    double*   _ptr_array_hh_1_n_Cl_E = _array_hh_1_n_Cl_E;
    double*   _ptr_array_hh_1_n_Cl_N = _array_hh_1_n_Cl_N;
    double*   _ptr_array_hh_1_n_K_E = _array_hh_1_n_K_E;
    double*   _ptr_array_hh_1_n_K_N = _array_hh_1_n_K_N;
    double*   _ptr_array_hh_1_n_Na_N = _array_hh_1_n_Na_N;
    double*   _ptr_array_hh_1_v = _array_hh_1_v;
    double*   _ptr_array_hh_1_n = _array_hh_1_n;
    double*   _ptr_array_hh_1_h = _array_hh_1_h;
    double*   _ptr_array_hh_1_m = _array_hh_1_m;
    double*   _ptr_array_svmon_clock_1_t = _array_svmon_clock_1_t;
    int32_t*   _ptr_array_svmon__indices = _array_svmon__indices;
    double* __restrict  _ptr_array_svmon_t = _array_svmon_t;


    _dynamic_array_svmon_t.push_back(_ptr_array_svmon_clock_1_t[0]);

    const size_t _new_size = _dynamic_array_svmon_t.size();
    // Resize the dynamic arrays
    _dynamic_array_svmon_C_Cl_N.resize(_new_size, _num_indices);
    _dynamic_array_svmon_C_K_N.resize(_new_size, _num_indices);
    _dynamic_array_svmon_C_Na_N.resize(_new_size, _num_indices);
    _dynamic_array_svmon_E_Cl.resize(_new_size, _num_indices);
    _dynamic_array_svmon_E_K.resize(_new_size, _num_indices);
    _dynamic_array_svmon_E_Na.resize(_new_size, _num_indices);
    _dynamic_array_svmon_I_Cl_L.resize(_new_size, _num_indices);
    _dynamic_array_svmon_I_K.resize(_new_size, _num_indices);
    _dynamic_array_svmon_I_KCC.resize(_new_size, _num_indices);
    _dynamic_array_svmon_I_Na.resize(_new_size, _num_indices);
    _dynamic_array_svmon_I_Na_L.resize(_new_size, _num_indices);
    _dynamic_array_svmon_I_NKP.resize(_new_size, _num_indices);
    _dynamic_array_svmon_v.resize(_new_size, _num_indices);

    // scalar code
    const size_t _vectorisation_idx = -1;
        
    const double _lio_1 = 1.0f*true/V_E;
    const double _lio_2 = 135.0 * mM;
    const double _lio_3 = 1.0f*true/V_N;
    const double _lio_4 = 133.0 * mM;
    const double _lio_5 = 155.0 * mM;
    const double _lio_6 = 1.0f*0.01485884101040119/mM;
    const double _lio_7 = 1.0f*(- 0.1)/V_T;
    const double _lio_8 = 1.0f*true/V_T;


    #pragma omp parallel for schedule(static)
    for (int _i = 0; _i < (int)_num_indices; _i++)
    {
        // vector code
        const size_t _idx = _ptr_array_svmon__indices[_i];
        const size_t _vectorisation_idx = _idx;
                
        const double _____source_I_NKP_hh_1_f_NaK_hh_1_sigma_hh_1_C_Na_E_hh_1_n_Na_E = _ptr_array_hh_1_n_Na_E[_idx];
        const double ____source_I_Cl_L_hh_1_E_Cl_hh_1_C_Cl_E_hh_1_n_Cl_E = _ptr_array_hh_1_n_Cl_E[_idx];
        const double ____source_I_Cl_L_hh_1_E_Cl_hh_1_C_Cl_N_hh_1_n_Cl_N = _ptr_array_hh_1_n_Cl_N[_idx];
        const double ____source_I_KCC_hh_1_E_Cl_hh_1_C_Cl_E_hh_1_n_Cl_E = _ptr_array_hh_1_n_Cl_E[_idx];
        const double ____source_I_KCC_hh_1_E_Cl_hh_1_C_Cl_N_hh_1_n_Cl_N = _ptr_array_hh_1_n_Cl_N[_idx];
        const double ____source_I_KCC_hh_1_E_K_hh_1_C_K_E_hh_1_n_K_E = _ptr_array_hh_1_n_K_E[_idx];
        const double ____source_I_KCC_hh_1_E_K_hh_1_C_K_N_hh_1_n_K_N = _ptr_array_hh_1_n_K_N[_idx];
        const double ____source_I_K_hh_1_E_K_hh_1_C_K_E_hh_1_n_K_E = _ptr_array_hh_1_n_K_E[_idx];
        const double ____source_I_K_hh_1_E_K_hh_1_C_K_N_hh_1_n_K_N = _ptr_array_hh_1_n_K_N[_idx];
        const double ____source_I_Na_L_hh_1_E_Na_hh_1_C_Na_E_hh_1_n_Na_E = _ptr_array_hh_1_n_Na_E[_idx];
        const double ____source_I_Na_L_hh_1_E_Na_hh_1_C_Na_N_hh_1_n_Na_N = _ptr_array_hh_1_n_Na_N[_idx];
        const double ____source_I_Na_hh_1_E_Na_hh_1_C_Na_E_hh_1_n_Na_E = _ptr_array_hh_1_n_Na_E[_idx];
        const double ____source_I_Na_hh_1_E_Na_hh_1_C_Na_N_hh_1_n_Na_N = _ptr_array_hh_1_n_Na_N[_idx];
        const double ___source_E_Cl_hh_1_C_Cl_E_hh_1_n_Cl_E = _ptr_array_hh_1_n_Cl_E[_idx];
        const double ___source_E_Cl_hh_1_C_Cl_N_hh_1_n_Cl_N = _ptr_array_hh_1_n_Cl_N[_idx];
        const double ___source_E_K_hh_1_C_K_E_hh_1_n_K_E = _ptr_array_hh_1_n_K_E[_idx];
        const double ___source_E_K_hh_1_C_K_N_hh_1_n_K_N = _ptr_array_hh_1_n_K_N[_idx];
        const double ___source_E_Na_hh_1_C_Na_E_hh_1_n_Na_E = _ptr_array_hh_1_n_Na_E[_idx];
        const double ___source_E_Na_hh_1_C_Na_N_hh_1_n_Na_N = _ptr_array_hh_1_n_Na_N[_idx];
        const double ___source_I_NKP_hh_1_C_K_E_hh_1_n_K_E = _ptr_array_hh_1_n_K_E[_idx];
        const double ___source_I_NKP_hh_1_C_Na_N_hh_1_n_Na_N = _ptr_array_hh_1_n_Na_N[_idx];
        const double ___source_I_NKP_hh_1_f_NaK_hh_1_v = _ptr_array_hh_1_v[_idx];
        const double __source_C_Cl_N_hh_1_n_Cl_N = _ptr_array_hh_1_n_Cl_N[_idx];
        const double __source_C_K_N_hh_1_n_K_N = _ptr_array_hh_1_n_K_N[_idx];
        const double __source_C_Na_N_hh_1_n_Na_N = _ptr_array_hh_1_n_Na_N[_idx];
        const double __source_I_Cl_L_hh_1_v = _ptr_array_hh_1_v[_idx];
        const double __source_I_K_hh_1_n = _ptr_array_hh_1_n[_idx];
        const double __source_I_K_hh_1_v = _ptr_array_hh_1_v[_idx];
        const double __source_I_Na_L_hh_1_v = _ptr_array_hh_1_v[_idx];
        const double __source_I_Na_hh_1_h = _ptr_array_hh_1_h[_idx];
        const double __source_I_Na_hh_1_m = _ptr_array_hh_1_m[_idx];
        const double __source_I_Na_hh_1_v = _ptr_array_hh_1_v[_idx];
        const double _source_v = _ptr_array_hh_1_v[_idx];
        const double _to_record_v = _source_v;
        const double __source_E_K_hh_1_C_K_E = _clip(C0_K_E + (_lio_1 * ___source_E_K_hh_1_C_K_E_hh_1_n_K_E), false, _lio_2);
        const double __source_E_K_hh_1_C_K_N = _clip(C0_K_N + (_lio_3 * ___source_E_K_hh_1_C_K_N_hh_1_n_K_N), false, _lio_2);
        const double _source_E_K = volt * NernstPotential(__source_E_K_hh_1_C_K_E, __source_E_K_hh_1_C_K_N, true, T_exp);
        const double _to_record_E_K = _source_E_K;
        const double __source_E_Cl_hh_1_C_Cl_E = _clip(C0_Cl_E + (_lio_1 * ___source_E_Cl_hh_1_C_Cl_E_hh_1_n_Cl_E), false, _lio_4);
        const double __source_E_Cl_hh_1_C_Cl_N = _clip(C0_Cl_N + (_lio_3 * ___source_E_Cl_hh_1_C_Cl_N_hh_1_n_Cl_N), false, _lio_4);
        const double _source_E_Cl = volt * NernstPotential(__source_E_Cl_hh_1_C_Cl_E, __source_E_Cl_hh_1_C_Cl_N, - true, T_exp);
        const double _to_record_E_Cl = _source_E_Cl;
        const double __source_E_Na_hh_1_C_Na_E = _clip(C0_Na_E + (_lio_1 * ___source_E_Na_hh_1_C_Na_E_hh_1_n_Na_E), false, _lio_5);
        const double __source_E_Na_hh_1_C_Na_N = _clip(C0_Na_N + (_lio_3 * ___source_E_Na_hh_1_C_Na_N_hh_1_n_Na_N), false, _lio_5);
        const double _source_E_Na = volt * NernstPotential(__source_E_Na_hh_1_C_Na_E, __source_E_Na_hh_1_C_Na_N, true, T_exp);
        const double _to_record_E_Na = _source_E_Na;
        const double ___source_I_Na_hh_1_E_Na_hh_1_C_Na_E = _clip(C0_Na_E + (_lio_1 * ____source_I_Na_hh_1_E_Na_hh_1_C_Na_E_hh_1_n_Na_E), false, _lio_5);
        const double ___source_I_Na_hh_1_E_Na_hh_1_C_Na_N = _clip(C0_Na_N + (_lio_3 * ____source_I_Na_hh_1_E_Na_hh_1_C_Na_N_hh_1_n_Na_N), false, _lio_5);
        const double __source_I_Na_hh_1_E_Na = volt * NernstPotential(___source_I_Na_hh_1_E_Na_hh_1_C_Na_E, ___source_I_Na_hh_1_E_Na_hh_1_C_Na_N, true, T_exp);
        const double _source_I_Na = g_Na * (((_brian_pow(__source_I_Na_hh_1_m, 3)) * __source_I_Na_hh_1_h) * (__source_I_Na_hh_1_v - __source_I_Na_hh_1_E_Na));
        const double _to_record_I_Na = _source_I_Na;
        const double ___source_I_Na_L_hh_1_E_Na_hh_1_C_Na_E = _clip(C0_Na_E + (_lio_1 * ____source_I_Na_L_hh_1_E_Na_hh_1_C_Na_E_hh_1_n_Na_E), false, _lio_5);
        const double ___source_I_Na_L_hh_1_E_Na_hh_1_C_Na_N = _clip(C0_Na_N + (_lio_3 * ____source_I_Na_L_hh_1_E_Na_hh_1_C_Na_N_hh_1_n_Na_N), false, _lio_5);
        const double __source_I_Na_L_hh_1_E_Na = volt * NernstPotential(___source_I_Na_L_hh_1_E_Na_hh_1_C_Na_E, ___source_I_Na_L_hh_1_E_Na_hh_1_C_Na_N, true, T_exp);
        const double _source_I_Na_L = g_Na_L * (__source_I_Na_L_hh_1_v - __source_I_Na_L_hh_1_E_Na);
        const double _to_record_I_Na_L = _source_I_Na_L;
        const double ___source_I_K_hh_1_E_K_hh_1_C_K_E = _clip(C0_K_E + (_lio_1 * ____source_I_K_hh_1_E_K_hh_1_C_K_E_hh_1_n_K_E), false, _lio_2);
        const double ___source_I_K_hh_1_E_K_hh_1_C_K_N = _clip(C0_K_N + (_lio_3 * ____source_I_K_hh_1_E_K_hh_1_C_K_N_hh_1_n_K_N), false, _lio_2);
        const double __source_I_K_hh_1_E_K = volt * NernstPotential(___source_I_K_hh_1_E_K_hh_1_C_K_E, ___source_I_K_hh_1_E_K_hh_1_C_K_N, true, T_exp);
        const double _source_I_K = g_K * (__source_I_K_hh_1_n * (__source_I_K_hh_1_v - __source_I_K_hh_1_E_K));
        const double _to_record_I_K = _source_I_K;
        const double ___source_I_Cl_L_hh_1_E_Cl_hh_1_C_Cl_E = _clip(C0_Cl_E + (_lio_1 * ____source_I_Cl_L_hh_1_E_Cl_hh_1_C_Cl_E_hh_1_n_Cl_E), false, _lio_4);
        const double ___source_I_Cl_L_hh_1_E_Cl_hh_1_C_Cl_N = _clip(C0_Cl_N + (_lio_3 * ____source_I_Cl_L_hh_1_E_Cl_hh_1_C_Cl_N_hh_1_n_Cl_N), false, _lio_4);
        const double __source_I_Cl_L_hh_1_E_Cl = volt * NernstPotential(___source_I_Cl_L_hh_1_E_Cl_hh_1_C_Cl_E, ___source_I_Cl_L_hh_1_E_Cl_hh_1_C_Cl_N, - true, T_exp);
        const double _source_I_Cl_L = g_Cl_L * (__source_I_Cl_L_hh_1_v - __source_I_Cl_L_hh_1_E_Cl);
        const double _to_record_I_Cl_L = _source_I_Cl_L;
        const double ___source_I_KCC_hh_1_E_Cl_hh_1_C_Cl_E = _clip(C0_Cl_E + (_lio_1 * ____source_I_KCC_hh_1_E_Cl_hh_1_C_Cl_E_hh_1_n_Cl_E), false, _lio_4);
        const double ___source_I_KCC_hh_1_E_Cl_hh_1_C_Cl_N = _clip(C0_Cl_N + (_lio_3 * ____source_I_KCC_hh_1_E_Cl_hh_1_C_Cl_N_hh_1_n_Cl_N), false, _lio_4);
        const double ___source_I_KCC_hh_1_E_K_hh_1_C_K_E = _clip(C0_K_E + (_lio_1 * ____source_I_KCC_hh_1_E_K_hh_1_C_K_E_hh_1_n_K_E), false, _lio_2);
        const double ___source_I_KCC_hh_1_E_K_hh_1_C_K_N = _clip(C0_K_N + (_lio_3 * ____source_I_KCC_hh_1_E_K_hh_1_C_K_N_hh_1_n_K_N), false, _lio_2);
        const double __source_I_KCC_hh_1_E_Cl = volt * NernstPotential(___source_I_KCC_hh_1_E_Cl_hh_1_C_Cl_E, ___source_I_KCC_hh_1_E_Cl_hh_1_C_Cl_N, - true, T_exp);
        const double __source_I_KCC_hh_1_E_K = volt * NernstPotential(___source_I_KCC_hh_1_E_K_hh_1_C_K_E, ___source_I_KCC_hh_1_E_K_hh_1_C_K_N, true, T_exp);
        const double _source_I_KCC = g_KCC * (__source_I_KCC_hh_1_E_K - __source_I_KCC_hh_1_E_Cl);
        const double _to_record_I_KCC = _source_I_KCC;
        const double ____source_I_NKP_hh_1_f_NaK_hh_1_sigma_hh_1_C_Na_E = _clip(C0_Na_E + (_lio_1 * _____source_I_NKP_hh_1_f_NaK_hh_1_sigma_hh_1_C_Na_E_hh_1_n_Na_E), false, _lio_5);
        const double __source_I_NKP_hh_1_C_K_E = _clip(C0_K_E + (_lio_1 * ___source_I_NKP_hh_1_C_K_E_hh_1_n_K_E), false, _lio_2);
        const double __source_I_NKP_hh_1_C_Na_N = _clip(C0_Na_N + (_lio_3 * ___source_I_NKP_hh_1_C_Na_N_hh_1_n_Na_N), false, _lio_5);
        const double ___source_I_NKP_hh_1_f_NaK_hh_1_sigma = 0.14285714285714285 * ((- true) + expm1(_lio_6 * ____source_I_NKP_hh_1_f_NaK_hh_1_sigma_hh_1_C_Na_E));
        const double __source_I_NKP_hh_1_f_NaK = 1.0f*true/(true + ((0.1245 * exp(_lio_7 * ___source_I_NKP_hh_1_f_NaK_hh_1_v)) + (0.0365 * (___source_I_NKP_hh_1_f_NaK_hh_1_sigma * expm1(_lio_8 * (- ___source_I_NKP_hh_1_f_NaK_hh_1_v))))));
        const double _source_I_NKP = 1.0f*(I_NKP_max * (__source_I_NKP_hh_1_f_NaK * Hill(__source_I_NKP_hh_1_C_K_E, zeta_K, true)))/(true + (_brian_pow(1.0f*zeta_Na/__source_I_NKP_hh_1_C_Na_N, 1.5)));
        const double _to_record_I_NKP = _source_I_NKP;
        const double _source_C_Cl_N = _clip(C0_Cl_N + (_lio_3 * __source_C_Cl_N_hh_1_n_Cl_N), false, _lio_4);
        const double _to_record_C_Cl_N = _source_C_Cl_N;
        const double _source_C_Na_N = _clip(C0_Na_N + (_lio_3 * __source_C_Na_N_hh_1_n_Na_N), false, _lio_5);
        const double _to_record_C_Na_N = _source_C_Na_N;
        const double _source_C_K_N = _clip(C0_K_N + (_lio_3 * __source_C_K_N_hh_1_n_K_N), false, _lio_2);
        const double _to_record_C_K_N = _source_C_K_N;


        _dynamic_array_svmon_C_Cl_N(_new_size-1, _i) = _to_record_C_Cl_N;
        _dynamic_array_svmon_C_K_N(_new_size-1, _i) = _to_record_C_K_N;
        _dynamic_array_svmon_C_Na_N(_new_size-1, _i) = _to_record_C_Na_N;
        _dynamic_array_svmon_E_Cl(_new_size-1, _i) = _to_record_E_Cl;
        _dynamic_array_svmon_E_K(_new_size-1, _i) = _to_record_E_K;
        _dynamic_array_svmon_E_Na(_new_size-1, _i) = _to_record_E_Na;
        _dynamic_array_svmon_I_Cl_L(_new_size-1, _i) = _to_record_I_Cl_L;
        _dynamic_array_svmon_I_K(_new_size-1, _i) = _to_record_I_K;
        _dynamic_array_svmon_I_KCC(_new_size-1, _i) = _to_record_I_KCC;
        _dynamic_array_svmon_I_Na(_new_size-1, _i) = _to_record_I_Na;
        _dynamic_array_svmon_I_Na_L(_new_size-1, _i) = _to_record_I_Na_L;
        _dynamic_array_svmon_I_NKP(_new_size-1, _i) = _to_record_I_NKP;
        _dynamic_array_svmon_v(_new_size-1, _i) = _to_record_v;
    }

    _ptr_array_svmon_N[0] = _new_size;


}


