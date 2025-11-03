#include "code_objects/hh_stateupdater_codeobject.h"
#include "objects.h"
#include "brianlib/common_math.h"
#include "brianlib/stdint_compat.h"
#include<cmath>
#include<ctime>
#include<iostream>
#include<fstream>
#include<climits>
#include "gsl_const_mksa.h"

////// SUPPORT CODE ///////
namespace {
        
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



void _run_hh_stateupdater_codeobject()
{
    using namespace brian;


    ///// CONSTANTS ///////////
    const double C0_Cl_E = 130.0;
const double C0_Cl_N = 5.0;
const double C0_K_E = 3.0;
const double C0_K_N = 130.0;
const double C0_Na_E = 145.0;
const double C0_Na_N = 10.0;
const double F = 96485.33212;
const double I_dc = 0.1584893192461114;
const double I_ramp = false;
const int64_t N = true;
const double S_N = 1.5e-10;
const double T_adj = 3.7910938669360985;
const int64_t T_exp = 37;
const double T_ramp = true;
const double U_h = - 0.066;
const double U_m = - 0.038;
const double U_n = 0.0187;
const double V_E = 5e-16;
const double V_N = 5e-16;
const double W_h = 0.006;
const double W_m = 0.006;
const double W_n = 0.0097;
const double c_m = 0.009999999999999998;
const size_t _numdt = 1;
const double g_Cl = 3.3799999999999994;
const double g_K = 6380.0;
const double g_KCC = 1999.9999999999998;
const double g_Na = 20400.0;
const size_t _numh = 1;
const size_t _numm = 1;
const double mM = true;
const double ms = 0.001;
const double mvolt = 0.001;
const size_t _numn = 1;
const size_t _numn_Cl_E = 1;
const size_t _numn_Cl_N = 1;
const size_t _numn_K_E = 1;
const size_t _numn_K_N = 1;
const size_t _numn_Na_E = 1;
const size_t _numn_Na_N = 1;
const double second = true;
const size_t _numt = 1;
const size_t _numv = 1;
const double volt = true;
    ///// POINTERS ////////////
        
    double*   _ptr_array_hh_clock_dt = _array_hh_clock_dt;
    double*   _ptr_array_hh_h = _array_hh_h;
    double*   _ptr_array_hh_m = _array_hh_m;
    double*   _ptr_array_hh_n = _array_hh_n;
    double*   _ptr_array_hh_n_Cl_E = _array_hh_n_Cl_E;
    double*   _ptr_array_hh_n_Cl_N = _array_hh_n_Cl_N;
    double*   _ptr_array_hh_n_K_E = _array_hh_n_K_E;
    double*   _ptr_array_hh_n_K_N = _array_hh_n_K_N;
    double*   _ptr_array_hh_n_Na_E = _array_hh_n_Na_E;
    double*   _ptr_array_hh_n_Na_N = _array_hh_n_Na_N;
    double*   _ptr_array_hh_clock_t = _array_hh_clock_t;
    double*   _ptr_array_hh_v = _array_hh_v;


    //// MAIN CODE ////////////
    // scalar code
    const size_t _vectorisation_idx = -1;
        
    const double dt = _ptr_array_hh_clock_dt[0];
    const double t = _ptr_array_hh_clock_t[0];
    const double _lio_1 = 1000.0 * (T_adj * dt);
    const double _lio_2 = 1.0f*true/(mvolt * second);
    const double _lio_3 = (- 0.015) * U_h;
    const double _lio_4 = 1.0f*true/W_h;
    const double _lio_5 = - U_h;
    const double _lio_6 = 0.015 * U_h;
    const double _lio_7 = 0.182 * U_m;
    const double _lio_8 = (- 0.124) * U_m;
    const double _lio_9 = 1.0f*true/W_m;
    const double _lio_10 = - U_m;
    const double _lio_11 = 1.0f*(0.25 * (T_adj * dt))/ms;
    const double _lio_12 = 1.0f*true/W_n;
    const double _lio_13 = 1.0f*true/mvolt;
    const double _lio_14 = (- 1.28137743543271) * mvolt;
    const double _lio_15 = 1.0f*((- S_N) * dt)/F;
    const double _lio_16 = 1.0f*true/V_E;
    const double _lio_17 = 133.0 * mM;
    const double _lio_18 = 1.0f*true/V_N;
    const double _lio_19 = - volt;
    const double _lio_20 = 135.0 * mM;
    const double _lio_21 = 1.0f*(S_N * dt)/F;
    const double _lio_22 = 1.0f*((S_N * dt) * g_Na)/F;
    const double _lio_23 = 155.0 * mM;
    const double _lio_24 = 1.0f*(((- S_N) * dt) * g_Na)/F;
    const double _lio_25 = 1.0f*dt/c_m;
    const double _lio_26 = I_dc + (1.0f*(I_ramp * t)/T_ramp);
    const double _lio_27 = false - (1.28137743543271 * mvolt);
    const double _lio_28 = I_dc + (1.0f*(I_ramp * ((0.5 * dt) + t))/T_ramp);
    const double _lio_29 = I_dc + (1.0f*(I_ramp * (dt + t))/T_ramp);


    const int _N = N;
    #pragma omp parallel for schedule(static)
    for(int _idx=0; _idx<_N; _idx++)
    {
        // vector code
        const size_t _vectorisation_idx = _idx;
                
        double h = _ptr_array_hh_h[_idx];
        double m = _ptr_array_hh_m[_idx];
        double n = _ptr_array_hh_n[_idx];
        double n_Cl_E = _ptr_array_hh_n_Cl_E[_idx];
        double n_Cl_N = _ptr_array_hh_n_Cl_N[_idx];
        double n_K_E = _ptr_array_hh_n_K_E[_idx];
        double n_K_N = _ptr_array_hh_n_K_N[_idx];
        double n_Na_E = _ptr_array_hh_n_Na_E[_idx];
        double n_Na_N = _ptr_array_hh_n_Na_N[_idx];
        double v = _ptr_array_hh_v[_idx];
        const double __k_1_h = _lio_1 * (((- h) + (1.0f*(_lio_2 * (_lio_3 + (0.015 * v)))/(((1.0f*(_lio_2 * (_lio_3 + (0.015 * v)))/expm1(_lio_4 * (_lio_5 + v))) + (1.0f*(_lio_2 * (_lio_6 - (0.015 * v)))/expm1(_lio_4 * (U_h - v)))) * expm1(_lio_4 * (_lio_5 + v))))) * ((1.0f*(_lio_2 * (_lio_3 + (0.015 * v)))/expm1(_lio_4 * (_lio_5 + v))) + (1.0f*(_lio_2 * (_lio_6 - (0.015 * v)))/expm1(_lio_4 * (U_h - v)))));
        const double __k_1_m = _lio_1 * (((- m) + (1.0f*(_lio_2 * (_lio_7 - (0.182 * v)))/(((1.0f*(_lio_2 * (_lio_8 + (0.124 * v)))/expm1(_lio_9 * (_lio_10 + v))) + (1.0f*(_lio_2 * (_lio_7 - (0.182 * v)))/expm1(_lio_9 * (U_m - v)))) * expm1(_lio_9 * (U_m - v))))) * ((1.0f*(_lio_2 * (_lio_8 + (0.124 * v)))/expm1(_lio_9 * (_lio_10 + v))) + (1.0f*(_lio_2 * (_lio_7 - (0.182 * v)))/expm1(_lio_9 * (U_m - v)))));
        const double __k_1_n = _lio_11 * (((- n) + (1.0f*true/(true + exp(_lio_12 * (U_n - v))))) * (true + exp(_lio_13 * (_lio_14 - (0.0226551880380607 * v)))));
        const double __k_1_n_Cl_E = _lio_15 * ((g_Cl * (v - (volt * NernstPotential(_clip(C0_Cl_E + (_lio_16 * n_Cl_E), false, _lio_17), _clip(C0_Cl_N + (_lio_18 * n_Cl_N), false, _lio_17), true, T_exp)))) + (g_KCC * ((_lio_19 * NernstPotential(_clip(C0_Cl_E + (_lio_16 * n_Cl_E), false, _lio_17), _clip(C0_Cl_N + (_lio_18 * n_Cl_N), false, _lio_17), true, T_exp)) + (volt * NernstPotential(_clip(C0_K_E + (_lio_16 * n_K_E), false, _lio_20), _clip(C0_K_N + (_lio_18 * n_K_N), false, _lio_20), true, T_exp)))));
        const double __k_1_n_Cl_N = _lio_21 * ((g_Cl * (v - (volt * NernstPotential(_clip(C0_Cl_E + (_lio_16 * n_Cl_E), false, _lio_17), _clip(C0_Cl_N + (_lio_18 * n_Cl_N), false, _lio_17), true, T_exp)))) + (g_KCC * ((_lio_19 * NernstPotential(_clip(C0_Cl_E + (_lio_16 * n_Cl_E), false, _lio_17), _clip(C0_Cl_N + (_lio_18 * n_Cl_N), false, _lio_17), true, T_exp)) + (volt * NernstPotential(_clip(C0_K_E + (_lio_16 * n_K_E), false, _lio_20), _clip(C0_K_N + (_lio_18 * n_K_N), false, _lio_20), true, T_exp)))));
        const double __k_1_n_K_E = _lio_21 * ((g_K * (n * (v - (volt * NernstPotential(_clip(C0_K_E + (_lio_16 * n_K_E), false, _lio_20), _clip(C0_K_N + (_lio_18 * n_K_N), false, _lio_20), true, T_exp))))) - (g_KCC * ((_lio_19 * NernstPotential(_clip(C0_Cl_E + (_lio_16 * n_Cl_E), false, _lio_17), _clip(C0_Cl_N + (_lio_18 * n_Cl_N), false, _lio_17), true, T_exp)) + (volt * NernstPotential(_clip(C0_K_E + (_lio_16 * n_K_E), false, _lio_20), _clip(C0_K_N + (_lio_18 * n_K_N), false, _lio_20), true, T_exp)))));
        const double __k_1_n_K_N = _lio_15 * ((g_K * (n * (v - (volt * NernstPotential(_clip(C0_K_E + (_lio_16 * n_K_E), false, _lio_20), _clip(C0_K_N + (_lio_18 * n_K_N), false, _lio_20), true, T_exp))))) - (g_KCC * ((_lio_19 * NernstPotential(_clip(C0_Cl_E + (_lio_16 * n_Cl_E), false, _lio_17), _clip(C0_Cl_N + (_lio_18 * n_Cl_N), false, _lio_17), true, T_exp)) + (volt * NernstPotential(_clip(C0_K_E + (_lio_16 * n_K_E), false, _lio_20), _clip(C0_K_N + (_lio_18 * n_K_N), false, _lio_20), true, T_exp)))));
        const double __k_1_n_Na_E = _lio_22 * ((h * (_brian_pow(m, 3))) * (v - (volt * NernstPotential(_clip(C0_Na_E + (_lio_16 * n_Na_E), false, _lio_23), _clip(C0_Na_N + (_lio_18 * n_Na_N), false, _lio_23), true, T_exp))));
        const double __k_1_n_Na_N = _lio_24 * ((h * (_brian_pow(m, 3))) * (v - (volt * NernstPotential(_clip(C0_Na_E + (_lio_16 * n_Na_E), false, _lio_23), _clip(C0_Na_N + (_lio_18 * n_Na_N), false, _lio_23), true, T_exp))));
        const double __k_1_v = _lio_25 * (_lio_26 - (((g_Cl * (v - (volt * NernstPotential(_clip(C0_Cl_E + (_lio_16 * n_Cl_E), false, _lio_17), _clip(C0_Cl_N + (_lio_18 * n_Cl_N), false, _lio_17), true, T_exp)))) + (g_K * (n * (v - (volt * NernstPotential(_clip(C0_K_E + (_lio_16 * n_K_E), false, _lio_20), _clip(C0_K_N + (_lio_18 * n_K_N), false, _lio_20), true, T_exp)))))) + (g_Na * ((h * (_brian_pow(m, 3))) * (v - (volt * NernstPotential(_clip(C0_Na_E + (_lio_16 * n_Na_E), false, _lio_23), _clip(C0_Na_N + (_lio_18 * n_Na_N), false, _lio_23), true, T_exp)))))));
        const double __k_2_h = _lio_1 * (((1.0f*(_lio_2 * (_lio_3 + ((0.0075 * __k_1_v) + (0.015 * v))))/expm1(_lio_4 * (_lio_5 + ((0.5 * __k_1_v) + v)))) + (1.0f*(_lio_2 * (_lio_6 - ((0.0075 * __k_1_v) + (0.015 * v))))/expm1(_lio_4 * (U_h - ((0.5 * __k_1_v) + v))))) * (((0.5 * (- __k_1_h)) + (1.0f*(_lio_2 * (_lio_3 + ((0.0075 * __k_1_v) + (0.015 * v))))/(((1.0f*(_lio_2 * (_lio_3 + ((0.0075 * __k_1_v) + (0.015 * v))))/expm1(_lio_4 * (_lio_5 + ((0.5 * __k_1_v) + v)))) + (1.0f*(_lio_2 * (_lio_6 - ((0.0075 * __k_1_v) + (0.015 * v))))/expm1(_lio_4 * (U_h - ((0.5 * __k_1_v) + v))))) * expm1(_lio_4 * (_lio_5 + ((0.5 * __k_1_v) + v)))))) - h));
        const double __k_2_m = _lio_1 * (((1.0f*(_lio_2 * (_lio_8 + ((0.062 * __k_1_v) + (0.124 * v))))/expm1(_lio_9 * (_lio_10 + ((0.5 * __k_1_v) + v)))) + (1.0f*(_lio_2 * (_lio_7 - ((0.091 * __k_1_v) + (0.182 * v))))/expm1(_lio_9 * (U_m - ((0.5 * __k_1_v) + v))))) * (((0.5 * (- __k_1_m)) + (1.0f*(_lio_2 * (_lio_7 - ((0.091 * __k_1_v) + (0.182 * v))))/(((1.0f*(_lio_2 * (_lio_8 + ((0.062 * __k_1_v) + (0.124 * v))))/expm1(_lio_9 * (_lio_10 + ((0.5 * __k_1_v) + v)))) + (1.0f*(_lio_2 * (_lio_7 - ((0.091 * __k_1_v) + (0.182 * v))))/expm1(_lio_9 * (U_m - ((0.5 * __k_1_v) + v))))) * expm1(_lio_9 * (U_m - ((0.5 * __k_1_v) + v)))))) - m));
        const double __k_2_n = _lio_11 * ((true + exp(_lio_13 * ((_lio_27 + ((- 0.0113275940190304) * __k_1_v)) - (0.0226551880380607 * v)))) * (((0.5 * (- __k_1_n)) + (1.0f*true/(true + exp(_lio_12 * (U_n - ((0.5 * __k_1_v) + v)))))) - n));
        const double __k_2_n_Cl_E = _lio_15 * ((g_Cl * (((0.5 * __k_1_v) + v) - (volt * NernstPotential(_clip(C0_Cl_E + (_lio_16 * ((0.5 * __k_1_n_Cl_E) + n_Cl_E)), false, _lio_17), _clip(C0_Cl_N + (_lio_18 * ((0.5 * __k_1_n_Cl_N) + n_Cl_N)), false, _lio_17), true, T_exp)))) + (g_KCC * ((_lio_19 * NernstPotential(_clip(C0_Cl_E + (_lio_16 * ((0.5 * __k_1_n_Cl_E) + n_Cl_E)), false, _lio_17), _clip(C0_Cl_N + (_lio_18 * ((0.5 * __k_1_n_Cl_N) + n_Cl_N)), false, _lio_17), true, T_exp)) + (volt * NernstPotential(_clip(C0_K_E + (_lio_16 * ((0.5 * __k_1_n_K_E) + n_K_E)), false, _lio_20), _clip(C0_K_N + (_lio_18 * ((0.5 * __k_1_n_K_N) + n_K_N)), false, _lio_20), true, T_exp)))));
        const double __k_2_n_Cl_N = _lio_21 * ((g_Cl * (((0.5 * __k_1_v) + v) - (volt * NernstPotential(_clip(C0_Cl_E + (_lio_16 * ((0.5 * __k_1_n_Cl_E) + n_Cl_E)), false, _lio_17), _clip(C0_Cl_N + (_lio_18 * ((0.5 * __k_1_n_Cl_N) + n_Cl_N)), false, _lio_17), true, T_exp)))) + (g_KCC * ((_lio_19 * NernstPotential(_clip(C0_Cl_E + (_lio_16 * ((0.5 * __k_1_n_Cl_E) + n_Cl_E)), false, _lio_17), _clip(C0_Cl_N + (_lio_18 * ((0.5 * __k_1_n_Cl_N) + n_Cl_N)), false, _lio_17), true, T_exp)) + (volt * NernstPotential(_clip(C0_K_E + (_lio_16 * ((0.5 * __k_1_n_K_E) + n_K_E)), false, _lio_20), _clip(C0_K_N + (_lio_18 * ((0.5 * __k_1_n_K_N) + n_K_N)), false, _lio_20), true, T_exp)))));
        const double __k_2_n_K_E = _lio_21 * ((g_K * (((0.5 * __k_1_n) + n) * (((0.5 * __k_1_v) + v) - (volt * NernstPotential(_clip(C0_K_E + (_lio_16 * ((0.5 * __k_1_n_K_E) + n_K_E)), false, _lio_20), _clip(C0_K_N + (_lio_18 * ((0.5 * __k_1_n_K_N) + n_K_N)), false, _lio_20), true, T_exp))))) - (g_KCC * ((_lio_19 * NernstPotential(_clip(C0_Cl_E + (_lio_16 * ((0.5 * __k_1_n_Cl_E) + n_Cl_E)), false, _lio_17), _clip(C0_Cl_N + (_lio_18 * ((0.5 * __k_1_n_Cl_N) + n_Cl_N)), false, _lio_17), true, T_exp)) + (volt * NernstPotential(_clip(C0_K_E + (_lio_16 * ((0.5 * __k_1_n_K_E) + n_K_E)), false, _lio_20), _clip(C0_K_N + (_lio_18 * ((0.5 * __k_1_n_K_N) + n_K_N)), false, _lio_20), true, T_exp)))));
        const double __k_2_n_K_N = _lio_15 * ((g_K * (((0.5 * __k_1_n) + n) * (((0.5 * __k_1_v) + v) - (volt * NernstPotential(_clip(C0_K_E + (_lio_16 * ((0.5 * __k_1_n_K_E) + n_K_E)), false, _lio_20), _clip(C0_K_N + (_lio_18 * ((0.5 * __k_1_n_K_N) + n_K_N)), false, _lio_20), true, T_exp))))) - (g_KCC * ((_lio_19 * NernstPotential(_clip(C0_Cl_E + (_lio_16 * ((0.5 * __k_1_n_Cl_E) + n_Cl_E)), false, _lio_17), _clip(C0_Cl_N + (_lio_18 * ((0.5 * __k_1_n_Cl_N) + n_Cl_N)), false, _lio_17), true, T_exp)) + (volt * NernstPotential(_clip(C0_K_E + (_lio_16 * ((0.5 * __k_1_n_K_E) + n_K_E)), false, _lio_20), _clip(C0_K_N + (_lio_18 * ((0.5 * __k_1_n_K_N) + n_K_N)), false, _lio_20), true, T_exp)))));
        const double __k_2_n_Na_E = _lio_22 * ((((0.5 * __k_1_h) + h) * (_brian_pow((0.5 * __k_1_m) + m, 3))) * (((0.5 * __k_1_v) + v) - (volt * NernstPotential(_clip(C0_Na_E + (_lio_16 * ((0.5 * __k_1_n_Na_E) + n_Na_E)), false, _lio_23), _clip(C0_Na_N + (_lio_18 * ((0.5 * __k_1_n_Na_N) + n_Na_N)), false, _lio_23), true, T_exp))));
        const double __k_2_n_Na_N = _lio_24 * ((((0.5 * __k_1_h) + h) * (_brian_pow((0.5 * __k_1_m) + m, 3))) * (((0.5 * __k_1_v) + v) - (volt * NernstPotential(_clip(C0_Na_E + (_lio_16 * ((0.5 * __k_1_n_Na_E) + n_Na_E)), false, _lio_23), _clip(C0_Na_N + (_lio_18 * ((0.5 * __k_1_n_Na_N) + n_Na_N)), false, _lio_23), true, T_exp))));
        const double __k_2_v = _lio_25 * (_lio_28 - (((g_Cl * (((0.5 * __k_1_v) + v) - (volt * NernstPotential(_clip(C0_Cl_E + (_lio_16 * ((0.5 * __k_1_n_Cl_E) + n_Cl_E)), false, _lio_17), _clip(C0_Cl_N + (_lio_18 * ((0.5 * __k_1_n_Cl_N) + n_Cl_N)), false, _lio_17), true, T_exp)))) + (g_K * (((0.5 * __k_1_n) + n) * (((0.5 * __k_1_v) + v) - (volt * NernstPotential(_clip(C0_K_E + (_lio_16 * ((0.5 * __k_1_n_K_E) + n_K_E)), false, _lio_20), _clip(C0_K_N + (_lio_18 * ((0.5 * __k_1_n_K_N) + n_K_N)), false, _lio_20), true, T_exp)))))) + (g_Na * ((((0.5 * __k_1_h) + h) * (_brian_pow((0.5 * __k_1_m) + m, 3))) * (((0.5 * __k_1_v) + v) - (volt * NernstPotential(_clip(C0_Na_E + (_lio_16 * ((0.5 * __k_1_n_Na_E) + n_Na_E)), false, _lio_23), _clip(C0_Na_N + (_lio_18 * ((0.5 * __k_1_n_Na_N) + n_Na_N)), false, _lio_23), true, T_exp)))))));
        const double __k_3_h = _lio_1 * (((1.0f*(_lio_2 * (_lio_3 + ((0.0075 * __k_2_v) + (0.015 * v))))/expm1(_lio_4 * (_lio_5 + ((0.5 * __k_2_v) + v)))) + (1.0f*(_lio_2 * (_lio_6 - ((0.0075 * __k_2_v) + (0.015 * v))))/expm1(_lio_4 * (U_h - ((0.5 * __k_2_v) + v))))) * (((0.5 * (- __k_2_h)) + (1.0f*(_lio_2 * (_lio_3 + ((0.0075 * __k_2_v) + (0.015 * v))))/(((1.0f*(_lio_2 * (_lio_3 + ((0.0075 * __k_2_v) + (0.015 * v))))/expm1(_lio_4 * (_lio_5 + ((0.5 * __k_2_v) + v)))) + (1.0f*(_lio_2 * (_lio_6 - ((0.0075 * __k_2_v) + (0.015 * v))))/expm1(_lio_4 * (U_h - ((0.5 * __k_2_v) + v))))) * expm1(_lio_4 * (_lio_5 + ((0.5 * __k_2_v) + v)))))) - h));
        const double __k_3_m = _lio_1 * (((1.0f*(_lio_2 * (_lio_8 + ((0.062 * __k_2_v) + (0.124 * v))))/expm1(_lio_9 * (_lio_10 + ((0.5 * __k_2_v) + v)))) + (1.0f*(_lio_2 * (_lio_7 - ((0.091 * __k_2_v) + (0.182 * v))))/expm1(_lio_9 * (U_m - ((0.5 * __k_2_v) + v))))) * (((0.5 * (- __k_2_m)) + (1.0f*(_lio_2 * (_lio_7 - ((0.091 * __k_2_v) + (0.182 * v))))/(((1.0f*(_lio_2 * (_lio_8 + ((0.062 * __k_2_v) + (0.124 * v))))/expm1(_lio_9 * (_lio_10 + ((0.5 * __k_2_v) + v)))) + (1.0f*(_lio_2 * (_lio_7 - ((0.091 * __k_2_v) + (0.182 * v))))/expm1(_lio_9 * (U_m - ((0.5 * __k_2_v) + v))))) * expm1(_lio_9 * (U_m - ((0.5 * __k_2_v) + v)))))) - m));
        const double __k_3_n = _lio_11 * ((true + exp(_lio_13 * ((_lio_27 + ((- 0.0113275940190304) * __k_2_v)) - (0.0226551880380607 * v)))) * (((0.5 * (- __k_2_n)) + (1.0f*true/(true + exp(_lio_12 * (U_n - ((0.5 * __k_2_v) + v)))))) - n));
        const double __k_3_n_Cl_E = _lio_15 * ((g_Cl * (((0.5 * __k_2_v) + v) - (volt * NernstPotential(_clip(C0_Cl_E + (_lio_16 * ((0.5 * __k_2_n_Cl_E) + n_Cl_E)), false, _lio_17), _clip(C0_Cl_N + (_lio_18 * ((0.5 * __k_2_n_Cl_N) + n_Cl_N)), false, _lio_17), true, T_exp)))) + (g_KCC * ((_lio_19 * NernstPotential(_clip(C0_Cl_E + (_lio_16 * ((0.5 * __k_2_n_Cl_E) + n_Cl_E)), false, _lio_17), _clip(C0_Cl_N + (_lio_18 * ((0.5 * __k_2_n_Cl_N) + n_Cl_N)), false, _lio_17), true, T_exp)) + (volt * NernstPotential(_clip(C0_K_E + (_lio_16 * ((0.5 * __k_2_n_K_E) + n_K_E)), false, _lio_20), _clip(C0_K_N + (_lio_18 * ((0.5 * __k_2_n_K_N) + n_K_N)), false, _lio_20), true, T_exp)))));
        const double __k_3_n_Cl_N = _lio_21 * ((g_Cl * (((0.5 * __k_2_v) + v) - (volt * NernstPotential(_clip(C0_Cl_E + (_lio_16 * ((0.5 * __k_2_n_Cl_E) + n_Cl_E)), false, _lio_17), _clip(C0_Cl_N + (_lio_18 * ((0.5 * __k_2_n_Cl_N) + n_Cl_N)), false, _lio_17), true, T_exp)))) + (g_KCC * ((_lio_19 * NernstPotential(_clip(C0_Cl_E + (_lio_16 * ((0.5 * __k_2_n_Cl_E) + n_Cl_E)), false, _lio_17), _clip(C0_Cl_N + (_lio_18 * ((0.5 * __k_2_n_Cl_N) + n_Cl_N)), false, _lio_17), true, T_exp)) + (volt * NernstPotential(_clip(C0_K_E + (_lio_16 * ((0.5 * __k_2_n_K_E) + n_K_E)), false, _lio_20), _clip(C0_K_N + (_lio_18 * ((0.5 * __k_2_n_K_N) + n_K_N)), false, _lio_20), true, T_exp)))));
        const double __k_3_n_K_E = _lio_21 * ((g_K * (((0.5 * __k_2_n) + n) * (((0.5 * __k_2_v) + v) - (volt * NernstPotential(_clip(C0_K_E + (_lio_16 * ((0.5 * __k_2_n_K_E) + n_K_E)), false, _lio_20), _clip(C0_K_N + (_lio_18 * ((0.5 * __k_2_n_K_N) + n_K_N)), false, _lio_20), true, T_exp))))) - (g_KCC * ((_lio_19 * NernstPotential(_clip(C0_Cl_E + (_lio_16 * ((0.5 * __k_2_n_Cl_E) + n_Cl_E)), false, _lio_17), _clip(C0_Cl_N + (_lio_18 * ((0.5 * __k_2_n_Cl_N) + n_Cl_N)), false, _lio_17), true, T_exp)) + (volt * NernstPotential(_clip(C0_K_E + (_lio_16 * ((0.5 * __k_2_n_K_E) + n_K_E)), false, _lio_20), _clip(C0_K_N + (_lio_18 * ((0.5 * __k_2_n_K_N) + n_K_N)), false, _lio_20), true, T_exp)))));
        const double __k_3_n_K_N = _lio_15 * ((g_K * (((0.5 * __k_2_n) + n) * (((0.5 * __k_2_v) + v) - (volt * NernstPotential(_clip(C0_K_E + (_lio_16 * ((0.5 * __k_2_n_K_E) + n_K_E)), false, _lio_20), _clip(C0_K_N + (_lio_18 * ((0.5 * __k_2_n_K_N) + n_K_N)), false, _lio_20), true, T_exp))))) - (g_KCC * ((_lio_19 * NernstPotential(_clip(C0_Cl_E + (_lio_16 * ((0.5 * __k_2_n_Cl_E) + n_Cl_E)), false, _lio_17), _clip(C0_Cl_N + (_lio_18 * ((0.5 * __k_2_n_Cl_N) + n_Cl_N)), false, _lio_17), true, T_exp)) + (volt * NernstPotential(_clip(C0_K_E + (_lio_16 * ((0.5 * __k_2_n_K_E) + n_K_E)), false, _lio_20), _clip(C0_K_N + (_lio_18 * ((0.5 * __k_2_n_K_N) + n_K_N)), false, _lio_20), true, T_exp)))));
        const double __k_3_n_Na_E = _lio_22 * ((((0.5 * __k_2_h) + h) * (_brian_pow((0.5 * __k_2_m) + m, 3))) * (((0.5 * __k_2_v) + v) - (volt * NernstPotential(_clip(C0_Na_E + (_lio_16 * ((0.5 * __k_2_n_Na_E) + n_Na_E)), false, _lio_23), _clip(C0_Na_N + (_lio_18 * ((0.5 * __k_2_n_Na_N) + n_Na_N)), false, _lio_23), true, T_exp))));
        const double __k_3_n_Na_N = _lio_24 * ((((0.5 * __k_2_h) + h) * (_brian_pow((0.5 * __k_2_m) + m, 3))) * (((0.5 * __k_2_v) + v) - (volt * NernstPotential(_clip(C0_Na_E + (_lio_16 * ((0.5 * __k_2_n_Na_E) + n_Na_E)), false, _lio_23), _clip(C0_Na_N + (_lio_18 * ((0.5 * __k_2_n_Na_N) + n_Na_N)), false, _lio_23), true, T_exp))));
        const double __k_3_v = _lio_25 * (_lio_28 - (((g_Cl * (((0.5 * __k_2_v) + v) - (volt * NernstPotential(_clip(C0_Cl_E + (_lio_16 * ((0.5 * __k_2_n_Cl_E) + n_Cl_E)), false, _lio_17), _clip(C0_Cl_N + (_lio_18 * ((0.5 * __k_2_n_Cl_N) + n_Cl_N)), false, _lio_17), true, T_exp)))) + (g_K * (((0.5 * __k_2_n) + n) * (((0.5 * __k_2_v) + v) - (volt * NernstPotential(_clip(C0_K_E + (_lio_16 * ((0.5 * __k_2_n_K_E) + n_K_E)), false, _lio_20), _clip(C0_K_N + (_lio_18 * ((0.5 * __k_2_n_K_N) + n_K_N)), false, _lio_20), true, T_exp)))))) + (g_Na * ((((0.5 * __k_2_h) + h) * (_brian_pow((0.5 * __k_2_m) + m, 3))) * (((0.5 * __k_2_v) + v) - (volt * NernstPotential(_clip(C0_Na_E + (_lio_16 * ((0.5 * __k_2_n_Na_E) + n_Na_E)), false, _lio_23), _clip(C0_Na_N + (_lio_18 * ((0.5 * __k_2_n_Na_N) + n_Na_N)), false, _lio_23), true, T_exp)))))));
        const double __k_4_h = _lio_1 * (((1.0f*(_lio_2 * (_lio_3 + ((0.015 * __k_3_v) + (0.015 * v))))/expm1(_lio_4 * (_lio_5 + (__k_3_v + v)))) + (1.0f*(_lio_2 * (_lio_6 - ((0.015 * __k_3_v) + (0.015 * v))))/expm1(_lio_4 * (U_h - (__k_3_v + v))))) * (((- __k_3_h) + (1.0f*(_lio_2 * (_lio_3 + ((0.015 * __k_3_v) + (0.015 * v))))/(((1.0f*(_lio_2 * (_lio_3 + ((0.015 * __k_3_v) + (0.015 * v))))/expm1(_lio_4 * (_lio_5 + (__k_3_v + v)))) + (1.0f*(_lio_2 * (_lio_6 - ((0.015 * __k_3_v) + (0.015 * v))))/expm1(_lio_4 * (U_h - (__k_3_v + v))))) * expm1(_lio_4 * (_lio_5 + (__k_3_v + v)))))) - h));
        const double __k_4_m = _lio_1 * (((1.0f*(_lio_2 * (_lio_8 + ((0.124 * __k_3_v) + (0.124 * v))))/expm1(_lio_9 * (_lio_10 + (__k_3_v + v)))) + (1.0f*(_lio_2 * (_lio_7 - ((0.182 * __k_3_v) + (0.182 * v))))/expm1(_lio_9 * (U_m - (__k_3_v + v))))) * (((- __k_3_m) + (1.0f*(_lio_2 * (_lio_7 - ((0.182 * __k_3_v) + (0.182 * v))))/(((1.0f*(_lio_2 * (_lio_8 + ((0.124 * __k_3_v) + (0.124 * v))))/expm1(_lio_9 * (_lio_10 + (__k_3_v + v)))) + (1.0f*(_lio_2 * (_lio_7 - ((0.182 * __k_3_v) + (0.182 * v))))/expm1(_lio_9 * (U_m - (__k_3_v + v))))) * expm1(_lio_9 * (U_m - (__k_3_v + v)))))) - m));
        const double __k_4_n = _lio_11 * ((true + exp(_lio_13 * ((_lio_27 + ((- 0.0226551880380607) * __k_3_v)) - (0.0226551880380607 * v)))) * (((- __k_3_n) + (1.0f*true/(true + exp(_lio_12 * (U_n - (__k_3_v + v)))))) - n));
        const double __k_4_n_Cl_E = _lio_15 * ((g_Cl * ((__k_3_v + v) - (volt * NernstPotential(_clip(C0_Cl_E + (_lio_16 * (__k_3_n_Cl_E + n_Cl_E)), false, _lio_17), _clip(C0_Cl_N + (_lio_18 * (__k_3_n_Cl_N + n_Cl_N)), false, _lio_17), true, T_exp)))) + (g_KCC * ((_lio_19 * NernstPotential(_clip(C0_Cl_E + (_lio_16 * (__k_3_n_Cl_E + n_Cl_E)), false, _lio_17), _clip(C0_Cl_N + (_lio_18 * (__k_3_n_Cl_N + n_Cl_N)), false, _lio_17), true, T_exp)) + (volt * NernstPotential(_clip(C0_K_E + (_lio_16 * (__k_3_n_K_E + n_K_E)), false, _lio_20), _clip(C0_K_N + (_lio_18 * (__k_3_n_K_N + n_K_N)), false, _lio_20), true, T_exp)))));
        const double __k_4_n_Cl_N = _lio_21 * ((g_Cl * ((__k_3_v + v) - (volt * NernstPotential(_clip(C0_Cl_E + (_lio_16 * (__k_3_n_Cl_E + n_Cl_E)), false, _lio_17), _clip(C0_Cl_N + (_lio_18 * (__k_3_n_Cl_N + n_Cl_N)), false, _lio_17), true, T_exp)))) + (g_KCC * ((_lio_19 * NernstPotential(_clip(C0_Cl_E + (_lio_16 * (__k_3_n_Cl_E + n_Cl_E)), false, _lio_17), _clip(C0_Cl_N + (_lio_18 * (__k_3_n_Cl_N + n_Cl_N)), false, _lio_17), true, T_exp)) + (volt * NernstPotential(_clip(C0_K_E + (_lio_16 * (__k_3_n_K_E + n_K_E)), false, _lio_20), _clip(C0_K_N + (_lio_18 * (__k_3_n_K_N + n_K_N)), false, _lio_20), true, T_exp)))));
        const double __k_4_n_K_E = _lio_21 * ((g_K * ((__k_3_n + n) * ((__k_3_v + v) - (volt * NernstPotential(_clip(C0_K_E + (_lio_16 * (__k_3_n_K_E + n_K_E)), false, _lio_20), _clip(C0_K_N + (_lio_18 * (__k_3_n_K_N + n_K_N)), false, _lio_20), true, T_exp))))) - (g_KCC * ((_lio_19 * NernstPotential(_clip(C0_Cl_E + (_lio_16 * (__k_3_n_Cl_E + n_Cl_E)), false, _lio_17), _clip(C0_Cl_N + (_lio_18 * (__k_3_n_Cl_N + n_Cl_N)), false, _lio_17), true, T_exp)) + (volt * NernstPotential(_clip(C0_K_E + (_lio_16 * (__k_3_n_K_E + n_K_E)), false, _lio_20), _clip(C0_K_N + (_lio_18 * (__k_3_n_K_N + n_K_N)), false, _lio_20), true, T_exp)))));
        const double __k_4_n_K_N = _lio_15 * ((g_K * ((__k_3_n + n) * ((__k_3_v + v) - (volt * NernstPotential(_clip(C0_K_E + (_lio_16 * (__k_3_n_K_E + n_K_E)), false, _lio_20), _clip(C0_K_N + (_lio_18 * (__k_3_n_K_N + n_K_N)), false, _lio_20), true, T_exp))))) - (g_KCC * ((_lio_19 * NernstPotential(_clip(C0_Cl_E + (_lio_16 * (__k_3_n_Cl_E + n_Cl_E)), false, _lio_17), _clip(C0_Cl_N + (_lio_18 * (__k_3_n_Cl_N + n_Cl_N)), false, _lio_17), true, T_exp)) + (volt * NernstPotential(_clip(C0_K_E + (_lio_16 * (__k_3_n_K_E + n_K_E)), false, _lio_20), _clip(C0_K_N + (_lio_18 * (__k_3_n_K_N + n_K_N)), false, _lio_20), true, T_exp)))));
        const double __k_4_n_Na_E = _lio_22 * (((__k_3_h + h) * (_brian_pow(__k_3_m + m, 3))) * ((__k_3_v + v) - (volt * NernstPotential(_clip(C0_Na_E + (_lio_16 * (__k_3_n_Na_E + n_Na_E)), false, _lio_23), _clip(C0_Na_N + (_lio_18 * (__k_3_n_Na_N + n_Na_N)), false, _lio_23), true, T_exp))));
        const double __k_4_n_Na_N = _lio_24 * (((__k_3_h + h) * (_brian_pow(__k_3_m + m, 3))) * ((__k_3_v + v) - (volt * NernstPotential(_clip(C0_Na_E + (_lio_16 * (__k_3_n_Na_E + n_Na_E)), false, _lio_23), _clip(C0_Na_N + (_lio_18 * (__k_3_n_Na_N + n_Na_N)), false, _lio_23), true, T_exp))));
        const double __k_4_v = _lio_25 * (_lio_29 - (((g_Cl * ((__k_3_v + v) - (volt * NernstPotential(_clip(C0_Cl_E + (_lio_16 * (__k_3_n_Cl_E + n_Cl_E)), false, _lio_17), _clip(C0_Cl_N + (_lio_18 * (__k_3_n_Cl_N + n_Cl_N)), false, _lio_17), true, T_exp)))) + (g_K * ((__k_3_n + n) * ((__k_3_v + v) - (volt * NernstPotential(_clip(C0_K_E + (_lio_16 * (__k_3_n_K_E + n_K_E)), false, _lio_20), _clip(C0_K_N + (_lio_18 * (__k_3_n_K_N + n_K_N)), false, _lio_20), true, T_exp)))))) + (g_Na * (((__k_3_h + h) * (_brian_pow(__k_3_m + m, 3))) * ((__k_3_v + v) - (volt * NernstPotential(_clip(C0_Na_E + (_lio_16 * (__k_3_n_Na_E + n_Na_E)), false, _lio_23), _clip(C0_Na_N + (_lio_18 * (__k_3_n_Na_N + n_Na_N)), false, _lio_23), true, T_exp)))))));
        const double _h = ((((0.16666666666666666 * __k_1_h) + (0.3333333333333333 * __k_2_h)) + (0.3333333333333333 * __k_3_h)) + (0.16666666666666666 * __k_4_h)) + h;
        const double _m = ((((0.16666666666666666 * __k_1_m) + (0.3333333333333333 * __k_2_m)) + (0.3333333333333333 * __k_3_m)) + (0.16666666666666666 * __k_4_m)) + m;
        const double _n = ((((0.16666666666666666 * __k_1_n) + (0.3333333333333333 * __k_2_n)) + (0.3333333333333333 * __k_3_n)) + (0.16666666666666666 * __k_4_n)) + n;
        const double _n_Cl_E = ((((0.16666666666666666 * __k_1_n_Cl_E) + (0.3333333333333333 * __k_2_n_Cl_E)) + (0.3333333333333333 * __k_3_n_Cl_E)) + (0.16666666666666666 * __k_4_n_Cl_E)) + n_Cl_E;
        const double _n_Cl_N = ((((0.16666666666666666 * __k_1_n_Cl_N) + (0.3333333333333333 * __k_2_n_Cl_N)) + (0.3333333333333333 * __k_3_n_Cl_N)) + (0.16666666666666666 * __k_4_n_Cl_N)) + n_Cl_N;
        const double _n_K_E = ((((0.16666666666666666 * __k_1_n_K_E) + (0.3333333333333333 * __k_2_n_K_E)) + (0.3333333333333333 * __k_3_n_K_E)) + (0.16666666666666666 * __k_4_n_K_E)) + n_K_E;
        const double _n_K_N = ((((0.16666666666666666 * __k_1_n_K_N) + (0.3333333333333333 * __k_2_n_K_N)) + (0.3333333333333333 * __k_3_n_K_N)) + (0.16666666666666666 * __k_4_n_K_N)) + n_K_N;
        const double _n_Na_E = ((((0.16666666666666666 * __k_1_n_Na_E) + (0.3333333333333333 * __k_2_n_Na_E)) + (0.3333333333333333 * __k_3_n_Na_E)) + (0.16666666666666666 * __k_4_n_Na_E)) + n_Na_E;
        const double _n_Na_N = ((((0.16666666666666666 * __k_1_n_Na_N) + (0.3333333333333333 * __k_2_n_Na_N)) + (0.3333333333333333 * __k_3_n_Na_N)) + (0.16666666666666666 * __k_4_n_Na_N)) + n_Na_N;
        const double _v = ((((0.16666666666666666 * __k_1_v) + (0.3333333333333333 * __k_2_v)) + (0.3333333333333333 * __k_3_v)) + (0.16666666666666666 * __k_4_v)) + v;
        h = _h;
        m = _m;
        n = _n;
        n_Cl_E = _n_Cl_E;
        n_Cl_N = _n_Cl_N;
        n_K_E = _n_K_E;
        n_K_N = _n_K_N;
        n_Na_E = _n_Na_E;
        n_Na_N = _n_Na_N;
        v = _v;
        _ptr_array_hh_h[_idx] = h;
        _ptr_array_hh_m[_idx] = m;
        _ptr_array_hh_n[_idx] = n;
        _ptr_array_hh_n_Cl_E[_idx] = n_Cl_E;
        _ptr_array_hh_n_Cl_N[_idx] = n_Cl_N;
        _ptr_array_hh_n_K_E[_idx] = n_K_E;
        _ptr_array_hh_n_K_N[_idx] = n_K_N;
        _ptr_array_hh_n_Na_E[_idx] = n_Na_E;
        _ptr_array_hh_n_Na_N[_idx] = n_Na_N;
        _ptr_array_hh_v[_idx] = v;

    }

}


