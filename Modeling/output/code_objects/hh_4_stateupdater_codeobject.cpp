#include "code_objects/hh_4_stateupdater_codeobject.h"
#include "objects.h"
#include "brianlib/common_math.h"
#include "brianlib/stdint_compat.h"
#include<cmath>
#include<ctime>
#include<iostream>
#include<fstream>
#include<climits>

////// SUPPORT CODE ///////
namespace {
        
    static double* _namespace_timedarray_5_values;
    static inline double _timedarray_5(const double t)
    {
        const double epsilon = 0.050000000000000003 / 67108864;
        int i = (int)((t/epsilon + 0.5)/67108864);
        if(i < 0)
           i = 0;
        if(i >= 4)
            i = 4-1;
        return _namespace_timedarray_5_values[i];
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



void _run_hh_4_stateupdater_codeobject()
{
    using namespace brian;


    ///// CONSTANTS ///////////
    const double E_Cl = - 0.08707812398837396;
const double E_K = - 0.10073080017019927;
const double E_Na = 0.07147113197557767;
const double I_max = false;
const int64_t N = true;
const double T_adj = 3.7910938669360985;
const double U_h = - 0.066;
const double U_m = - 0.038;
const double U_n = 0.0187;
const double W_h = 0.006;
const double W_m = 0.006;
const double W_n = 0.0097;
const double c_m = 0.009999999999999998;
const size_t _numdt = 1;
const double g_Cl_L = 0.49999999999999994;
const double g_K = 6929.999999999999;
const double g_Na = 20400.0;
const size_t _numh = 1;
const size_t _numm = 1;
const double ms = 0.001;
const double mvolt = 0.001;
const size_t _numn = 1;
const double second = true;
const size_t _numt = 1;
const size_t _numv = 1;
    ///// POINTERS ////////////
        
    double*   _ptr_array_hh_4_clock_dt = _array_hh_4_clock_dt;
    double*   _ptr_array_hh_4_h = _array_hh_4_h;
    double*   _ptr_array_hh_4_m = _array_hh_4_m;
    double*   _ptr_array_hh_4_n = _array_hh_4_n;
    double*   _ptr_array_hh_4_clock_t = _array_hh_4_clock_t;
    double*   _ptr_array_hh_4_v = _array_hh_4_v;
    _namespace_timedarray_5_values = _timedarray_5_values;


    //// MAIN CODE ////////////
    // scalar code
    const size_t _vectorisation_idx = -1;
        
    const double dt = _ptr_array_hh_4_clock_dt[0];
    const double t = _ptr_array_hh_4_clock_t[0];
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
    const double _lio_15 = 1.0f*dt/c_m;
    const double _lio_16 = I_max * _timedarray_5(t);
    const double _lio_17 = - E_Cl;
    const double _lio_18 = - E_K;
    const double _lio_19 = - E_Na;
    const double _lio_20 = false - (1.28137743543271 * mvolt);
    const double _lio_21 = I_max * _timedarray_5((0.5 * dt) + t);
    const double _lio_22 = I_max * _timedarray_5(dt + t);


    const int _N = N;
    #pragma omp parallel for schedule(static)
    for(int _idx=0; _idx<_N; _idx++)
    {
        // vector code
        const size_t _vectorisation_idx = _idx;
                
        double h = _ptr_array_hh_4_h[_idx];
        double m = _ptr_array_hh_4_m[_idx];
        double n = _ptr_array_hh_4_n[_idx];
        double v = _ptr_array_hh_4_v[_idx];
        const double __k_1_h = _lio_1 * (((- h) + (1.0f*(_lio_2 * (_lio_3 + (0.015 * v)))/(((1.0f*(_lio_2 * (_lio_3 + (0.015 * v)))/expm1(_lio_4 * (_lio_5 + v))) + (1.0f*(_lio_2 * (_lio_6 - (0.015 * v)))/expm1(_lio_4 * (U_h - v)))) * expm1(_lio_4 * (_lio_5 + v))))) * ((1.0f*(_lio_2 * (_lio_3 + (0.015 * v)))/expm1(_lio_4 * (_lio_5 + v))) + (1.0f*(_lio_2 * (_lio_6 - (0.015 * v)))/expm1(_lio_4 * (U_h - v)))));
        const double __k_1_m = _lio_1 * (((- m) + (1.0f*(_lio_2 * (_lio_7 - (0.182 * v)))/(((1.0f*(_lio_2 * (_lio_8 + (0.124 * v)))/expm1(_lio_9 * (_lio_10 + v))) + (1.0f*(_lio_2 * (_lio_7 - (0.182 * v)))/expm1(_lio_9 * (U_m - v)))) * expm1(_lio_9 * (U_m - v))))) * ((1.0f*(_lio_2 * (_lio_8 + (0.124 * v)))/expm1(_lio_9 * (_lio_10 + v))) + (1.0f*(_lio_2 * (_lio_7 - (0.182 * v)))/expm1(_lio_9 * (U_m - v)))));
        const double __k_1_n = _lio_11 * (((- n) + (1.0f*true/(true + exp(_lio_12 * (U_n - v))))) * (true + exp(_lio_13 * (_lio_14 - (0.0226551880380607 * v)))));
        const double __k_1_v = _lio_15 * (_lio_16 - (((g_Cl_L * (_lio_17 + v)) + (g_K * (n * (_lio_18 + v)))) + (g_Na * ((h * (_brian_pow(m, 3))) * (_lio_19 + v)))));
        const double __k_2_h = _lio_1 * (((1.0f*(_lio_2 * (_lio_3 + ((0.0075 * __k_1_v) + (0.015 * v))))/expm1(_lio_4 * (_lio_5 + ((0.5 * __k_1_v) + v)))) + (1.0f*(_lio_2 * (_lio_6 - ((0.0075 * __k_1_v) + (0.015 * v))))/expm1(_lio_4 * (U_h - ((0.5 * __k_1_v) + v))))) * (((0.5 * (- __k_1_h)) + (1.0f*(_lio_2 * (_lio_3 + ((0.0075 * __k_1_v) + (0.015 * v))))/(((1.0f*(_lio_2 * (_lio_3 + ((0.0075 * __k_1_v) + (0.015 * v))))/expm1(_lio_4 * (_lio_5 + ((0.5 * __k_1_v) + v)))) + (1.0f*(_lio_2 * (_lio_6 - ((0.0075 * __k_1_v) + (0.015 * v))))/expm1(_lio_4 * (U_h - ((0.5 * __k_1_v) + v))))) * expm1(_lio_4 * (_lio_5 + ((0.5 * __k_1_v) + v)))))) - h));
        const double __k_2_m = _lio_1 * (((1.0f*(_lio_2 * (_lio_8 + ((0.062 * __k_1_v) + (0.124 * v))))/expm1(_lio_9 * (_lio_10 + ((0.5 * __k_1_v) + v)))) + (1.0f*(_lio_2 * (_lio_7 - ((0.091 * __k_1_v) + (0.182 * v))))/expm1(_lio_9 * (U_m - ((0.5 * __k_1_v) + v))))) * (((0.5 * (- __k_1_m)) + (1.0f*(_lio_2 * (_lio_7 - ((0.091 * __k_1_v) + (0.182 * v))))/(((1.0f*(_lio_2 * (_lio_8 + ((0.062 * __k_1_v) + (0.124 * v))))/expm1(_lio_9 * (_lio_10 + ((0.5 * __k_1_v) + v)))) + (1.0f*(_lio_2 * (_lio_7 - ((0.091 * __k_1_v) + (0.182 * v))))/expm1(_lio_9 * (U_m - ((0.5 * __k_1_v) + v))))) * expm1(_lio_9 * (U_m - ((0.5 * __k_1_v) + v)))))) - m));
        const double __k_2_n = _lio_11 * ((true + exp(_lio_13 * ((_lio_20 + ((- 0.0113275940190304) * __k_1_v)) - (0.0226551880380607 * v)))) * (((0.5 * (- __k_1_n)) + (1.0f*true/(true + exp(_lio_12 * (U_n - ((0.5 * __k_1_v) + v)))))) - n));
        const double __k_2_v = _lio_15 * (_lio_21 - (((g_Cl_L * (_lio_17 + ((0.5 * __k_1_v) + v))) + (g_K * (((0.5 * __k_1_n) + n) * (_lio_18 + ((0.5 * __k_1_v) + v))))) + (g_Na * ((((0.5 * __k_1_h) + h) * (_brian_pow((0.5 * __k_1_m) + m, 3))) * (_lio_19 + ((0.5 * __k_1_v) + v))))));
        const double __k_3_h = _lio_1 * (((1.0f*(_lio_2 * (_lio_3 + ((0.0075 * __k_2_v) + (0.015 * v))))/expm1(_lio_4 * (_lio_5 + ((0.5 * __k_2_v) + v)))) + (1.0f*(_lio_2 * (_lio_6 - ((0.0075 * __k_2_v) + (0.015 * v))))/expm1(_lio_4 * (U_h - ((0.5 * __k_2_v) + v))))) * (((0.5 * (- __k_2_h)) + (1.0f*(_lio_2 * (_lio_3 + ((0.0075 * __k_2_v) + (0.015 * v))))/(((1.0f*(_lio_2 * (_lio_3 + ((0.0075 * __k_2_v) + (0.015 * v))))/expm1(_lio_4 * (_lio_5 + ((0.5 * __k_2_v) + v)))) + (1.0f*(_lio_2 * (_lio_6 - ((0.0075 * __k_2_v) + (0.015 * v))))/expm1(_lio_4 * (U_h - ((0.5 * __k_2_v) + v))))) * expm1(_lio_4 * (_lio_5 + ((0.5 * __k_2_v) + v)))))) - h));
        const double __k_3_m = _lio_1 * (((1.0f*(_lio_2 * (_lio_8 + ((0.062 * __k_2_v) + (0.124 * v))))/expm1(_lio_9 * (_lio_10 + ((0.5 * __k_2_v) + v)))) + (1.0f*(_lio_2 * (_lio_7 - ((0.091 * __k_2_v) + (0.182 * v))))/expm1(_lio_9 * (U_m - ((0.5 * __k_2_v) + v))))) * (((0.5 * (- __k_2_m)) + (1.0f*(_lio_2 * (_lio_7 - ((0.091 * __k_2_v) + (0.182 * v))))/(((1.0f*(_lio_2 * (_lio_8 + ((0.062 * __k_2_v) + (0.124 * v))))/expm1(_lio_9 * (_lio_10 + ((0.5 * __k_2_v) + v)))) + (1.0f*(_lio_2 * (_lio_7 - ((0.091 * __k_2_v) + (0.182 * v))))/expm1(_lio_9 * (U_m - ((0.5 * __k_2_v) + v))))) * expm1(_lio_9 * (U_m - ((0.5 * __k_2_v) + v)))))) - m));
        const double __k_3_n = _lio_11 * ((true + exp(_lio_13 * ((_lio_20 + ((- 0.0113275940190304) * __k_2_v)) - (0.0226551880380607 * v)))) * (((0.5 * (- __k_2_n)) + (1.0f*true/(true + exp(_lio_12 * (U_n - ((0.5 * __k_2_v) + v)))))) - n));
        const double __k_3_v = _lio_15 * (_lio_21 - (((g_Cl_L * (_lio_17 + ((0.5 * __k_2_v) + v))) + (g_K * (((0.5 * __k_2_n) + n) * (_lio_18 + ((0.5 * __k_2_v) + v))))) + (g_Na * ((((0.5 * __k_2_h) + h) * (_brian_pow((0.5 * __k_2_m) + m, 3))) * (_lio_19 + ((0.5 * __k_2_v) + v))))));
        const double __k_4_h = _lio_1 * (((1.0f*(_lio_2 * (_lio_3 + ((0.015 * __k_3_v) + (0.015 * v))))/expm1(_lio_4 * (_lio_5 + (__k_3_v + v)))) + (1.0f*(_lio_2 * (_lio_6 - ((0.015 * __k_3_v) + (0.015 * v))))/expm1(_lio_4 * (U_h - (__k_3_v + v))))) * (((- __k_3_h) + (1.0f*(_lio_2 * (_lio_3 + ((0.015 * __k_3_v) + (0.015 * v))))/(((1.0f*(_lio_2 * (_lio_3 + ((0.015 * __k_3_v) + (0.015 * v))))/expm1(_lio_4 * (_lio_5 + (__k_3_v + v)))) + (1.0f*(_lio_2 * (_lio_6 - ((0.015 * __k_3_v) + (0.015 * v))))/expm1(_lio_4 * (U_h - (__k_3_v + v))))) * expm1(_lio_4 * (_lio_5 + (__k_3_v + v)))))) - h));
        const double __k_4_m = _lio_1 * (((1.0f*(_lio_2 * (_lio_8 + ((0.124 * __k_3_v) + (0.124 * v))))/expm1(_lio_9 * (_lio_10 + (__k_3_v + v)))) + (1.0f*(_lio_2 * (_lio_7 - ((0.182 * __k_3_v) + (0.182 * v))))/expm1(_lio_9 * (U_m - (__k_3_v + v))))) * (((- __k_3_m) + (1.0f*(_lio_2 * (_lio_7 - ((0.182 * __k_3_v) + (0.182 * v))))/(((1.0f*(_lio_2 * (_lio_8 + ((0.124 * __k_3_v) + (0.124 * v))))/expm1(_lio_9 * (_lio_10 + (__k_3_v + v)))) + (1.0f*(_lio_2 * (_lio_7 - ((0.182 * __k_3_v) + (0.182 * v))))/expm1(_lio_9 * (U_m - (__k_3_v + v))))) * expm1(_lio_9 * (U_m - (__k_3_v + v)))))) - m));
        const double __k_4_n = _lio_11 * ((true + exp(_lio_13 * ((_lio_20 + ((- 0.0226551880380607) * __k_3_v)) - (0.0226551880380607 * v)))) * (((- __k_3_n) + (1.0f*true/(true + exp(_lio_12 * (U_n - (__k_3_v + v)))))) - n));
        const double __k_4_v = _lio_15 * (_lio_22 - (((g_Cl_L * (_lio_17 + (__k_3_v + v))) + (g_K * ((__k_3_n + n) * (_lio_18 + (__k_3_v + v))))) + (g_Na * (((__k_3_h + h) * (_brian_pow(__k_3_m + m, 3))) * (_lio_19 + (__k_3_v + v))))));
        const double _h = ((((0.16666666666666666 * __k_1_h) + (0.3333333333333333 * __k_2_h)) + (0.3333333333333333 * __k_3_h)) + (0.16666666666666666 * __k_4_h)) + h;
        const double _m = ((((0.16666666666666666 * __k_1_m) + (0.3333333333333333 * __k_2_m)) + (0.3333333333333333 * __k_3_m)) + (0.16666666666666666 * __k_4_m)) + m;
        const double _n = ((((0.16666666666666666 * __k_1_n) + (0.3333333333333333 * __k_2_n)) + (0.3333333333333333 * __k_3_n)) + (0.16666666666666666 * __k_4_n)) + n;
        const double _v = ((((0.16666666666666666 * __k_1_v) + (0.3333333333333333 * __k_2_v)) + (0.3333333333333333 * __k_3_v)) + (0.16666666666666666 * __k_4_v)) + v;
        h = _h;
        m = _m;
        n = _n;
        v = _v;
        _ptr_array_hh_4_h[_idx] = h;
        _ptr_array_hh_4_m[_idx] = m;
        _ptr_array_hh_4_n[_idx] = n;
        _ptr_array_hh_4_v[_idx] = v;

    }

}


