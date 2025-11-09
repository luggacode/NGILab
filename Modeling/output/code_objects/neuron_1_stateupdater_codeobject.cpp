#include "code_objects/neuron_1_stateupdater_codeobject.h"
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



void _run_neuron_1_stateupdater_codeobject()
{
    using namespace brian;


    ///// CONSTANTS ///////////
    const double C = 0.009999999999999998;
const double EK = - 0.1;
const double ENa = 0.07100000000000001;
const double El = - 0.08700000000000001;
const size_t _numIin = 1;
const int64_t N = true;
const double Q10 = 2.3;
const int64_t Texp = 37;
const double Uh = - 0.066;
const double Um = - 0.038;
const double Un = 0.0187;
const double Wh = 0.006;
const double Wm = 0.006;
const double Wn = 0.0097;
const size_t _numdt = 1;
const double gK = 69300.0;
const double gNa = 20400.0;
const double gl = 0.3379999999999999;
const size_t _numh = 1;
const size_t _numm = 1;
const double mV = 0.001;
const double mvolt = 0.001;
const size_t _numn = 1;
const double second = true;
const size_t _numv = 1;
    ///// POINTERS ////////////
        
    double*   _ptr_array_neuron_1_Iin = _array_neuron_1_Iin;
    double*   _ptr_array_neuron_1_clock_dt = _array_neuron_1_clock_dt;
    double*   _ptr_array_neuron_1_h = _array_neuron_1_h;
    double*   _ptr_array_neuron_1_m = _array_neuron_1_m;
    double*   _ptr_array_neuron_1_n = _array_neuron_1_n;
    double*   _ptr_array_neuron_1_v = _array_neuron_1_v;


    //// MAIN CODE ////////////
    // scalar code
    const size_t _vectorisation_idx = -1;
        
    const double dt = _ptr_array_neuron_1_clock_dt[0];
    const double _lio_1 = ((- 225.0) * (_brian_pow(Q10, (0.1 * Texp) - 2.1))) * (_brian_pow(Uh, 2));
    const double _lio_2 = ((((- 15.0) * Uh) * (_brian_pow(mvolt, 2))) * (_brian_pow(second, 2))) * exp(1.0f*Uh/Wh);
    const double _lio_3 = (mvolt * second) * exp(1.0f*Um/Wh);
    const double _lio_4 = (((- mvolt) * second) * exp(1.0f*Uh/Wh)) * exp(1.0f*Um/Wh);
    const double _lio_5 = 1.0f*true/Wh;
    const double _lio_6 = 15.0 * ((Uh * (_brian_pow(mvolt, 2))) * (_brian_pow(second, 2)));
    const double _lio_7 = 15.0 * (((Uh * (_brian_pow(mvolt, 2))) * (_brian_pow(second, 2))) * exp(1.0f*Uh/Wh));
    const double _lio_8 = mvolt * second;
    const double _lio_9 = ((- mvolt) * second) * exp(1.0f*Uh/Wh);
    const double _lio_10 = (mvolt * second) * exp(1.0f*(- Um)/Wh);
    const double _lio_11 = 15.0 * (((_brian_pow(mvolt, 2)) * (_brian_pow(second, 2))) * exp(1.0f*Uh/Wh));
    const double _lio_12 = 15.0 * ((_brian_pow(mvolt, 2)) * (_brian_pow(second, 2)));
    const double _lio_13 = 1.0f*2.0/Wh;
    const double _lio_14 = 225.0 * ((_brian_pow(Q10, (0.1 * Texp) - 2.1)) * (_brian_pow(Uh, 2)));
    const double _lio_15 = (((- 15.0) * Uh) * (_brian_pow(mvolt, 2))) * (_brian_pow(second, 2));
    const double _lio_16 = (mvolt * second) * exp(1.0f*(2.0 * Um)/Wh);
    const double _lio_17 = (((- mvolt) * second) * exp(1.0f*Uh/Wh)) * exp(1.0f*(2.0 * Um)/Wh);
    const double _lio_18 = 30.0 * ((Uh * (_brian_pow(mvolt, 2))) * (_brian_pow(second, 2)));
    const double _lio_19 = 30.0 * ((_brian_pow(mvolt, 2)) * (_brian_pow(second, 2)));
    const double _lio_20 = 450.0 * ((_brian_pow(Q10, (0.1 * Texp) - 2.1)) * Uh);
    const double _lio_21 = 225.0 * (_brian_pow(Q10, (0.1 * Texp) - 2.1));
    const double _lio_22 = 15.0 * ((_brian_pow(Q10, (0.1 * Texp) - 2.1)) * Uh);
    const double _lio_23 = 15.0 * (_brian_pow(Q10, (0.1 * Texp) - 2.1));
    const double _lio_24 = ((- 22568.0) * (_brian_pow(Q10, (0.1 * Texp) - 2.1))) * (_brian_pow(Um, 2));
    const double _lio_25 = 182.0 * ((Um * (_brian_pow(mvolt, 2))) * (_brian_pow(second, 2)));
    const double _lio_26 = 1.0f*true/Wm;
    const double _lio_27 = (mvolt * second) * exp(1.0f*Um/Wm);
    const double _lio_28 = ((- mvolt) * second) * exp(1.0f*(2.0 * Um)/Wm);
    const double _lio_29 = 124.0 * ((Um * (_brian_pow(mvolt, 2))) * (_brian_pow(second, 2)));
    const double _lio_30 = (mvolt * second) * exp(1.0f*(- Um)/Wm);
    const double _lio_31 = 1.0f*2.0/Wm;
    const double _lio_32 = 182.0 * (((Um * (_brian_pow(mvolt, 2))) * (_brian_pow(second, 2))) * exp(1.0f*Um/Wm));
    const double _lio_33 = ((- mvolt) * second) * exp(1.0f*Um/Wm);
    const double _lio_34 = 182.0 * ((_brian_pow(mvolt, 2)) * (_brian_pow(second, 2)));
    const double _lio_35 = 124.0 * (((_brian_pow(mvolt, 2)) * (_brian_pow(second, 2))) * exp(1.0f*Um/Wm));
    const double _lio_36 = 124.0 * ((_brian_pow(mvolt, 2)) * (_brian_pow(second, 2)));
    const double _lio_37 = 124.0 * (((Um * (_brian_pow(mvolt, 2))) * (_brian_pow(second, 2))) * exp(1.0f*Um/Wm));
    const double _lio_38 = 182.0 * (((_brian_pow(mvolt, 2)) * (_brian_pow(second, 2))) * exp(1.0f*Um/Wm));
    const double _lio_39 = 33124.0 * ((_brian_pow(Q10, (0.1 * Texp) - 2.1)) * (_brian_pow(Um, 2)));
    const double _lio_40 = ((((- 182.0) * Um) * (_brian_pow(mvolt, 2))) * (_brian_pow(second, 2))) * exp(1.0f*(2.0 * Um)/Wm);
    const double _lio_41 = 124.0 * (((Um * (_brian_pow(mvolt, 2))) * (_brian_pow(second, 2))) * exp(1.0f*(2.0 * Um)/Wm));
    const double _lio_42 = 1.0f*3.0/Wm;
    const double _lio_43 = 364.0 * (((Um * (_brian_pow(mvolt, 2))) * (_brian_pow(second, 2))) * exp(1.0f*Um/Wm));
    const double _lio_44 = 182.0 * (((_brian_pow(mvolt, 2)) * (_brian_pow(second, 2))) * exp(1.0f*(2.0 * Um)/Wm));
    const double _lio_45 = 248.0 * (((_brian_pow(mvolt, 2)) * (_brian_pow(second, 2))) * exp(1.0f*Um/Wm));
    const double _lio_46 = 248.0 * (((Um * (_brian_pow(mvolt, 2))) * (_brian_pow(second, 2))) * exp(1.0f*Um/Wm));
    const double _lio_47 = 124.0 * (((_brian_pow(mvolt, 2)) * (_brian_pow(second, 2))) * exp(1.0f*(2.0 * Um)/Wm));
    const double _lio_48 = 364.0 * (((_brian_pow(mvolt, 2)) * (_brian_pow(second, 2))) * exp(1.0f*Um/Wm));
    const double _lio_49 = 45136.0 * ((_brian_pow(Q10, (0.1 * Texp) - 2.1)) * Um);
    const double _lio_50 = 33124.0 * (_brian_pow(Q10, (0.1 * Texp) - 2.1));
    const double _lio_51 = 66248.0 * ((_brian_pow(Q10, (0.1 * Texp) - 2.1)) * Um);
    const double _lio_52 = 22568.0 * (_brian_pow(Q10, (0.1 * Texp) - 2.1));
    const double _lio_53 = 182.0 * ((_brian_pow(Q10, (0.1 * Texp) - 2.1)) * Um);
    const double _lio_54 = 124.0 * (_brian_pow(Q10, (0.1 * Texp) - 2.1));
    const double _lio_55 = 124.0 * ((_brian_pow(Q10, (0.1 * Texp) - 2.1)) * Um);
    const double _lio_56 = 182.0 * (_brian_pow(Q10, (0.1 * Texp) - 2.1));
    const double _lio_57 = (((- true) * (_brian_pow(Q10, 2.1 - (0.1 * Texp)))) * (_brian_pow(Q10, (0.1 * Texp) - 2.1))) * second;
    const double _lio_58 = second * exp(1.0f*Un/Wn);
    const double _lio_59 = 1.0f*true/Wn;
    const double _lio_60 = 1.0f*(((- 900.399319534103) * (_brian_pow(Q10, (0.1 * Texp) - 2.1))) * dt)/second;
    const double _lio_61 = 1.0f*0.0226551880380607/mV;
    const double _lio_62 = 1.0f*(El * gl)/C;
    const double _lio_63 = 1.0f*(EK * gK)/C;
    const double _lio_64 = 1.0f*(ENa * gNa)/C;
    const double _lio_65 = 1.0f*true/C;
    const double _lio_66 = false - (1.0f*gl/C);
    const double _lio_67 = 1.0f*(- gK)/C;
    const double _lio_68 = 1.0f*gNa/C;


    const int _N = N;
    #pragma omp parallel for schedule(static)
    for(int _idx=0; _idx<_N; _idx++)
    {
        // vector code
        const size_t _vectorisation_idx = _idx;
                
        const double Iin = _ptr_array_neuron_1_Iin[_idx];
        double h = _ptr_array_neuron_1_h[_idx];
        double m = _ptr_array_neuron_1_m[_idx];
        double n = _ptr_array_neuron_1_n[_idx];
        double v = _ptr_array_neuron_1_v[_idx];
        const double _BA_h = 1.0f*(((((1.0f*_lio_1/((((((((((1.0f*_lio_2/(_lio_3 + (_lio_4 * exp(_lio_5 * (- v))))) + (1.0f*(_lio_6 * exp(_lio_5 * v))/(_lio_3 + (_lio_4 * exp(_lio_5 * (- v)))))) + (1.0f*_lio_7/(_lio_3 - (_lio_8 * exp(_lio_5 * v))))) + (1.0f*_lio_7/(_lio_9 + (_lio_8 * exp(_lio_5 * v))))) + (1.0f*_lio_6/(_lio_8 - (_lio_10 * exp(_lio_5 * v))))) + (1.0f*(_lio_11 * v)/(_lio_3 + (_lio_4 * exp(_lio_5 * (- v)))))) + (1.0f*(_lio_12 * v)/(_lio_8 + (_lio_9 * exp(_lio_5 * (- v)))))) + (1.0f*(_lio_11 * v)/((_lio_8 * exp(_lio_5 * v)) - (_lio_10 * exp(_lio_13 * v))))) + (1.0f*(_lio_12 * (v * exp(_lio_5 * v)))/(_lio_3 - (_lio_8 * exp(_lio_5 * v))))) - (((((((1.0f*_lio_6/(_lio_8 + (_lio_9 * exp(_lio_5 * (- v))))) + (1.0f*_lio_7/((_lio_8 * exp(_lio_5 * v)) - (_lio_10 * exp(_lio_13 * v))))) + (1.0f*(_lio_6 * exp(_lio_5 * v))/(_lio_3 - (_lio_8 * exp(_lio_5 * v))))) + (1.0f*(_lio_12 * (v * exp(_lio_5 * v)))/(_lio_3 + (_lio_4 * exp(_lio_5 * (- v)))))) + (1.0f*(_lio_11 * v)/(_lio_3 - (_lio_8 * exp(_lio_5 * v))))) + (1.0f*(_lio_11 * v)/(_lio_9 + (_lio_8 * exp(_lio_5 * v))))) + (1.0f*(_lio_12 * v)/(_lio_8 - (_lio_10 * exp(_lio_5 * v))))))) + (1.0f*_lio_14/((((((((1.0f*(_lio_15 * exp(_lio_13 * v))/(_lio_16 + (_lio_17 * exp(_lio_5 * (- v))))) + (1.0f*(_lio_18 * exp(_lio_5 * v))/(_lio_3 + (_lio_4 * exp(_lio_5 * (- v)))))) + (1.0f*(_lio_6 * exp(_lio_13 * v))/(_lio_16 - (_lio_3 * exp(_lio_5 * v))))) + (1.0f*_lio_6/(_lio_8 - (_lio_10 * exp(_lio_5 * v))))) + (1.0f*(_lio_12 * (v * exp(_lio_13 * v)))/(_lio_16 + (_lio_17 * exp(_lio_5 * (- v)))))) + (1.0f*(_lio_12 * v)/(_lio_8 + (_lio_9 * exp(_lio_5 * (- v)))))) + (1.0f*(_lio_19 * (v * exp(_lio_5 * v)))/(_lio_3 - (_lio_8 * exp(_lio_5 * v))))) - (((((1.0f*_lio_6/(_lio_8 + (_lio_9 * exp(_lio_5 * (- v))))) + (1.0f*(_lio_18 * exp(_lio_5 * v))/(_lio_3 - (_lio_8 * exp(_lio_5 * v))))) + (1.0f*(_lio_19 * (v * exp(_lio_5 * v)))/(_lio_3 + (_lio_4 * exp(_lio_5 * (- v)))))) + (1.0f*(_lio_12 * (v * exp(_lio_13 * v)))/(_lio_16 - (_lio_3 * exp(_lio_5 * v))))) + (1.0f*(_lio_12 * v)/(_lio_8 - (_lio_10 * exp(_lio_5 * v)))))))) + (1.0f*(_lio_20 * v)/((((((((((1.0f*_lio_2/(_lio_3 + (_lio_4 * exp(_lio_5 * (- v))))) + (1.0f*(_lio_6 * exp(_lio_5 * v))/(_lio_3 + (_lio_4 * exp(_lio_5 * (- v)))))) + (1.0f*_lio_7/(_lio_3 - (_lio_8 * exp(_lio_5 * v))))) + (1.0f*_lio_7/(_lio_9 + (_lio_8 * exp(_lio_5 * v))))) + (1.0f*_lio_6/(_lio_8 - (_lio_10 * exp(_lio_5 * v))))) + (1.0f*(_lio_11 * v)/(_lio_3 + (_lio_4 * exp(_lio_5 * (- v)))))) + (1.0f*(_lio_12 * v)/(_lio_8 + (_lio_9 * exp(_lio_5 * (- v)))))) + (1.0f*(_lio_11 * v)/((_lio_8 * exp(_lio_5 * v)) - (_lio_10 * exp(_lio_13 * v))))) + (1.0f*(_lio_12 * (v * exp(_lio_5 * v)))/(_lio_3 - (_lio_8 * exp(_lio_5 * v))))) - (((((((1.0f*_lio_6/(_lio_8 + (_lio_9 * exp(_lio_5 * (- v))))) + (1.0f*_lio_7/((_lio_8 * exp(_lio_5 * v)) - (_lio_10 * exp(_lio_13 * v))))) + (1.0f*(_lio_6 * exp(_lio_5 * v))/(_lio_3 - (_lio_8 * exp(_lio_5 * v))))) + (1.0f*(_lio_12 * (v * exp(_lio_5 * v)))/(_lio_3 + (_lio_4 * exp(_lio_5 * (- v)))))) + (1.0f*(_lio_11 * v)/(_lio_3 - (_lio_8 * exp(_lio_5 * v))))) + (1.0f*(_lio_11 * v)/(_lio_9 + (_lio_8 * exp(_lio_5 * v))))) + (1.0f*(_lio_12 * v)/(_lio_8 - (_lio_10 * exp(_lio_5 * v)))))))) + (1.0f*(_lio_21 * (_brian_pow(v, 2)))/((((((((1.0f*(_lio_15 * exp(_lio_13 * v))/(_lio_16 + (_lio_17 * exp(_lio_5 * (- v))))) + (1.0f*(_lio_18 * exp(_lio_5 * v))/(_lio_3 + (_lio_4 * exp(_lio_5 * (- v)))))) + (1.0f*(_lio_6 * exp(_lio_13 * v))/(_lio_16 - (_lio_3 * exp(_lio_5 * v))))) + (1.0f*_lio_6/(_lio_8 - (_lio_10 * exp(_lio_5 * v))))) + (1.0f*(_lio_12 * (v * exp(_lio_13 * v)))/(_lio_16 + (_lio_17 * exp(_lio_5 * (- v)))))) + (1.0f*(_lio_12 * v)/(_lio_8 + (_lio_9 * exp(_lio_5 * (- v)))))) + (1.0f*(_lio_19 * (v * exp(_lio_5 * v)))/(_lio_3 - (_lio_8 * exp(_lio_5 * v))))) - (((((1.0f*_lio_6/(_lio_8 + (_lio_9 * exp(_lio_5 * (- v))))) + (1.0f*(_lio_18 * exp(_lio_5 * v))/(_lio_3 - (_lio_8 * exp(_lio_5 * v))))) + (1.0f*(_lio_19 * (v * exp(_lio_5 * v)))/(_lio_3 + (_lio_4 * exp(_lio_5 * (- v)))))) + (1.0f*(_lio_12 * (v * exp(_lio_13 * v)))/(_lio_16 - (_lio_3 * exp(_lio_5 * v))))) + (1.0f*(_lio_12 * v)/(_lio_8 - (_lio_10 * exp(_lio_5 * v)))))))) - ((1.0f*(_lio_20 * v)/((((((((1.0f*(_lio_15 * exp(_lio_13 * v))/(_lio_16 + (_lio_17 * exp(_lio_5 * (- v))))) + (1.0f*(_lio_18 * exp(_lio_5 * v))/(_lio_3 + (_lio_4 * exp(_lio_5 * (- v)))))) + (1.0f*(_lio_6 * exp(_lio_13 * v))/(_lio_16 - (_lio_3 * exp(_lio_5 * v))))) + (1.0f*_lio_6/(_lio_8 - (_lio_10 * exp(_lio_5 * v))))) + (1.0f*(_lio_12 * (v * exp(_lio_13 * v)))/(_lio_16 + (_lio_17 * exp(_lio_5 * (- v)))))) + (1.0f*(_lio_12 * v)/(_lio_8 + (_lio_9 * exp(_lio_5 * (- v)))))) + (1.0f*(_lio_19 * (v * exp(_lio_5 * v)))/(_lio_3 - (_lio_8 * exp(_lio_5 * v))))) - (((((1.0f*_lio_6/(_lio_8 + (_lio_9 * exp(_lio_5 * (- v))))) + (1.0f*(_lio_18 * exp(_lio_5 * v))/(_lio_3 - (_lio_8 * exp(_lio_5 * v))))) + (1.0f*(_lio_19 * (v * exp(_lio_5 * v)))/(_lio_3 + (_lio_4 * exp(_lio_5 * (- v)))))) + (1.0f*(_lio_12 * (v * exp(_lio_13 * v)))/(_lio_16 - (_lio_3 * exp(_lio_5 * v))))) + (1.0f*(_lio_12 * v)/(_lio_8 - (_lio_10 * exp(_lio_5 * v))))))) + (1.0f*(_lio_21 * (_brian_pow(v, 2)))/((((((((((1.0f*_lio_2/(_lio_3 + (_lio_4 * exp(_lio_5 * (- v))))) + (1.0f*(_lio_6 * exp(_lio_5 * v))/(_lio_3 + (_lio_4 * exp(_lio_5 * (- v)))))) + (1.0f*_lio_7/(_lio_3 - (_lio_8 * exp(_lio_5 * v))))) + (1.0f*_lio_7/(_lio_9 + (_lio_8 * exp(_lio_5 * v))))) + (1.0f*_lio_6/(_lio_8 - (_lio_10 * exp(_lio_5 * v))))) + (1.0f*(_lio_11 * v)/(_lio_3 + (_lio_4 * exp(_lio_5 * (- v)))))) + (1.0f*(_lio_12 * v)/(_lio_8 + (_lio_9 * exp(_lio_5 * (- v)))))) + (1.0f*(_lio_11 * v)/((_lio_8 * exp(_lio_5 * v)) - (_lio_10 * exp(_lio_13 * v))))) + (1.0f*(_lio_12 * (v * exp(_lio_5 * v)))/(_lio_3 - (_lio_8 * exp(_lio_5 * v))))) - (((((((1.0f*_lio_6/(_lio_8 + (_lio_9 * exp(_lio_5 * (- v))))) + (1.0f*_lio_7/((_lio_8 * exp(_lio_5 * v)) - (_lio_10 * exp(_lio_13 * v))))) + (1.0f*(_lio_6 * exp(_lio_5 * v))/(_lio_3 - (_lio_8 * exp(_lio_5 * v))))) + (1.0f*(_lio_12 * (v * exp(_lio_5 * v)))/(_lio_3 + (_lio_4 * exp(_lio_5 * (- v)))))) + (1.0f*(_lio_11 * v)/(_lio_3 - (_lio_8 * exp(_lio_5 * v))))) + (1.0f*(_lio_11 * v)/(_lio_9 + (_lio_8 * exp(_lio_5 * v))))) + (1.0f*(_lio_12 * v)/(_lio_8 - (_lio_10 * exp(_lio_5 * v)))))))))/(((1.0f*_lio_22/(_lio_8 + (_lio_9 * exp(_lio_5 * (- v))))) + (1.0f*(_lio_23 * v)/(_lio_8 - (_lio_10 * exp(_lio_5 * v))))) - ((1.0f*_lio_22/(_lio_8 - (_lio_10 * exp(_lio_5 * v)))) + (1.0f*(_lio_23 * v)/(_lio_8 + (_lio_9 * exp(_lio_5 * (- v)))))));
        const double _h = (- _BA_h) + ((_BA_h + h) * exp(dt * (((1.0f*_lio_22/(_lio_8 + (_lio_9 * exp(_lio_5 * (- v))))) + (1.0f*(_lio_23 * v)/(_lio_8 - (_lio_10 * exp(_lio_5 * v))))) - ((1.0f*_lio_22/(_lio_8 - (_lio_10 * exp(_lio_5 * v)))) + (1.0f*(_lio_23 * v)/(_lio_8 + (_lio_9 * exp(_lio_5 * (- v)))))))));
        const double _BA_m = 1.0f*(((((1.0f*_lio_24/(((((((((1.0f*(_lio_25 * exp(_lio_26 * v))/(_lio_27 + (_lio_28 * exp(_lio_26 * (- v))))) + (1.0f*(_lio_29 * exp(_lio_26 * v))/((_lio_8 * exp(_lio_26 * v)) - (_lio_30 * exp(_lio_31 * v))))) + (1.0f*_lio_32/(_lio_33 + (_lio_8 * exp(_lio_26 * v))))) + (1.0f*_lio_29/(_lio_8 - (_lio_30 * exp(_lio_26 * v))))) + (1.0f*(_lio_34 * v)/(_lio_8 + (_lio_33 * exp(_lio_26 * (- v)))))) + (1.0f*(_lio_35 * v)/((_lio_8 * exp(_lio_26 * v)) - (_lio_30 * exp(_lio_31 * v))))) + (1.0f*(_lio_36 * (v * exp(_lio_26 * v)))/(_lio_27 - (_lio_8 * exp(_lio_26 * v))))) + (1.0f*(_lio_34 * (v * exp(_lio_26 * v)))/(_lio_33 + (_lio_8 * exp(_lio_26 * v))))) - ((((((((1.0f*_lio_25/(_lio_8 + (_lio_33 * exp(_lio_26 * (- v))))) + (1.0f*_lio_37/((_lio_8 * exp(_lio_26 * v)) - (_lio_30 * exp(_lio_31 * v))))) + (1.0f*(_lio_29 * exp(_lio_26 * v))/(_lio_27 - (_lio_8 * exp(_lio_26 * v))))) + (1.0f*(_lio_25 * exp(_lio_26 * v))/(_lio_33 + (_lio_8 * exp(_lio_26 * v))))) + (1.0f*(_lio_34 * (v * exp(_lio_26 * v)))/(_lio_27 + (_lio_28 * exp(_lio_26 * (- v)))))) + (1.0f*(_lio_36 * (v * exp(_lio_26 * v)))/((_lio_8 * exp(_lio_26 * v)) - (_lio_30 * exp(_lio_31 * v))))) + (1.0f*(_lio_38 * v)/(_lio_33 + (_lio_8 * exp(_lio_26 * v))))) + (1.0f*(_lio_36 * v)/(_lio_8 - (_lio_30 * exp(_lio_26 * v))))))) + (1.0f*_lio_39/((((((((1.0f*_lio_40/((_lio_33 * exp(_lio_26 * v)) + (_lio_8 * exp(_lio_31 * v)))) + (1.0f*_lio_41/((_lio_8 * exp(_lio_31 * v)) - (_lio_30 * exp(_lio_42 * v))))) + (1.0f*_lio_43/(_lio_33 + (_lio_8 * exp(_lio_26 * v))))) + (1.0f*_lio_29/(_lio_8 - (_lio_30 * exp(_lio_26 * v))))) + (1.0f*(_lio_44 * v)/((_lio_33 * exp(_lio_26 * v)) + (_lio_8 * exp(_lio_31 * v))))) + (1.0f*(_lio_34 * v)/(_lio_8 + (_lio_33 * exp(_lio_26 * (- v)))))) + (1.0f*(_lio_45 * v)/((_lio_8 * exp(_lio_26 * v)) - (_lio_30 * exp(_lio_31 * v))))) - (((((1.0f*_lio_25/(_lio_8 + (_lio_33 * exp(_lio_26 * (- v))))) + (1.0f*_lio_46/((_lio_8 * exp(_lio_26 * v)) - (_lio_30 * exp(_lio_31 * v))))) + (1.0f*(_lio_47 * v)/((_lio_8 * exp(_lio_31 * v)) - (_lio_30 * exp(_lio_42 * v))))) + (1.0f*(_lio_48 * v)/(_lio_33 + (_lio_8 * exp(_lio_26 * v))))) + (1.0f*(_lio_36 * v)/(_lio_8 - (_lio_30 * exp(_lio_26 * v)))))))) + (1.0f*(_lio_49 * v)/(((((((((1.0f*(_lio_25 * exp(_lio_26 * v))/(_lio_27 + (_lio_28 * exp(_lio_26 * (- v))))) + (1.0f*(_lio_29 * exp(_lio_26 * v))/((_lio_8 * exp(_lio_26 * v)) - (_lio_30 * exp(_lio_31 * v))))) + (1.0f*_lio_32/(_lio_33 + (_lio_8 * exp(_lio_26 * v))))) + (1.0f*_lio_29/(_lio_8 - (_lio_30 * exp(_lio_26 * v))))) + (1.0f*(_lio_34 * v)/(_lio_8 + (_lio_33 * exp(_lio_26 * (- v)))))) + (1.0f*(_lio_35 * v)/((_lio_8 * exp(_lio_26 * v)) - (_lio_30 * exp(_lio_31 * v))))) + (1.0f*(_lio_36 * (v * exp(_lio_26 * v)))/(_lio_27 - (_lio_8 * exp(_lio_26 * v))))) + (1.0f*(_lio_34 * (v * exp(_lio_26 * v)))/(_lio_33 + (_lio_8 * exp(_lio_26 * v))))) - ((((((((1.0f*_lio_25/(_lio_8 + (_lio_33 * exp(_lio_26 * (- v))))) + (1.0f*_lio_37/((_lio_8 * exp(_lio_26 * v)) - (_lio_30 * exp(_lio_31 * v))))) + (1.0f*(_lio_29 * exp(_lio_26 * v))/(_lio_27 - (_lio_8 * exp(_lio_26 * v))))) + (1.0f*(_lio_25 * exp(_lio_26 * v))/(_lio_33 + (_lio_8 * exp(_lio_26 * v))))) + (1.0f*(_lio_34 * (v * exp(_lio_26 * v)))/(_lio_27 + (_lio_28 * exp(_lio_26 * (- v)))))) + (1.0f*(_lio_36 * (v * exp(_lio_26 * v)))/((_lio_8 * exp(_lio_26 * v)) - (_lio_30 * exp(_lio_31 * v))))) + (1.0f*(_lio_38 * v)/(_lio_33 + (_lio_8 * exp(_lio_26 * v))))) + (1.0f*(_lio_36 * v)/(_lio_8 - (_lio_30 * exp(_lio_26 * v)))))))) + (1.0f*(_lio_50 * (_brian_pow(v, 2)))/((((((((1.0f*_lio_40/((_lio_33 * exp(_lio_26 * v)) + (_lio_8 * exp(_lio_31 * v)))) + (1.0f*_lio_41/((_lio_8 * exp(_lio_31 * v)) - (_lio_30 * exp(_lio_42 * v))))) + (1.0f*_lio_43/(_lio_33 + (_lio_8 * exp(_lio_26 * v))))) + (1.0f*_lio_29/(_lio_8 - (_lio_30 * exp(_lio_26 * v))))) + (1.0f*(_lio_44 * v)/((_lio_33 * exp(_lio_26 * v)) + (_lio_8 * exp(_lio_31 * v))))) + (1.0f*(_lio_34 * v)/(_lio_8 + (_lio_33 * exp(_lio_26 * (- v)))))) + (1.0f*(_lio_45 * v)/((_lio_8 * exp(_lio_26 * v)) - (_lio_30 * exp(_lio_31 * v))))) - (((((1.0f*_lio_25/(_lio_8 + (_lio_33 * exp(_lio_26 * (- v))))) + (1.0f*_lio_46/((_lio_8 * exp(_lio_26 * v)) - (_lio_30 * exp(_lio_31 * v))))) + (1.0f*(_lio_47 * v)/((_lio_8 * exp(_lio_31 * v)) - (_lio_30 * exp(_lio_42 * v))))) + (1.0f*(_lio_48 * v)/(_lio_33 + (_lio_8 * exp(_lio_26 * v))))) + (1.0f*(_lio_36 * v)/(_lio_8 - (_lio_30 * exp(_lio_26 * v)))))))) - ((1.0f*(_lio_51 * v)/((((((((1.0f*_lio_40/((_lio_33 * exp(_lio_26 * v)) + (_lio_8 * exp(_lio_31 * v)))) + (1.0f*_lio_41/((_lio_8 * exp(_lio_31 * v)) - (_lio_30 * exp(_lio_42 * v))))) + (1.0f*_lio_43/(_lio_33 + (_lio_8 * exp(_lio_26 * v))))) + (1.0f*_lio_29/(_lio_8 - (_lio_30 * exp(_lio_26 * v))))) + (1.0f*(_lio_44 * v)/((_lio_33 * exp(_lio_26 * v)) + (_lio_8 * exp(_lio_31 * v))))) + (1.0f*(_lio_34 * v)/(_lio_8 + (_lio_33 * exp(_lio_26 * (- v)))))) + (1.0f*(_lio_45 * v)/((_lio_8 * exp(_lio_26 * v)) - (_lio_30 * exp(_lio_31 * v))))) - (((((1.0f*_lio_25/(_lio_8 + (_lio_33 * exp(_lio_26 * (- v))))) + (1.0f*_lio_46/((_lio_8 * exp(_lio_26 * v)) - (_lio_30 * exp(_lio_31 * v))))) + (1.0f*(_lio_47 * v)/((_lio_8 * exp(_lio_31 * v)) - (_lio_30 * exp(_lio_42 * v))))) + (1.0f*(_lio_48 * v)/(_lio_33 + (_lio_8 * exp(_lio_26 * v))))) + (1.0f*(_lio_36 * v)/(_lio_8 - (_lio_30 * exp(_lio_26 * v))))))) + (1.0f*(_lio_52 * (_brian_pow(v, 2)))/(((((((((1.0f*(_lio_25 * exp(_lio_26 * v))/(_lio_27 + (_lio_28 * exp(_lio_26 * (- v))))) + (1.0f*(_lio_29 * exp(_lio_26 * v))/((_lio_8 * exp(_lio_26 * v)) - (_lio_30 * exp(_lio_31 * v))))) + (1.0f*_lio_32/(_lio_33 + (_lio_8 * exp(_lio_26 * v))))) + (1.0f*_lio_29/(_lio_8 - (_lio_30 * exp(_lio_26 * v))))) + (1.0f*(_lio_34 * v)/(_lio_8 + (_lio_33 * exp(_lio_26 * (- v)))))) + (1.0f*(_lio_35 * v)/((_lio_8 * exp(_lio_26 * v)) - (_lio_30 * exp(_lio_31 * v))))) + (1.0f*(_lio_36 * (v * exp(_lio_26 * v)))/(_lio_27 - (_lio_8 * exp(_lio_26 * v))))) + (1.0f*(_lio_34 * (v * exp(_lio_26 * v)))/(_lio_33 + (_lio_8 * exp(_lio_26 * v))))) - ((((((((1.0f*_lio_25/(_lio_8 + (_lio_33 * exp(_lio_26 * (- v))))) + (1.0f*_lio_37/((_lio_8 * exp(_lio_26 * v)) - (_lio_30 * exp(_lio_31 * v))))) + (1.0f*(_lio_29 * exp(_lio_26 * v))/(_lio_27 - (_lio_8 * exp(_lio_26 * v))))) + (1.0f*(_lio_25 * exp(_lio_26 * v))/(_lio_33 + (_lio_8 * exp(_lio_26 * v))))) + (1.0f*(_lio_34 * (v * exp(_lio_26 * v)))/(_lio_27 + (_lio_28 * exp(_lio_26 * (- v)))))) + (1.0f*(_lio_36 * (v * exp(_lio_26 * v)))/((_lio_8 * exp(_lio_26 * v)) - (_lio_30 * exp(_lio_31 * v))))) + (1.0f*(_lio_38 * v)/(_lio_33 + (_lio_8 * exp(_lio_26 * v))))) + (1.0f*(_lio_36 * v)/(_lio_8 - (_lio_30 * exp(_lio_26 * v)))))))))/(((1.0f*_lio_53/(_lio_8 + (_lio_33 * exp(_lio_26 * (- v))))) + (1.0f*(_lio_54 * v)/(_lio_8 - (_lio_30 * exp(_lio_26 * v))))) - ((1.0f*_lio_55/(_lio_8 - (_lio_30 * exp(_lio_26 * v)))) + (1.0f*(_lio_56 * v)/(_lio_8 + (_lio_33 * exp(_lio_26 * (- v)))))));
        const double _m = (- _BA_m) + ((_BA_m + m) * exp(dt * (((1.0f*_lio_53/(_lio_8 + (_lio_33 * exp(_lio_26 * (- v))))) + (1.0f*(_lio_54 * v)/(_lio_8 - (_lio_30 * exp(_lio_26 * v))))) - ((1.0f*_lio_55/(_lio_8 - (_lio_30 * exp(_lio_26 * v)))) + (1.0f*(_lio_56 * v)/(_lio_8 + (_lio_33 * exp(_lio_26 * (- v)))))))));
        const double _BA_n = 1.0f*_lio_57/(second + (_lio_58 * exp(_lio_59 * (- v))));
        const double _n = (- _BA_n) + ((_BA_n + n) * exp(_lio_60 * exp(_lio_61 * v)));
        const double _BA_v = 1.0f*(_lio_62 + (((_lio_63 * n) + (_lio_64 * (h * (_brian_pow(m, 3))))) + (_lio_65 * Iin)))/((_lio_66 + (_lio_67 * n)) - (_lio_68 * (h * (_brian_pow(m, 3)))));
        const double _v = (- _BA_v) + ((_BA_v + v) * exp(dt * ((_lio_66 + (_lio_67 * n)) - (_lio_68 * (h * (_brian_pow(m, 3)))))));
        h = _h;
        m = _m;
        n = _n;
        v = _v;
        _ptr_array_neuron_1_h[_idx] = h;
        _ptr_array_neuron_1_m[_idx] = m;
        _ptr_array_neuron_1_n[_idx] = n;
        _ptr_array_neuron_1_v[_idx] = v;

    }

}


