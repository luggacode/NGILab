#include "code_objects/neumon_1_codeobject.h"
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



void _run_neumon_1_codeobject()
{
    using namespace brian;


    ///// CONSTANTS ///////////
    const size_t _numN = 1;
const size_t _num_clock_t = 1;
const size_t _num_indices = 1;
const size_t _num_source_h = 1;
const size_t _num_source_m = 1;
const size_t _num_source_n = 1;
const size_t _num_source_v = 1;
double* const _array_neumon_1_t = _dynamic_array_neumon_1_t.empty()? 0 : &_dynamic_array_neumon_1_t[0];
const size_t _numt = _dynamic_array_neumon_1_t.size();
    ///// POINTERS ////////////
        
    int32_t*   _ptr_array_neumon_1_N = _array_neumon_1_N;
    double*   _ptr_array_neuron_1_clock_t = _array_neuron_1_clock_t;
    int32_t*   _ptr_array_neumon_1__indices = _array_neumon_1__indices;
    double*   _ptr_array_neuron_1_h = _array_neuron_1_h;
    double*   _ptr_array_neuron_1_m = _array_neuron_1_m;
    double*   _ptr_array_neuron_1_n = _array_neuron_1_n;
    double*   _ptr_array_neuron_1_v = _array_neuron_1_v;
    double* __restrict  _ptr_array_neumon_1_t = _array_neumon_1_t;


    _dynamic_array_neumon_1_t.push_back(_ptr_array_neuron_1_clock_t[0]);

    const size_t _new_size = _dynamic_array_neumon_1_t.size();
    // Resize the dynamic arrays
    _dynamic_array_neumon_1_h.resize(_new_size, _num_indices);
    _dynamic_array_neumon_1_m.resize(_new_size, _num_indices);
    _dynamic_array_neumon_1_n.resize(_new_size, _num_indices);
    _dynamic_array_neumon_1_v.resize(_new_size, _num_indices);

    // scalar code
    const size_t _vectorisation_idx = -1;
        


    #pragma omp parallel for schedule(static)
    for (int _i = 0; _i < (int)_num_indices; _i++)
    {
        // vector code
        const size_t _idx = _ptr_array_neumon_1__indices[_i];
        const size_t _vectorisation_idx = _idx;
                
        const double _source_h = _ptr_array_neuron_1_h[_idx];
        const double _source_m = _ptr_array_neuron_1_m[_idx];
        const double _source_n = _ptr_array_neuron_1_n[_idx];
        const double _source_v = _ptr_array_neuron_1_v[_idx];
        const double _to_record_v = _source_v;
        const double _to_record_m = _source_m;
        const double _to_record_h = _source_h;
        const double _to_record_n = _source_n;


        _dynamic_array_neumon_1_h(_new_size-1, _i) = _to_record_h;
        _dynamic_array_neumon_1_m(_new_size-1, _i) = _to_record_m;
        _dynamic_array_neumon_1_n(_new_size-1, _i) = _to_record_n;
        _dynamic_array_neumon_1_v(_new_size-1, _i) = _to_record_v;
    }

    _ptr_array_neumon_1_N[0] = _new_size;


}


