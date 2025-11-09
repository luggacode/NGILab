
#ifndef _BRIAN_OBJECTS_H
#define _BRIAN_OBJECTS_H

#include "synapses_classes.h"
#include "brianlib/clocks.h"
#include "brianlib/dynamic_array.h"
#include "brianlib/stdint_compat.h"
#include "network.h"
#include "randomkit.h"
#include<vector>
#include <omp.h>

namespace brian {

extern std::string results_dir;
// In OpenMP we need one state per thread
extern std::vector< rk_state* > _mersenne_twister_states;

//////////////// clocks ///////////////////
extern Clock hh_1_clock;
extern Clock svmon_clock_1;

//////////////// networks /////////////////
extern Network network_1;



void set_variable_by_name(std::string, std::string);

//////////////// dynamic arrays ///////////
extern std::vector<double> _dynamic_array_svmon_t;

//////////////// arrays ///////////////////
extern double *_array_defaultclock_dt;
extern const int _num__array_defaultclock_dt;
extern double *_array_defaultclock_t;
extern const int _num__array_defaultclock_t;
extern int64_t *_array_defaultclock_timestep;
extern const int _num__array_defaultclock_timestep;
extern int32_t *_array_hh_1__spikespace;
extern const int _num__array_hh_1__spikespace;
extern double *_array_hh_1_clock_dt;
extern const int _num__array_hh_1_clock_dt;
extern double *_array_hh_1_clock_t;
extern const int _num__array_hh_1_clock_t;
extern int64_t *_array_hh_1_clock_timestep;
extern const int _num__array_hh_1_clock_timestep;
extern double *_array_hh_1_h;
extern const int _num__array_hh_1_h;
extern int32_t *_array_hh_1_i;
extern const int _num__array_hh_1_i;
extern double *_array_hh_1_m;
extern const int _num__array_hh_1_m;
extern double *_array_hh_1_n;
extern const int _num__array_hh_1_n;
extern double *_array_hh_1_n_Cl_E;
extern const int _num__array_hh_1_n_Cl_E;
extern double *_array_hh_1_n_Cl_N;
extern const int _num__array_hh_1_n_Cl_N;
extern double *_array_hh_1_n_K_E;
extern const int _num__array_hh_1_n_K_E;
extern double *_array_hh_1_n_K_N;
extern const int _num__array_hh_1_n_K_N;
extern double *_array_hh_1_n_Na_E;
extern const int _num__array_hh_1_n_Na_E;
extern double *_array_hh_1_n_Na_N;
extern const int _num__array_hh_1_n_Na_N;
extern double *_array_hh_1_v;
extern const int _num__array_hh_1_v;
extern int32_t *_array_svmon__indices;
extern const int _num__array_svmon__indices;
extern double *_array_svmon_C_Cl_N;
extern const int _num__array_svmon_C_Cl_N;
extern double *_array_svmon_C_K_N;
extern const int _num__array_svmon_C_K_N;
extern double *_array_svmon_C_Na_N;
extern const int _num__array_svmon_C_Na_N;
extern double *_array_svmon_clock_1_dt;
extern const int _num__array_svmon_clock_1_dt;
extern double *_array_svmon_clock_1_t;
extern const int _num__array_svmon_clock_1_t;
extern int64_t *_array_svmon_clock_1_timestep;
extern const int _num__array_svmon_clock_1_timestep;
extern double *_array_svmon_E_Cl;
extern const int _num__array_svmon_E_Cl;
extern double *_array_svmon_E_K;
extern const int _num__array_svmon_E_K;
extern double *_array_svmon_E_Na;
extern const int _num__array_svmon_E_Na;
extern double *_array_svmon_I_Cl_L;
extern const int _num__array_svmon_I_Cl_L;
extern double *_array_svmon_I_K;
extern const int _num__array_svmon_I_K;
extern double *_array_svmon_I_KCC;
extern const int _num__array_svmon_I_KCC;
extern double *_array_svmon_I_Na;
extern const int _num__array_svmon_I_Na;
extern double *_array_svmon_I_Na_L;
extern const int _num__array_svmon_I_Na_L;
extern double *_array_svmon_I_NKP;
extern const int _num__array_svmon_I_NKP;
extern int32_t *_array_svmon_N;
extern const int _num__array_svmon_N;
extern double *_array_svmon_v;
extern const int _num__array_svmon_v;

//////////////// dynamic arrays 2d /////////
extern DynamicArray2D<double> _dynamic_array_svmon_C_Cl_N;
extern DynamicArray2D<double> _dynamic_array_svmon_C_K_N;
extern DynamicArray2D<double> _dynamic_array_svmon_C_Na_N;
extern DynamicArray2D<double> _dynamic_array_svmon_E_Cl;
extern DynamicArray2D<double> _dynamic_array_svmon_E_K;
extern DynamicArray2D<double> _dynamic_array_svmon_E_Na;
extern DynamicArray2D<double> _dynamic_array_svmon_I_Cl_L;
extern DynamicArray2D<double> _dynamic_array_svmon_I_K;
extern DynamicArray2D<double> _dynamic_array_svmon_I_KCC;
extern DynamicArray2D<double> _dynamic_array_svmon_I_Na;
extern DynamicArray2D<double> _dynamic_array_svmon_I_Na_L;
extern DynamicArray2D<double> _dynamic_array_svmon_I_NKP;
extern DynamicArray2D<double> _dynamic_array_svmon_v;

/////////////// static arrays /////////////
extern double *_timedarray_1_values;
extern const int _num__timedarray_1_values;

//////////////// synapses /////////////////

// Profiling information for each code object
}

void _init_arrays();
void _load_arrays();
void _write_arrays();
void _dealloc_arrays();

#endif


