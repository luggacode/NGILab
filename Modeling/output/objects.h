
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
extern Clock hh_4_clock;
extern Clock svmon_clock_3;

//////////////// networks /////////////////
extern Network network_3;



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
extern int32_t *_array_hh_4__spikespace;
extern const int _num__array_hh_4__spikespace;
extern double *_array_hh_4_clock_dt;
extern const int _num__array_hh_4_clock_dt;
extern double *_array_hh_4_clock_t;
extern const int _num__array_hh_4_clock_t;
extern int64_t *_array_hh_4_clock_timestep;
extern const int _num__array_hh_4_clock_timestep;
extern double *_array_hh_4_h;
extern const int _num__array_hh_4_h;
extern int32_t *_array_hh_4_i;
extern const int _num__array_hh_4_i;
extern double *_array_hh_4_m;
extern const int _num__array_hh_4_m;
extern double *_array_hh_4_n;
extern const int _num__array_hh_4_n;
extern double *_array_hh_4_v;
extern const int _num__array_hh_4_v;
extern int32_t *_array_svmon__indices;
extern const int _num__array_svmon__indices;
extern double *_array_svmon_clock_3_dt;
extern const int _num__array_svmon_clock_3_dt;
extern double *_array_svmon_clock_3_t;
extern const int _num__array_svmon_clock_3_t;
extern int64_t *_array_svmon_clock_3_timestep;
extern const int _num__array_svmon_clock_3_timestep;
extern double *_array_svmon_h;
extern const int _num__array_svmon_h;
extern double *_array_svmon_m;
extern const int _num__array_svmon_m;
extern int32_t *_array_svmon_N;
extern const int _num__array_svmon_N;
extern double *_array_svmon_n;
extern const int _num__array_svmon_n;
extern double *_array_svmon_v;
extern const int _num__array_svmon_v;

//////////////// dynamic arrays 2d /////////
extern DynamicArray2D<double> _dynamic_array_svmon_h;
extern DynamicArray2D<double> _dynamic_array_svmon_m;
extern DynamicArray2D<double> _dynamic_array_svmon_n;
extern DynamicArray2D<double> _dynamic_array_svmon_v;

/////////////// static arrays /////////////
extern double *_timedarray_5_values;
extern const int _num__timedarray_5_values;

//////////////// synapses /////////////////

// Profiling information for each code object
}

void _init_arrays();
void _load_arrays();
void _write_arrays();
void _dealloc_arrays();

#endif


