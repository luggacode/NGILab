
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
extern Clock NGexc_clock_1;
extern Clock NGinh_clock_1;
extern Clock svmon_exc_add_clock_1;
extern Clock svmon_inh_add_clock_1;

//////////////// networks /////////////////
extern Network network_1;



void set_variable_by_name(std::string, std::string);

//////////////// dynamic arrays ///////////
extern std::vector<int32_t> _dynamic_array_spikemonitor_2_i;
extern std::vector<double> _dynamic_array_spikemonitor_2_t;
extern std::vector<int32_t> _dynamic_array_spikemonitor_3_i;
extern std::vector<double> _dynamic_array_spikemonitor_3_t;
extern std::vector<double> _dynamic_array_svmon_exc_add_t;
extern std::vector<double> _dynamic_array_svmon_inh_add_t;

//////////////// arrays ///////////////////
extern double *_array_defaultclock_dt;
extern const int _num__array_defaultclock_dt;
extern double *_array_defaultclock_t;
extern const int _num__array_defaultclock_t;
extern int64_t *_array_defaultclock_timestep;
extern const int _num__array_defaultclock_timestep;
extern int32_t *_array_NGexc__spikespace;
extern const int _num__array_NGexc__spikespace;
extern double *_array_NGexc_clock_1_dt;
extern const int _num__array_NGexc_clock_1_dt;
extern double *_array_NGexc_clock_1_t;
extern const int _num__array_NGexc_clock_1_t;
extern int64_t *_array_NGexc_clock_1_timestep;
extern const int _num__array_NGexc_clock_1_timestep;
extern double *_array_NGexc_h;
extern const int _num__array_NGexc_h;
extern int32_t *_array_NGexc_i;
extern const int _num__array_NGexc_i;
extern double *_array_NGexc_I_max;
extern const int _num__array_NGexc_I_max;
extern double *_array_NGexc_lastspike;
extern const int _num__array_NGexc_lastspike;
extern double *_array_NGexc_m;
extern const int _num__array_NGexc_m;
extern double *_array_NGexc_mask_noise;
extern const int _num__array_NGexc_mask_noise;
extern double *_array_NGexc_n;
extern const int _num__array_NGexc_n;
extern char *_array_NGexc_not_refractory;
extern const int _num__array_NGexc_not_refractory;
extern double *_array_NGexc_p;
extern const int _num__array_NGexc_p;
extern int32_t *_array_NGexc_subgroup__sub_idx;
extern const int _num__array_NGexc_subgroup__sub_idx;
extern int32_t *_array_NGexc_subgroup__sub_idx_1;
extern const int _num__array_NGexc_subgroup__sub_idx_1;
extern int32_t *_array_NGexc_subgroup__sub_idx_2;
extern const int _num__array_NGexc_subgroup__sub_idx_2;
extern double *_array_NGexc_v;
extern const int _num__array_NGexc_v;
extern double *_array_NGexc_x;
extern const int _num__array_NGexc_x;
extern double *_array_NGexc_y;
extern const int _num__array_NGexc_y;
extern double *_array_NGexc_z;
extern const int _num__array_NGexc_z;
extern int32_t *_array_NGinh__spikespace;
extern const int _num__array_NGinh__spikespace;
extern double *_array_NGinh_clock_1_dt;
extern const int _num__array_NGinh_clock_1_dt;
extern double *_array_NGinh_clock_1_t;
extern const int _num__array_NGinh_clock_1_t;
extern int64_t *_array_NGinh_clock_1_timestep;
extern const int _num__array_NGinh_clock_1_timestep;
extern double *_array_NGinh_h;
extern const int _num__array_NGinh_h;
extern int32_t *_array_NGinh_i;
extern const int _num__array_NGinh_i;
extern double *_array_NGinh_I_max;
extern const int _num__array_NGinh_I_max;
extern double *_array_NGinh_lastspike;
extern const int _num__array_NGinh_lastspike;
extern double *_array_NGinh_m;
extern const int _num__array_NGinh_m;
extern double *_array_NGinh_mask_noise;
extern const int _num__array_NGinh_mask_noise;
extern double *_array_NGinh_n;
extern const int _num__array_NGinh_n;
extern char *_array_NGinh_not_refractory;
extern const int _num__array_NGinh_not_refractory;
extern double *_array_NGinh_p;
extern const int _num__array_NGinh_p;
extern int32_t *_array_NGinh_subgroup__sub_idx;
extern const int _num__array_NGinh_subgroup__sub_idx;
extern int32_t *_array_NGinh_subgroup__sub_idx_1;
extern const int _num__array_NGinh_subgroup__sub_idx_1;
extern int32_t *_array_NGinh_subgroup__sub_idx_2;
extern const int _num__array_NGinh_subgroup__sub_idx_2;
extern double *_array_NGinh_v;
extern const int _num__array_NGinh_v;
extern double *_array_NGinh_x;
extern const int _num__array_NGinh_x;
extern double *_array_NGinh_y;
extern const int _num__array_NGinh_y;
extern double *_array_NGinh_z;
extern const int _num__array_NGinh_z;
extern int32_t *_array_spikemonitor_2__source_idx;
extern const int _num__array_spikemonitor_2__source_idx;
extern int32_t *_array_spikemonitor_2_count;
extern const int _num__array_spikemonitor_2_count;
extern int32_t *_array_spikemonitor_2_N;
extern const int _num__array_spikemonitor_2_N;
extern int32_t *_array_spikemonitor_3__source_idx;
extern const int _num__array_spikemonitor_3__source_idx;
extern int32_t *_array_spikemonitor_3_count;
extern const int _num__array_spikemonitor_3_count;
extern int32_t *_array_spikemonitor_3_N;
extern const int _num__array_spikemonitor_3_N;
extern int32_t *_array_svmon_exc_add__indices;
extern const int _num__array_svmon_exc_add__indices;
extern double *_array_svmon_exc_add_a_m;
extern const int _num__array_svmon_exc_add_a_m;
extern double *_array_svmon_exc_add_b_m;
extern const int _num__array_svmon_exc_add_b_m;
extern double *_array_svmon_exc_add_clock_1_dt;
extern const int _num__array_svmon_exc_add_clock_1_dt;
extern double *_array_svmon_exc_add_clock_1_t;
extern const int _num__array_svmon_exc_add_clock_1_t;
extern int64_t *_array_svmon_exc_add_clock_1_timestep;
extern const int _num__array_svmon_exc_add_clock_1_timestep;
extern double *_array_svmon_exc_add_h;
extern const int _num__array_svmon_exc_add_h;
extern double *_array_svmon_exc_add_I_inj;
extern const int _num__array_svmon_exc_add_I_inj;
extern double *_array_svmon_exc_add_I_Kd;
extern const int _num__array_svmon_exc_add_I_Kd;
extern double *_array_svmon_exc_add_I_L;
extern const int _num__array_svmon_exc_add_I_L;
extern double *_array_svmon_exc_add_I_M;
extern const int _num__array_svmon_exc_add_I_M;
extern double *_array_svmon_exc_add_I_max;
extern const int _num__array_svmon_exc_add_I_max;
extern double *_array_svmon_exc_add_I_Na;
extern const int _num__array_svmon_exc_add_I_Na;
extern double *_array_svmon_exc_add_m;
extern const int _num__array_svmon_exc_add_m;
extern int32_t *_array_svmon_exc_add_N;
extern const int _num__array_svmon_exc_add_N;
extern double *_array_svmon_exc_add_n;
extern const int _num__array_svmon_exc_add_n;
extern double *_array_svmon_exc_add_p;
extern const int _num__array_svmon_exc_add_p;
extern double *_array_svmon_exc_add_p_inf;
extern const int _num__array_svmon_exc_add_p_inf;
extern double *_array_svmon_exc_add_tau_h;
extern const int _num__array_svmon_exc_add_tau_h;
extern double *_array_svmon_exc_add_tau_m;
extern const int _num__array_svmon_exc_add_tau_m;
extern double *_array_svmon_exc_add_tau_n;
extern const int _num__array_svmon_exc_add_tau_n;
extern double *_array_svmon_exc_add_tau_p;
extern const int _num__array_svmon_exc_add_tau_p;
extern double *_array_svmon_exc_add_v;
extern const int _num__array_svmon_exc_add_v;
extern int32_t *_array_svmon_inh_add__indices;
extern const int _num__array_svmon_inh_add__indices;
extern double *_array_svmon_inh_add_a_m;
extern const int _num__array_svmon_inh_add_a_m;
extern double *_array_svmon_inh_add_b_m;
extern const int _num__array_svmon_inh_add_b_m;
extern double *_array_svmon_inh_add_clock_1_dt;
extern const int _num__array_svmon_inh_add_clock_1_dt;
extern double *_array_svmon_inh_add_clock_1_t;
extern const int _num__array_svmon_inh_add_clock_1_t;
extern int64_t *_array_svmon_inh_add_clock_1_timestep;
extern const int _num__array_svmon_inh_add_clock_1_timestep;
extern double *_array_svmon_inh_add_h;
extern const int _num__array_svmon_inh_add_h;
extern double *_array_svmon_inh_add_I_inj;
extern const int _num__array_svmon_inh_add_I_inj;
extern double *_array_svmon_inh_add_I_Kd;
extern const int _num__array_svmon_inh_add_I_Kd;
extern double *_array_svmon_inh_add_I_L;
extern const int _num__array_svmon_inh_add_I_L;
extern double *_array_svmon_inh_add_I_M;
extern const int _num__array_svmon_inh_add_I_M;
extern double *_array_svmon_inh_add_I_max;
extern const int _num__array_svmon_inh_add_I_max;
extern double *_array_svmon_inh_add_I_Na;
extern const int _num__array_svmon_inh_add_I_Na;
extern double *_array_svmon_inh_add_m;
extern const int _num__array_svmon_inh_add_m;
extern int32_t *_array_svmon_inh_add_N;
extern const int _num__array_svmon_inh_add_N;
extern double *_array_svmon_inh_add_n;
extern const int _num__array_svmon_inh_add_n;
extern double *_array_svmon_inh_add_p;
extern const int _num__array_svmon_inh_add_p;
extern double *_array_svmon_inh_add_p_inf;
extern const int _num__array_svmon_inh_add_p_inf;
extern double *_array_svmon_inh_add_tau_h;
extern const int _num__array_svmon_inh_add_tau_h;
extern double *_array_svmon_inh_add_tau_m;
extern const int _num__array_svmon_inh_add_tau_m;
extern double *_array_svmon_inh_add_tau_n;
extern const int _num__array_svmon_inh_add_tau_n;
extern double *_array_svmon_inh_add_tau_p;
extern const int _num__array_svmon_inh_add_tau_p;
extern double *_array_svmon_inh_add_v;
extern const int _num__array_svmon_inh_add_v;

//////////////// dynamic arrays 2d /////////
extern DynamicArray2D<double> _dynamic_array_svmon_exc_add_a_m;
extern DynamicArray2D<double> _dynamic_array_svmon_exc_add_b_m;
extern DynamicArray2D<double> _dynamic_array_svmon_exc_add_h;
extern DynamicArray2D<double> _dynamic_array_svmon_exc_add_I_inj;
extern DynamicArray2D<double> _dynamic_array_svmon_exc_add_I_Kd;
extern DynamicArray2D<double> _dynamic_array_svmon_exc_add_I_L;
extern DynamicArray2D<double> _dynamic_array_svmon_exc_add_I_M;
extern DynamicArray2D<double> _dynamic_array_svmon_exc_add_I_max;
extern DynamicArray2D<double> _dynamic_array_svmon_exc_add_I_Na;
extern DynamicArray2D<double> _dynamic_array_svmon_exc_add_m;
extern DynamicArray2D<double> _dynamic_array_svmon_exc_add_n;
extern DynamicArray2D<double> _dynamic_array_svmon_exc_add_p;
extern DynamicArray2D<double> _dynamic_array_svmon_exc_add_p_inf;
extern DynamicArray2D<double> _dynamic_array_svmon_exc_add_tau_h;
extern DynamicArray2D<double> _dynamic_array_svmon_exc_add_tau_m;
extern DynamicArray2D<double> _dynamic_array_svmon_exc_add_tau_n;
extern DynamicArray2D<double> _dynamic_array_svmon_exc_add_tau_p;
extern DynamicArray2D<double> _dynamic_array_svmon_exc_add_v;
extern DynamicArray2D<double> _dynamic_array_svmon_inh_add_a_m;
extern DynamicArray2D<double> _dynamic_array_svmon_inh_add_b_m;
extern DynamicArray2D<double> _dynamic_array_svmon_inh_add_h;
extern DynamicArray2D<double> _dynamic_array_svmon_inh_add_I_inj;
extern DynamicArray2D<double> _dynamic_array_svmon_inh_add_I_Kd;
extern DynamicArray2D<double> _dynamic_array_svmon_inh_add_I_L;
extern DynamicArray2D<double> _dynamic_array_svmon_inh_add_I_M;
extern DynamicArray2D<double> _dynamic_array_svmon_inh_add_I_max;
extern DynamicArray2D<double> _dynamic_array_svmon_inh_add_I_Na;
extern DynamicArray2D<double> _dynamic_array_svmon_inh_add_m;
extern DynamicArray2D<double> _dynamic_array_svmon_inh_add_n;
extern DynamicArray2D<double> _dynamic_array_svmon_inh_add_p;
extern DynamicArray2D<double> _dynamic_array_svmon_inh_add_p_inf;
extern DynamicArray2D<double> _dynamic_array_svmon_inh_add_tau_h;
extern DynamicArray2D<double> _dynamic_array_svmon_inh_add_tau_m;
extern DynamicArray2D<double> _dynamic_array_svmon_inh_add_tau_n;
extern DynamicArray2D<double> _dynamic_array_svmon_inh_add_tau_p;
extern DynamicArray2D<double> _dynamic_array_svmon_inh_add_v;

/////////////// static arrays /////////////
extern int32_t *_static_array__index__array_NGexc_x;
extern const int _num__static_array__index__array_NGexc_x;
extern int32_t *_static_array__index__array_NGexc_y;
extern const int _num__static_array__index__array_NGexc_y;
extern int32_t *_static_array__index__array_NGexc_z;
extern const int _num__static_array__index__array_NGexc_z;
extern int32_t *_static_array__index__array_NGinh_x;
extern const int _num__static_array__index__array_NGinh_x;
extern int32_t *_static_array__index__array_NGinh_y;
extern const int _num__static_array__index__array_NGinh_y;
extern int32_t *_static_array__index__array_NGinh_z;
extern const int _num__static_array__index__array_NGinh_z;
extern double *_static_array__value__array_NGexc_x;
extern const int _num__static_array__value__array_NGexc_x;
extern double *_static_array__value__array_NGexc_y;
extern const int _num__static_array__value__array_NGexc_y;
extern double *_static_array__value__array_NGexc_z;
extern const int _num__static_array__value__array_NGexc_z;
extern double *_static_array__value__array_NGinh_x;
extern const int _num__static_array__value__array_NGinh_x;
extern double *_static_array__value__array_NGinh_y;
extern const int _num__static_array__value__array_NGinh_y;
extern double *_static_array__value__array_NGinh_z;
extern const int _num__static_array__value__array_NGinh_z;
extern double *_timedarray_3_values;
extern const int _num__timedarray_3_values;
extern double *_timedarray_4_values;
extern const int _num__timedarray_4_values;

//////////////// synapses /////////////////

// Profiling information for each code object
}

void _init_arrays();
void _load_arrays();
void _write_arrays();
void _dealloc_arrays();

#endif


