

#include "objects.h"
#include "synapses_classes.h"
#include "brianlib/clocks.h"
#include "brianlib/dynamic_array.h"
#include "brianlib/stdint_compat.h"
#include "network.h"
#include "randomkit.h"
#include<vector>
#include<iostream>
#include<fstream>
#include<map>
#include<tuple>
#include<cstdlib>
#include<string>

namespace brian {

std::string results_dir = "results/";  // can be overwritten by --results_dir command line arg
std::vector< rk_state* > _mersenne_twister_states;

//////////////// networks /////////////////
Network network_1;

void set_variable_from_value(std::string varname, char* var_pointer, size_t size, char value) {
    #ifdef DEBUG
    std::cout << "Setting '" << varname << "' to " << (value == 1 ? "True" : "False") << std::endl;
    #endif
    std::fill(var_pointer, var_pointer+size, value);
}

template<class T> void set_variable_from_value(std::string varname, T* var_pointer, size_t size, T value) {
    #ifdef DEBUG
    std::cout << "Setting '" << varname << "' to " << value << std::endl;
    #endif
    std::fill(var_pointer, var_pointer+size, value);
}

template<class T> void set_variable_from_file(std::string varname, T* var_pointer, size_t data_size, std::string filename) {
    ifstream f;
    streampos size;
    #ifdef DEBUG
    std::cout << "Setting '" << varname << "' from file '" << filename << "'" << std::endl;
    #endif
    f.open(filename, ios::in | ios::binary | ios::ate);
    size = f.tellg();
    if (size != data_size) {
        std::cerr << "Error reading '" << filename << "': file size " << size << " does not match expected size " << data_size << std::endl;
        return;
    }
    f.seekg(0, ios::beg);
    if (f.is_open())
        f.read(reinterpret_cast<char *>(var_pointer), data_size);
    else
        std::cerr << "Could not read '" << filename << "'" << std::endl;
    if (f.fail())
        std::cerr << "Error reading '" << filename << "'" << std::endl;
}

//////////////// set arrays by name ///////
void set_variable_by_name(std::string name, std::string s_value) {
	size_t var_size;
	size_t data_size;
	std::for_each(s_value.begin(), s_value.end(), [](char& c) // modify in-place
    {
        c = std::tolower(static_cast<unsigned char>(c));
    });
    if (s_value == "true")
        s_value = "1";
    else if (s_value == "false")
        s_value = "0";
	// non-dynamic arrays
    if (name == "NGexc._spikespace") {
        var_size = 2;
        data_size = 2*sizeof(int32_t);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<int32_t>(name, _array_NGexc__spikespace, var_size, (int32_t)atoi(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_NGexc__spikespace, data_size, s_value);
        }
        return;
    }
    if (name == "NGexc.h") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_NGexc_h, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_NGexc_h, data_size, s_value);
        }
        return;
    }
    if (name == "NGexc.I_max") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_NGexc_I_max, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_NGexc_I_max, data_size, s_value);
        }
        return;
    }
    if (name == "NGexc.lastspike") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_NGexc_lastspike, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_NGexc_lastspike, data_size, s_value);
        }
        return;
    }
    if (name == "NGexc.m") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_NGexc_m, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_NGexc_m, data_size, s_value);
        }
        return;
    }
    if (name == "NGexc.mask_noise") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_NGexc_mask_noise, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_NGexc_mask_noise, data_size, s_value);
        }
        return;
    }
    if (name == "NGexc.n") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_NGexc_n, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_NGexc_n, data_size, s_value);
        }
        return;
    }
    if (name == "NGexc.not_refractory") {
        var_size = 1;
        data_size = 1*sizeof(char);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value(name, _array_NGexc_not_refractory, var_size, (char)atoi(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_NGexc_not_refractory, data_size, s_value);
        }
        return;
    }
    if (name == "NGexc.p") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_NGexc_p, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_NGexc_p, data_size, s_value);
        }
        return;
    }
    if (name == "NGexc.v") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_NGexc_v, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_NGexc_v, data_size, s_value);
        }
        return;
    }
    if (name == "NGexc.x") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_NGexc_x, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_NGexc_x, data_size, s_value);
        }
        return;
    }
    if (name == "NGexc.y") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_NGexc_y, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_NGexc_y, data_size, s_value);
        }
        return;
    }
    if (name == "NGexc.z") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_NGexc_z, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_NGexc_z, data_size, s_value);
        }
        return;
    }
    if (name == "NGinh._spikespace") {
        var_size = 2;
        data_size = 2*sizeof(int32_t);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<int32_t>(name, _array_NGinh__spikespace, var_size, (int32_t)atoi(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_NGinh__spikespace, data_size, s_value);
        }
        return;
    }
    if (name == "NGinh.h") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_NGinh_h, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_NGinh_h, data_size, s_value);
        }
        return;
    }
    if (name == "NGinh.I_max") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_NGinh_I_max, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_NGinh_I_max, data_size, s_value);
        }
        return;
    }
    if (name == "NGinh.lastspike") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_NGinh_lastspike, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_NGinh_lastspike, data_size, s_value);
        }
        return;
    }
    if (name == "NGinh.m") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_NGinh_m, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_NGinh_m, data_size, s_value);
        }
        return;
    }
    if (name == "NGinh.mask_noise") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_NGinh_mask_noise, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_NGinh_mask_noise, data_size, s_value);
        }
        return;
    }
    if (name == "NGinh.n") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_NGinh_n, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_NGinh_n, data_size, s_value);
        }
        return;
    }
    if (name == "NGinh.not_refractory") {
        var_size = 1;
        data_size = 1*sizeof(char);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value(name, _array_NGinh_not_refractory, var_size, (char)atoi(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_NGinh_not_refractory, data_size, s_value);
        }
        return;
    }
    if (name == "NGinh.p") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_NGinh_p, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_NGinh_p, data_size, s_value);
        }
        return;
    }
    if (name == "NGinh.v") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_NGinh_v, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_NGinh_v, data_size, s_value);
        }
        return;
    }
    if (name == "NGinh.x") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_NGinh_x, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_NGinh_x, data_size, s_value);
        }
        return;
    }
    if (name == "NGinh.y") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_NGinh_y, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_NGinh_y, data_size, s_value);
        }
        return;
    }
    if (name == "NGinh.z") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_NGinh_z, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_NGinh_z, data_size, s_value);
        }
        return;
    }
    // dynamic arrays (1d)
    if (name == "_timedarray_3.values") {
        var_size = 10;
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _timedarray_3_values, var_size, (double)atof(s_value.c_str()));


        } else {
            // set from file
            set_variable_from_file(name, _timedarray_3_values, data_size, s_value);
        }
        return;
    }
    if (name == "_timedarray_4.values") {
        var_size = 10;
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _timedarray_4_values, var_size, (double)atof(s_value.c_str()));


        } else {
            // set from file
            set_variable_from_file(name, _timedarray_4_values, data_size, s_value);
        }
        return;
    }
    std::cerr << "Cannot set unknown variable '" << name << "'." << std::endl;
    exit(1);
}
//////////////// arrays ///////////////////
double * _array_defaultclock_dt;
const int _num__array_defaultclock_dt = 1;
double * _array_defaultclock_t;
const int _num__array_defaultclock_t = 1;
int64_t * _array_defaultclock_timestep;
const int _num__array_defaultclock_timestep = 1;
int32_t * _array_NGexc__spikespace;
const int _num__array_NGexc__spikespace = 2;
double * _array_NGexc_clock_1_dt;
const int _num__array_NGexc_clock_1_dt = 1;
double * _array_NGexc_clock_1_t;
const int _num__array_NGexc_clock_1_t = 1;
int64_t * _array_NGexc_clock_1_timestep;
const int _num__array_NGexc_clock_1_timestep = 1;
double * _array_NGexc_h;
const int _num__array_NGexc_h = 1;
int32_t * _array_NGexc_i;
const int _num__array_NGexc_i = 1;
double * _array_NGexc_I_max;
const int _num__array_NGexc_I_max = 1;
double * _array_NGexc_lastspike;
const int _num__array_NGexc_lastspike = 1;
double * _array_NGexc_m;
const int _num__array_NGexc_m = 1;
double * _array_NGexc_mask_noise;
const int _num__array_NGexc_mask_noise = 1;
double * _array_NGexc_n;
const int _num__array_NGexc_n = 1;
char * _array_NGexc_not_refractory;
const int _num__array_NGexc_not_refractory = 1;
double * _array_NGexc_p;
const int _num__array_NGexc_p = 1;
int32_t * _array_NGexc_subgroup__sub_idx;
const int _num__array_NGexc_subgroup__sub_idx = 1;
int32_t * _array_NGexc_subgroup__sub_idx_1;
const int _num__array_NGexc_subgroup__sub_idx_1 = 1;
int32_t * _array_NGexc_subgroup__sub_idx_2;
const int _num__array_NGexc_subgroup__sub_idx_2 = 1;
double * _array_NGexc_v;
const int _num__array_NGexc_v = 1;
double * _array_NGexc_x;
const int _num__array_NGexc_x = 1;
double * _array_NGexc_y;
const int _num__array_NGexc_y = 1;
double * _array_NGexc_z;
const int _num__array_NGexc_z = 1;
int32_t * _array_NGinh__spikespace;
const int _num__array_NGinh__spikespace = 2;
double * _array_NGinh_clock_1_dt;
const int _num__array_NGinh_clock_1_dt = 1;
double * _array_NGinh_clock_1_t;
const int _num__array_NGinh_clock_1_t = 1;
int64_t * _array_NGinh_clock_1_timestep;
const int _num__array_NGinh_clock_1_timestep = 1;
double * _array_NGinh_h;
const int _num__array_NGinh_h = 1;
int32_t * _array_NGinh_i;
const int _num__array_NGinh_i = 1;
double * _array_NGinh_I_max;
const int _num__array_NGinh_I_max = 1;
double * _array_NGinh_lastspike;
const int _num__array_NGinh_lastspike = 1;
double * _array_NGinh_m;
const int _num__array_NGinh_m = 1;
double * _array_NGinh_mask_noise;
const int _num__array_NGinh_mask_noise = 1;
double * _array_NGinh_n;
const int _num__array_NGinh_n = 1;
char * _array_NGinh_not_refractory;
const int _num__array_NGinh_not_refractory = 1;
double * _array_NGinh_p;
const int _num__array_NGinh_p = 1;
int32_t * _array_NGinh_subgroup__sub_idx;
const int _num__array_NGinh_subgroup__sub_idx = 1;
int32_t * _array_NGinh_subgroup__sub_idx_1;
const int _num__array_NGinh_subgroup__sub_idx_1 = 1;
int32_t * _array_NGinh_subgroup__sub_idx_2;
const int _num__array_NGinh_subgroup__sub_idx_2 = 1;
double * _array_NGinh_v;
const int _num__array_NGinh_v = 1;
double * _array_NGinh_x;
const int _num__array_NGinh_x = 1;
double * _array_NGinh_y;
const int _num__array_NGinh_y = 1;
double * _array_NGinh_z;
const int _num__array_NGinh_z = 1;
int32_t * _array_spikemonitor_2__source_idx;
const int _num__array_spikemonitor_2__source_idx = 1;
int32_t * _array_spikemonitor_2_count;
const int _num__array_spikemonitor_2_count = 1;
int32_t * _array_spikemonitor_2_N;
const int _num__array_spikemonitor_2_N = 1;
int32_t * _array_spikemonitor_3__source_idx;
const int _num__array_spikemonitor_3__source_idx = 1;
int32_t * _array_spikemonitor_3_count;
const int _num__array_spikemonitor_3_count = 1;
int32_t * _array_spikemonitor_3_N;
const int _num__array_spikemonitor_3_N = 1;
int32_t * _array_svmon_exc_add__indices;
const int _num__array_svmon_exc_add__indices = 1;
double * _array_svmon_exc_add_a_m;
const int _num__array_svmon_exc_add_a_m = (0, 1);
double * _array_svmon_exc_add_b_m;
const int _num__array_svmon_exc_add_b_m = (0, 1);
double * _array_svmon_exc_add_clock_1_dt;
const int _num__array_svmon_exc_add_clock_1_dt = 1;
double * _array_svmon_exc_add_clock_1_t;
const int _num__array_svmon_exc_add_clock_1_t = 1;
int64_t * _array_svmon_exc_add_clock_1_timestep;
const int _num__array_svmon_exc_add_clock_1_timestep = 1;
double * _array_svmon_exc_add_h;
const int _num__array_svmon_exc_add_h = (0, 1);
double * _array_svmon_exc_add_I_inj;
const int _num__array_svmon_exc_add_I_inj = (0, 1);
double * _array_svmon_exc_add_I_Kd;
const int _num__array_svmon_exc_add_I_Kd = (0, 1);
double * _array_svmon_exc_add_I_L;
const int _num__array_svmon_exc_add_I_L = (0, 1);
double * _array_svmon_exc_add_I_M;
const int _num__array_svmon_exc_add_I_M = (0, 1);
double * _array_svmon_exc_add_I_max;
const int _num__array_svmon_exc_add_I_max = (0, 1);
double * _array_svmon_exc_add_I_Na;
const int _num__array_svmon_exc_add_I_Na = (0, 1);
double * _array_svmon_exc_add_m;
const int _num__array_svmon_exc_add_m = (0, 1);
int32_t * _array_svmon_exc_add_N;
const int _num__array_svmon_exc_add_N = 1;
double * _array_svmon_exc_add_n;
const int _num__array_svmon_exc_add_n = (0, 1);
double * _array_svmon_exc_add_p;
const int _num__array_svmon_exc_add_p = (0, 1);
double * _array_svmon_exc_add_p_inf;
const int _num__array_svmon_exc_add_p_inf = (0, 1);
double * _array_svmon_exc_add_tau_h;
const int _num__array_svmon_exc_add_tau_h = (0, 1);
double * _array_svmon_exc_add_tau_m;
const int _num__array_svmon_exc_add_tau_m = (0, 1);
double * _array_svmon_exc_add_tau_n;
const int _num__array_svmon_exc_add_tau_n = (0, 1);
double * _array_svmon_exc_add_tau_p;
const int _num__array_svmon_exc_add_tau_p = (0, 1);
double * _array_svmon_exc_add_v;
const int _num__array_svmon_exc_add_v = (0, 1);
int32_t * _array_svmon_inh_add__indices;
const int _num__array_svmon_inh_add__indices = 1;
double * _array_svmon_inh_add_a_m;
const int _num__array_svmon_inh_add_a_m = (0, 1);
double * _array_svmon_inh_add_b_m;
const int _num__array_svmon_inh_add_b_m = (0, 1);
double * _array_svmon_inh_add_clock_1_dt;
const int _num__array_svmon_inh_add_clock_1_dt = 1;
double * _array_svmon_inh_add_clock_1_t;
const int _num__array_svmon_inh_add_clock_1_t = 1;
int64_t * _array_svmon_inh_add_clock_1_timestep;
const int _num__array_svmon_inh_add_clock_1_timestep = 1;
double * _array_svmon_inh_add_h;
const int _num__array_svmon_inh_add_h = (0, 1);
double * _array_svmon_inh_add_I_inj;
const int _num__array_svmon_inh_add_I_inj = (0, 1);
double * _array_svmon_inh_add_I_Kd;
const int _num__array_svmon_inh_add_I_Kd = (0, 1);
double * _array_svmon_inh_add_I_L;
const int _num__array_svmon_inh_add_I_L = (0, 1);
double * _array_svmon_inh_add_I_M;
const int _num__array_svmon_inh_add_I_M = (0, 1);
double * _array_svmon_inh_add_I_max;
const int _num__array_svmon_inh_add_I_max = (0, 1);
double * _array_svmon_inh_add_I_Na;
const int _num__array_svmon_inh_add_I_Na = (0, 1);
double * _array_svmon_inh_add_m;
const int _num__array_svmon_inh_add_m = (0, 1);
int32_t * _array_svmon_inh_add_N;
const int _num__array_svmon_inh_add_N = 1;
double * _array_svmon_inh_add_n;
const int _num__array_svmon_inh_add_n = (0, 1);
double * _array_svmon_inh_add_p;
const int _num__array_svmon_inh_add_p = (0, 1);
double * _array_svmon_inh_add_p_inf;
const int _num__array_svmon_inh_add_p_inf = (0, 1);
double * _array_svmon_inh_add_tau_h;
const int _num__array_svmon_inh_add_tau_h = (0, 1);
double * _array_svmon_inh_add_tau_m;
const int _num__array_svmon_inh_add_tau_m = (0, 1);
double * _array_svmon_inh_add_tau_n;
const int _num__array_svmon_inh_add_tau_n = (0, 1);
double * _array_svmon_inh_add_tau_p;
const int _num__array_svmon_inh_add_tau_p = (0, 1);
double * _array_svmon_inh_add_v;
const int _num__array_svmon_inh_add_v = (0, 1);

//////////////// dynamic arrays 1d /////////
std::vector<int32_t> _dynamic_array_spikemonitor_2_i;
std::vector<double> _dynamic_array_spikemonitor_2_t;
std::vector<int32_t> _dynamic_array_spikemonitor_3_i;
std::vector<double> _dynamic_array_spikemonitor_3_t;
std::vector<double> _dynamic_array_svmon_exc_add_t;
std::vector<double> _dynamic_array_svmon_inh_add_t;

//////////////// dynamic arrays 2d /////////
DynamicArray2D<double> _dynamic_array_svmon_exc_add_a_m;
DynamicArray2D<double> _dynamic_array_svmon_exc_add_b_m;
DynamicArray2D<double> _dynamic_array_svmon_exc_add_h;
DynamicArray2D<double> _dynamic_array_svmon_exc_add_I_inj;
DynamicArray2D<double> _dynamic_array_svmon_exc_add_I_Kd;
DynamicArray2D<double> _dynamic_array_svmon_exc_add_I_L;
DynamicArray2D<double> _dynamic_array_svmon_exc_add_I_M;
DynamicArray2D<double> _dynamic_array_svmon_exc_add_I_max;
DynamicArray2D<double> _dynamic_array_svmon_exc_add_I_Na;
DynamicArray2D<double> _dynamic_array_svmon_exc_add_m;
DynamicArray2D<double> _dynamic_array_svmon_exc_add_n;
DynamicArray2D<double> _dynamic_array_svmon_exc_add_p;
DynamicArray2D<double> _dynamic_array_svmon_exc_add_p_inf;
DynamicArray2D<double> _dynamic_array_svmon_exc_add_tau_h;
DynamicArray2D<double> _dynamic_array_svmon_exc_add_tau_m;
DynamicArray2D<double> _dynamic_array_svmon_exc_add_tau_n;
DynamicArray2D<double> _dynamic_array_svmon_exc_add_tau_p;
DynamicArray2D<double> _dynamic_array_svmon_exc_add_v;
DynamicArray2D<double> _dynamic_array_svmon_inh_add_a_m;
DynamicArray2D<double> _dynamic_array_svmon_inh_add_b_m;
DynamicArray2D<double> _dynamic_array_svmon_inh_add_h;
DynamicArray2D<double> _dynamic_array_svmon_inh_add_I_inj;
DynamicArray2D<double> _dynamic_array_svmon_inh_add_I_Kd;
DynamicArray2D<double> _dynamic_array_svmon_inh_add_I_L;
DynamicArray2D<double> _dynamic_array_svmon_inh_add_I_M;
DynamicArray2D<double> _dynamic_array_svmon_inh_add_I_max;
DynamicArray2D<double> _dynamic_array_svmon_inh_add_I_Na;
DynamicArray2D<double> _dynamic_array_svmon_inh_add_m;
DynamicArray2D<double> _dynamic_array_svmon_inh_add_n;
DynamicArray2D<double> _dynamic_array_svmon_inh_add_p;
DynamicArray2D<double> _dynamic_array_svmon_inh_add_p_inf;
DynamicArray2D<double> _dynamic_array_svmon_inh_add_tau_h;
DynamicArray2D<double> _dynamic_array_svmon_inh_add_tau_m;
DynamicArray2D<double> _dynamic_array_svmon_inh_add_tau_n;
DynamicArray2D<double> _dynamic_array_svmon_inh_add_tau_p;
DynamicArray2D<double> _dynamic_array_svmon_inh_add_v;

/////////////// static arrays /////////////
int32_t * _static_array__index__array_NGexc_x;
const int _num__static_array__index__array_NGexc_x = 1;
int32_t * _static_array__index__array_NGexc_y;
const int _num__static_array__index__array_NGexc_y = 1;
int32_t * _static_array__index__array_NGexc_z;
const int _num__static_array__index__array_NGexc_z = 1;
int32_t * _static_array__index__array_NGinh_x;
const int _num__static_array__index__array_NGinh_x = 1;
int32_t * _static_array__index__array_NGinh_y;
const int _num__static_array__index__array_NGinh_y = 1;
int32_t * _static_array__index__array_NGinh_z;
const int _num__static_array__index__array_NGinh_z = 1;
double * _static_array__value__array_NGexc_x;
const int _num__static_array__value__array_NGexc_x = 1;
double * _static_array__value__array_NGexc_y;
const int _num__static_array__value__array_NGexc_y = 1;
double * _static_array__value__array_NGexc_z;
const int _num__static_array__value__array_NGexc_z = 1;
double * _static_array__value__array_NGinh_x;
const int _num__static_array__value__array_NGinh_x = 1;
double * _static_array__value__array_NGinh_y;
const int _num__static_array__value__array_NGinh_y = 1;
double * _static_array__value__array_NGinh_z;
const int _num__static_array__value__array_NGinh_z = 1;
double * _timedarray_3_values;
const int _num__timedarray_3_values = 10;
double * _timedarray_4_values;
const int _num__timedarray_4_values = 10;

//////////////// synapses /////////////////

//////////////// clocks ///////////////////
Clock NGexc_clock_1;  // attributes will be set in run.cpp
Clock NGinh_clock_1;  // attributes will be set in run.cpp
Clock svmon_exc_add_clock_1;  // attributes will be set in run.cpp
Clock svmon_inh_add_clock_1;  // attributes will be set in run.cpp

// Profiling information for each code object
}

void _init_arrays()
{
	using namespace brian;

    // Arrays initialized to 0
	_array_defaultclock_dt = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_defaultclock_dt[i] = 0;

	_array_defaultclock_t = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_defaultclock_t[i] = 0;

	_array_defaultclock_timestep = new int64_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_defaultclock_timestep[i] = 0;

	_array_NGexc__spikespace = new int32_t[2];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<2; i++) _array_NGexc__spikespace[i] = 0;

	_array_NGexc_clock_1_dt = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGexc_clock_1_dt[i] = 0;

	_array_NGexc_clock_1_t = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGexc_clock_1_t[i] = 0;

	_array_NGexc_clock_1_timestep = new int64_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGexc_clock_1_timestep[i] = 0;

	_array_NGexc_h = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGexc_h[i] = 0;

	_array_NGexc_i = new int32_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGexc_i[i] = 0;

	_array_NGexc_I_max = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGexc_I_max[i] = 0;

	_array_NGexc_lastspike = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGexc_lastspike[i] = 0;

	_array_NGexc_m = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGexc_m[i] = 0;

	_array_NGexc_mask_noise = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGexc_mask_noise[i] = 0;

	_array_NGexc_n = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGexc_n[i] = 0;

	_array_NGexc_not_refractory = new char[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGexc_not_refractory[i] = 0;

	_array_NGexc_p = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGexc_p[i] = 0;

	_array_NGexc_subgroup__sub_idx = new int32_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGexc_subgroup__sub_idx[i] = 0;

	_array_NGexc_subgroup__sub_idx_1 = new int32_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGexc_subgroup__sub_idx_1[i] = 0;

	_array_NGexc_subgroup__sub_idx_2 = new int32_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGexc_subgroup__sub_idx_2[i] = 0;

	_array_NGexc_v = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGexc_v[i] = 0;

	_array_NGexc_x = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGexc_x[i] = 0;

	_array_NGexc_y = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGexc_y[i] = 0;

	_array_NGexc_z = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGexc_z[i] = 0;

	_array_NGinh__spikespace = new int32_t[2];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<2; i++) _array_NGinh__spikespace[i] = 0;

	_array_NGinh_clock_1_dt = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGinh_clock_1_dt[i] = 0;

	_array_NGinh_clock_1_t = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGinh_clock_1_t[i] = 0;

	_array_NGinh_clock_1_timestep = new int64_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGinh_clock_1_timestep[i] = 0;

	_array_NGinh_h = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGinh_h[i] = 0;

	_array_NGinh_i = new int32_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGinh_i[i] = 0;

	_array_NGinh_I_max = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGinh_I_max[i] = 0;

	_array_NGinh_lastspike = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGinh_lastspike[i] = 0;

	_array_NGinh_m = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGinh_m[i] = 0;

	_array_NGinh_mask_noise = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGinh_mask_noise[i] = 0;

	_array_NGinh_n = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGinh_n[i] = 0;

	_array_NGinh_not_refractory = new char[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGinh_not_refractory[i] = 0;

	_array_NGinh_p = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGinh_p[i] = 0;

	_array_NGinh_subgroup__sub_idx = new int32_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGinh_subgroup__sub_idx[i] = 0;

	_array_NGinh_subgroup__sub_idx_1 = new int32_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGinh_subgroup__sub_idx_1[i] = 0;

	_array_NGinh_subgroup__sub_idx_2 = new int32_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGinh_subgroup__sub_idx_2[i] = 0;

	_array_NGinh_v = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGinh_v[i] = 0;

	_array_NGinh_x = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGinh_x[i] = 0;

	_array_NGinh_y = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGinh_y[i] = 0;

	_array_NGinh_z = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGinh_z[i] = 0;

	_array_spikemonitor_2__source_idx = new int32_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_spikemonitor_2__source_idx[i] = 0;

	_array_spikemonitor_2_count = new int32_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_spikemonitor_2_count[i] = 0;

	_array_spikemonitor_2_N = new int32_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_spikemonitor_2_N[i] = 0;

	_array_spikemonitor_3__source_idx = new int32_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_spikemonitor_3__source_idx[i] = 0;

	_array_spikemonitor_3_count = new int32_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_spikemonitor_3_count[i] = 0;

	_array_spikemonitor_3_N = new int32_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_spikemonitor_3_N[i] = 0;

	_array_svmon_exc_add__indices = new int32_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_svmon_exc_add__indices[i] = 0;

	_array_svmon_exc_add_clock_1_dt = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_svmon_exc_add_clock_1_dt[i] = 0;

	_array_svmon_exc_add_clock_1_t = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_svmon_exc_add_clock_1_t[i] = 0;

	_array_svmon_exc_add_clock_1_timestep = new int64_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_svmon_exc_add_clock_1_timestep[i] = 0;

	_array_svmon_exc_add_N = new int32_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_svmon_exc_add_N[i] = 0;

	_array_svmon_inh_add__indices = new int32_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_svmon_inh_add__indices[i] = 0;

	_array_svmon_inh_add_clock_1_dt = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_svmon_inh_add_clock_1_dt[i] = 0;

	_array_svmon_inh_add_clock_1_t = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_svmon_inh_add_clock_1_t[i] = 0;

	_array_svmon_inh_add_clock_1_timestep = new int64_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_svmon_inh_add_clock_1_timestep[i] = 0;

	_array_svmon_inh_add_N = new int32_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_svmon_inh_add_N[i] = 0;


	// Arrays initialized to an "arange"
	_array_NGexc_i = new int32_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGexc_i[i] = 0 + i;

	_array_NGexc_subgroup__sub_idx = new int32_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGexc_subgroup__sub_idx[i] = 0 + i;

	_array_NGexc_subgroup__sub_idx_1 = new int32_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGexc_subgroup__sub_idx_1[i] = 0 + i;

	_array_NGexc_subgroup__sub_idx_2 = new int32_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGexc_subgroup__sub_idx_2[i] = 0 + i;

	_array_NGinh_i = new int32_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGinh_i[i] = 0 + i;

	_array_NGinh_subgroup__sub_idx = new int32_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGinh_subgroup__sub_idx[i] = 0 + i;

	_array_NGinh_subgroup__sub_idx_1 = new int32_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGinh_subgroup__sub_idx_1[i] = 0 + i;

	_array_NGinh_subgroup__sub_idx_2 = new int32_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_NGinh_subgroup__sub_idx_2[i] = 0 + i;

	_array_spikemonitor_2__source_idx = new int32_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_spikemonitor_2__source_idx[i] = 0 + i;

	_array_spikemonitor_3__source_idx = new int32_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_spikemonitor_3__source_idx[i] = 0 + i;


	// static arrays
	_static_array__index__array_NGexc_x = new int32_t[1];
	_static_array__index__array_NGexc_y = new int32_t[1];
	_static_array__index__array_NGexc_z = new int32_t[1];
	_static_array__index__array_NGinh_x = new int32_t[1];
	_static_array__index__array_NGinh_y = new int32_t[1];
	_static_array__index__array_NGinh_z = new int32_t[1];
	_static_array__value__array_NGexc_x = new double[1];
	_static_array__value__array_NGexc_y = new double[1];
	_static_array__value__array_NGexc_z = new double[1];
	_static_array__value__array_NGinh_x = new double[1];
	_static_array__value__array_NGinh_y = new double[1];
	_static_array__value__array_NGinh_z = new double[1];
	_timedarray_3_values = new double[10];
	_timedarray_4_values = new double[10];

	// Random number generator states
	for (int i=0; i<1; i++)
	    _mersenne_twister_states.push_back(new rk_state());
}

void _load_arrays()
{
	using namespace brian;

	ifstream f_static_array__index__array_NGexc_x;
	f_static_array__index__array_NGexc_x.open("static_arrays/_static_array__index__array_NGexc_x", ios::in | ios::binary);
	if(f_static_array__index__array_NGexc_x.is_open())
	{
		f_static_array__index__array_NGexc_x.read(reinterpret_cast<char*>(_static_array__index__array_NGexc_x), 1*sizeof(int32_t));
	} else
	{
		std::cout << "Error opening static array _static_array__index__array_NGexc_x." << endl;
	}
	ifstream f_static_array__index__array_NGexc_y;
	f_static_array__index__array_NGexc_y.open("static_arrays/_static_array__index__array_NGexc_y", ios::in | ios::binary);
	if(f_static_array__index__array_NGexc_y.is_open())
	{
		f_static_array__index__array_NGexc_y.read(reinterpret_cast<char*>(_static_array__index__array_NGexc_y), 1*sizeof(int32_t));
	} else
	{
		std::cout << "Error opening static array _static_array__index__array_NGexc_y." << endl;
	}
	ifstream f_static_array__index__array_NGexc_z;
	f_static_array__index__array_NGexc_z.open("static_arrays/_static_array__index__array_NGexc_z", ios::in | ios::binary);
	if(f_static_array__index__array_NGexc_z.is_open())
	{
		f_static_array__index__array_NGexc_z.read(reinterpret_cast<char*>(_static_array__index__array_NGexc_z), 1*sizeof(int32_t));
	} else
	{
		std::cout << "Error opening static array _static_array__index__array_NGexc_z." << endl;
	}
	ifstream f_static_array__index__array_NGinh_x;
	f_static_array__index__array_NGinh_x.open("static_arrays/_static_array__index__array_NGinh_x", ios::in | ios::binary);
	if(f_static_array__index__array_NGinh_x.is_open())
	{
		f_static_array__index__array_NGinh_x.read(reinterpret_cast<char*>(_static_array__index__array_NGinh_x), 1*sizeof(int32_t));
	} else
	{
		std::cout << "Error opening static array _static_array__index__array_NGinh_x." << endl;
	}
	ifstream f_static_array__index__array_NGinh_y;
	f_static_array__index__array_NGinh_y.open("static_arrays/_static_array__index__array_NGinh_y", ios::in | ios::binary);
	if(f_static_array__index__array_NGinh_y.is_open())
	{
		f_static_array__index__array_NGinh_y.read(reinterpret_cast<char*>(_static_array__index__array_NGinh_y), 1*sizeof(int32_t));
	} else
	{
		std::cout << "Error opening static array _static_array__index__array_NGinh_y." << endl;
	}
	ifstream f_static_array__index__array_NGinh_z;
	f_static_array__index__array_NGinh_z.open("static_arrays/_static_array__index__array_NGinh_z", ios::in | ios::binary);
	if(f_static_array__index__array_NGinh_z.is_open())
	{
		f_static_array__index__array_NGinh_z.read(reinterpret_cast<char*>(_static_array__index__array_NGinh_z), 1*sizeof(int32_t));
	} else
	{
		std::cout << "Error opening static array _static_array__index__array_NGinh_z." << endl;
	}
	ifstream f_static_array__value__array_NGexc_x;
	f_static_array__value__array_NGexc_x.open("static_arrays/_static_array__value__array_NGexc_x", ios::in | ios::binary);
	if(f_static_array__value__array_NGexc_x.is_open())
	{
		f_static_array__value__array_NGexc_x.read(reinterpret_cast<char*>(_static_array__value__array_NGexc_x), 1*sizeof(double));
	} else
	{
		std::cout << "Error opening static array _static_array__value__array_NGexc_x." << endl;
	}
	ifstream f_static_array__value__array_NGexc_y;
	f_static_array__value__array_NGexc_y.open("static_arrays/_static_array__value__array_NGexc_y", ios::in | ios::binary);
	if(f_static_array__value__array_NGexc_y.is_open())
	{
		f_static_array__value__array_NGexc_y.read(reinterpret_cast<char*>(_static_array__value__array_NGexc_y), 1*sizeof(double));
	} else
	{
		std::cout << "Error opening static array _static_array__value__array_NGexc_y." << endl;
	}
	ifstream f_static_array__value__array_NGexc_z;
	f_static_array__value__array_NGexc_z.open("static_arrays/_static_array__value__array_NGexc_z", ios::in | ios::binary);
	if(f_static_array__value__array_NGexc_z.is_open())
	{
		f_static_array__value__array_NGexc_z.read(reinterpret_cast<char*>(_static_array__value__array_NGexc_z), 1*sizeof(double));
	} else
	{
		std::cout << "Error opening static array _static_array__value__array_NGexc_z." << endl;
	}
	ifstream f_static_array__value__array_NGinh_x;
	f_static_array__value__array_NGinh_x.open("static_arrays/_static_array__value__array_NGinh_x", ios::in | ios::binary);
	if(f_static_array__value__array_NGinh_x.is_open())
	{
		f_static_array__value__array_NGinh_x.read(reinterpret_cast<char*>(_static_array__value__array_NGinh_x), 1*sizeof(double));
	} else
	{
		std::cout << "Error opening static array _static_array__value__array_NGinh_x." << endl;
	}
	ifstream f_static_array__value__array_NGinh_y;
	f_static_array__value__array_NGinh_y.open("static_arrays/_static_array__value__array_NGinh_y", ios::in | ios::binary);
	if(f_static_array__value__array_NGinh_y.is_open())
	{
		f_static_array__value__array_NGinh_y.read(reinterpret_cast<char*>(_static_array__value__array_NGinh_y), 1*sizeof(double));
	} else
	{
		std::cout << "Error opening static array _static_array__value__array_NGinh_y." << endl;
	}
	ifstream f_static_array__value__array_NGinh_z;
	f_static_array__value__array_NGinh_z.open("static_arrays/_static_array__value__array_NGinh_z", ios::in | ios::binary);
	if(f_static_array__value__array_NGinh_z.is_open())
	{
		f_static_array__value__array_NGinh_z.read(reinterpret_cast<char*>(_static_array__value__array_NGinh_z), 1*sizeof(double));
	} else
	{
		std::cout << "Error opening static array _static_array__value__array_NGinh_z." << endl;
	}
	ifstream f_timedarray_3_values;
	f_timedarray_3_values.open("static_arrays/_timedarray_3_values", ios::in | ios::binary);
	if(f_timedarray_3_values.is_open())
	{
		f_timedarray_3_values.read(reinterpret_cast<char*>(_timedarray_3_values), 10*sizeof(double));
	} else
	{
		std::cout << "Error opening static array _timedarray_3_values." << endl;
	}
	ifstream f_timedarray_4_values;
	f_timedarray_4_values.open("static_arrays/_timedarray_4_values", ios::in | ios::binary);
	if(f_timedarray_4_values.is_open())
	{
		f_timedarray_4_values.read(reinterpret_cast<char*>(_timedarray_4_values), 10*sizeof(double));
	} else
	{
		std::cout << "Error opening static array _timedarray_4_values." << endl;
	}
}

void _write_arrays()
{
	using namespace brian;

	ofstream outfile__array_defaultclock_dt;
	outfile__array_defaultclock_dt.open(results_dir + "_array_defaultclock_dt_1978099143", ios::binary | ios::out);
	if(outfile__array_defaultclock_dt.is_open())
	{
		outfile__array_defaultclock_dt.write(reinterpret_cast<char*>(_array_defaultclock_dt), 1*sizeof(_array_defaultclock_dt[0]));
		outfile__array_defaultclock_dt.close();
	} else
	{
		std::cout << "Error writing output file for _array_defaultclock_dt." << endl;
	}
	ofstream outfile__array_defaultclock_t;
	outfile__array_defaultclock_t.open(results_dir + "_array_defaultclock_t_2669362164", ios::binary | ios::out);
	if(outfile__array_defaultclock_t.is_open())
	{
		outfile__array_defaultclock_t.write(reinterpret_cast<char*>(_array_defaultclock_t), 1*sizeof(_array_defaultclock_t[0]));
		outfile__array_defaultclock_t.close();
	} else
	{
		std::cout << "Error writing output file for _array_defaultclock_t." << endl;
	}
	ofstream outfile__array_defaultclock_timestep;
	outfile__array_defaultclock_timestep.open(results_dir + "_array_defaultclock_timestep_144223508", ios::binary | ios::out);
	if(outfile__array_defaultclock_timestep.is_open())
	{
		outfile__array_defaultclock_timestep.write(reinterpret_cast<char*>(_array_defaultclock_timestep), 1*sizeof(_array_defaultclock_timestep[0]));
		outfile__array_defaultclock_timestep.close();
	} else
	{
		std::cout << "Error writing output file for _array_defaultclock_timestep." << endl;
	}
	ofstream outfile__array_NGexc__spikespace;
	outfile__array_NGexc__spikespace.open(results_dir + "_array_NGexc__spikespace_1873077302", ios::binary | ios::out);
	if(outfile__array_NGexc__spikespace.is_open())
	{
		outfile__array_NGexc__spikespace.write(reinterpret_cast<char*>(_array_NGexc__spikespace), 2*sizeof(_array_NGexc__spikespace[0]));
		outfile__array_NGexc__spikespace.close();
	} else
	{
		std::cout << "Error writing output file for _array_NGexc__spikespace." << endl;
	}
	ofstream outfile__array_NGexc_clock_1_dt;
	outfile__array_NGexc_clock_1_dt.open(results_dir + "_array_NGexc_clock_1_dt_211645572", ios::binary | ios::out);
	if(outfile__array_NGexc_clock_1_dt.is_open())
	{
		outfile__array_NGexc_clock_1_dt.write(reinterpret_cast<char*>(_array_NGexc_clock_1_dt), 1*sizeof(_array_NGexc_clock_1_dt[0]));
		outfile__array_NGexc_clock_1_dt.close();
	} else
	{
		std::cout << "Error writing output file for _array_NGexc_clock_1_dt." << endl;
	}
	ofstream outfile__array_NGexc_clock_1_t;
	outfile__array_NGexc_clock_1_t.open(results_dir + "_array_NGexc_clock_1_t_968023293", ios::binary | ios::out);
	if(outfile__array_NGexc_clock_1_t.is_open())
	{
		outfile__array_NGexc_clock_1_t.write(reinterpret_cast<char*>(_array_NGexc_clock_1_t), 1*sizeof(_array_NGexc_clock_1_t[0]));
		outfile__array_NGexc_clock_1_t.close();
	} else
	{
		std::cout << "Error writing output file for _array_NGexc_clock_1_t." << endl;
	}
	ofstream outfile__array_NGexc_clock_1_timestep;
	outfile__array_NGexc_clock_1_timestep.open(results_dir + "_array_NGexc_clock_1_timestep_2205579459", ios::binary | ios::out);
	if(outfile__array_NGexc_clock_1_timestep.is_open())
	{
		outfile__array_NGexc_clock_1_timestep.write(reinterpret_cast<char*>(_array_NGexc_clock_1_timestep), 1*sizeof(_array_NGexc_clock_1_timestep[0]));
		outfile__array_NGexc_clock_1_timestep.close();
	} else
	{
		std::cout << "Error writing output file for _array_NGexc_clock_1_timestep." << endl;
	}
	ofstream outfile__array_NGexc_h;
	outfile__array_NGexc_h.open(results_dir + "_array_NGexc_h_1052152754", ios::binary | ios::out);
	if(outfile__array_NGexc_h.is_open())
	{
		outfile__array_NGexc_h.write(reinterpret_cast<char*>(_array_NGexc_h), 1*sizeof(_array_NGexc_h[0]));
		outfile__array_NGexc_h.close();
	} else
	{
		std::cout << "Error writing output file for _array_NGexc_h." << endl;
	}
	ofstream outfile__array_NGexc_i;
	outfile__array_NGexc_i.open(results_dir + "_array_NGexc_i_1236378404", ios::binary | ios::out);
	if(outfile__array_NGexc_i.is_open())
	{
		outfile__array_NGexc_i.write(reinterpret_cast<char*>(_array_NGexc_i), 1*sizeof(_array_NGexc_i[0]));
		outfile__array_NGexc_i.close();
	} else
	{
		std::cout << "Error writing output file for _array_NGexc_i." << endl;
	}
	ofstream outfile__array_NGexc_I_max;
	outfile__array_NGexc_I_max.open(results_dir + "_array_NGexc_I_max_392997354", ios::binary | ios::out);
	if(outfile__array_NGexc_I_max.is_open())
	{
		outfile__array_NGexc_I_max.write(reinterpret_cast<char*>(_array_NGexc_I_max), 1*sizeof(_array_NGexc_I_max[0]));
		outfile__array_NGexc_I_max.close();
	} else
	{
		std::cout << "Error writing output file for _array_NGexc_I_max." << endl;
	}
	ofstream outfile__array_NGexc_lastspike;
	outfile__array_NGexc_lastspike.open(results_dir + "_array_NGexc_lastspike_3629583682", ios::binary | ios::out);
	if(outfile__array_NGexc_lastspike.is_open())
	{
		outfile__array_NGexc_lastspike.write(reinterpret_cast<char*>(_array_NGexc_lastspike), 1*sizeof(_array_NGexc_lastspike[0]));
		outfile__array_NGexc_lastspike.close();
	} else
	{
		std::cout << "Error writing output file for _array_NGexc_lastspike." << endl;
	}
	ofstream outfile__array_NGexc_m;
	outfile__array_NGexc_m.open(results_dir + "_array_NGexc_m_1323067197", ios::binary | ios::out);
	if(outfile__array_NGexc_m.is_open())
	{
		outfile__array_NGexc_m.write(reinterpret_cast<char*>(_array_NGexc_m), 1*sizeof(_array_NGexc_m[0]));
		outfile__array_NGexc_m.close();
	} else
	{
		std::cout << "Error writing output file for _array_NGexc_m." << endl;
	}
	ofstream outfile__array_NGexc_mask_noise;
	outfile__array_NGexc_mask_noise.open(results_dir + "_array_NGexc_mask_noise_2622314553", ios::binary | ios::out);
	if(outfile__array_NGexc_mask_noise.is_open())
	{
		outfile__array_NGexc_mask_noise.write(reinterpret_cast<char*>(_array_NGexc_mask_noise), 1*sizeof(_array_NGexc_mask_noise[0]));
		outfile__array_NGexc_mask_noise.close();
	} else
	{
		std::cout << "Error writing output file for _array_NGexc_mask_noise." << endl;
	}
	ofstream outfile__array_NGexc_n;
	outfile__array_NGexc_n.open(results_dir + "_array_NGexc_n_3621074567", ios::binary | ios::out);
	if(outfile__array_NGexc_n.is_open())
	{
		outfile__array_NGexc_n.write(reinterpret_cast<char*>(_array_NGexc_n), 1*sizeof(_array_NGexc_n[0]));
		outfile__array_NGexc_n.close();
	} else
	{
		std::cout << "Error writing output file for _array_NGexc_n." << endl;
	}
	ofstream outfile__array_NGexc_not_refractory;
	outfile__array_NGexc_not_refractory.open(results_dir + "_array_NGexc_not_refractory_137558986", ios::binary | ios::out);
	if(outfile__array_NGexc_not_refractory.is_open())
	{
		outfile__array_NGexc_not_refractory.write(reinterpret_cast<char*>(_array_NGexc_not_refractory), 1*sizeof(_array_NGexc_not_refractory[0]));
		outfile__array_NGexc_not_refractory.close();
	} else
	{
		std::cout << "Error writing output file for _array_NGexc_not_refractory." << endl;
	}
	ofstream outfile__array_NGexc_p;
	outfile__array_NGexc_p.open(results_dir + "_array_NGexc_p_769264612", ios::binary | ios::out);
	if(outfile__array_NGexc_p.is_open())
	{
		outfile__array_NGexc_p.write(reinterpret_cast<char*>(_array_NGexc_p), 1*sizeof(_array_NGexc_p[0]));
		outfile__array_NGexc_p.close();
	} else
	{
		std::cout << "Error writing output file for _array_NGexc_p." << endl;
	}
	ofstream outfile__array_NGexc_subgroup__sub_idx;
	outfile__array_NGexc_subgroup__sub_idx.open(results_dir + "_array_NGexc_subgroup__sub_idx_3268875373", ios::binary | ios::out);
	if(outfile__array_NGexc_subgroup__sub_idx.is_open())
	{
		outfile__array_NGexc_subgroup__sub_idx.write(reinterpret_cast<char*>(_array_NGexc_subgroup__sub_idx), 1*sizeof(_array_NGexc_subgroup__sub_idx[0]));
		outfile__array_NGexc_subgroup__sub_idx.close();
	} else
	{
		std::cout << "Error writing output file for _array_NGexc_subgroup__sub_idx." << endl;
	}
	ofstream outfile__array_NGexc_subgroup__sub_idx_1;
	outfile__array_NGexc_subgroup__sub_idx_1.open(results_dir + "_array_NGexc_subgroup__sub_idx_1_3886895134", ios::binary | ios::out);
	if(outfile__array_NGexc_subgroup__sub_idx_1.is_open())
	{
		outfile__array_NGexc_subgroup__sub_idx_1.write(reinterpret_cast<char*>(_array_NGexc_subgroup__sub_idx_1), 1*sizeof(_array_NGexc_subgroup__sub_idx_1[0]));
		outfile__array_NGexc_subgroup__sub_idx_1.close();
	} else
	{
		std::cout << "Error writing output file for _array_NGexc_subgroup__sub_idx_1." << endl;
	}
	ofstream outfile__array_NGexc_subgroup__sub_idx_2;
	outfile__array_NGexc_subgroup__sub_idx_2.open(results_dir + "_array_NGexc_subgroup__sub_idx_2_2124677540", ios::binary | ios::out);
	if(outfile__array_NGexc_subgroup__sub_idx_2.is_open())
	{
		outfile__array_NGexc_subgroup__sub_idx_2.write(reinterpret_cast<char*>(_array_NGexc_subgroup__sub_idx_2), 1*sizeof(_array_NGexc_subgroup__sub_idx_2[0]));
		outfile__array_NGexc_subgroup__sub_idx_2.close();
	} else
	{
		std::cout << "Error writing output file for _array_NGexc_subgroup__sub_idx_2." << endl;
	}
	ofstream outfile__array_NGexc_v;
	outfile__array_NGexc_v.open(results_dir + "_array_NGexc_v_3300503249", ios::binary | ios::out);
	if(outfile__array_NGexc_v.is_open())
	{
		outfile__array_NGexc_v.write(reinterpret_cast<char*>(_array_NGexc_v), 1*sizeof(_array_NGexc_v[0]));
		outfile__array_NGexc_v.close();
	} else
	{
		std::cout << "Error writing output file for _array_NGexc_v." << endl;
	}
	ofstream outfile__array_NGexc_x;
	outfile__array_NGexc_x.open(results_dir + "_array_NGexc_x_587301846", ios::binary | ios::out);
	if(outfile__array_NGexc_x.is_open())
	{
		outfile__array_NGexc_x.write(reinterpret_cast<char*>(_array_NGexc_x), 1*sizeof(_array_NGexc_x[0]));
		outfile__array_NGexc_x.close();
	} else
	{
		std::cout << "Error writing output file for _array_NGexc_x." << endl;
	}
	ofstream outfile__array_NGexc_y;
	outfile__array_NGexc_y.open(results_dir + "_array_NGexc_y_1409725248", ios::binary | ios::out);
	if(outfile__array_NGexc_y.is_open())
	{
		outfile__array_NGexc_y.write(reinterpret_cast<char*>(_array_NGexc_y), 1*sizeof(_array_NGexc_y[0]));
		outfile__array_NGexc_y.close();
	} else
	{
		std::cout << "Error writing output file for _array_NGexc_y." << endl;
	}
	ofstream outfile__array_NGexc_z;
	outfile__array_NGexc_z.open(results_dir + "_array_NGexc_z_3440370426", ios::binary | ios::out);
	if(outfile__array_NGexc_z.is_open())
	{
		outfile__array_NGexc_z.write(reinterpret_cast<char*>(_array_NGexc_z), 1*sizeof(_array_NGexc_z[0]));
		outfile__array_NGexc_z.close();
	} else
	{
		std::cout << "Error writing output file for _array_NGexc_z." << endl;
	}
	ofstream outfile__array_NGinh__spikespace;
	outfile__array_NGinh__spikespace.open(results_dir + "_array_NGinh__spikespace_1407747499", ios::binary | ios::out);
	if(outfile__array_NGinh__spikespace.is_open())
	{
		outfile__array_NGinh__spikespace.write(reinterpret_cast<char*>(_array_NGinh__spikespace), 2*sizeof(_array_NGinh__spikespace[0]));
		outfile__array_NGinh__spikespace.close();
	} else
	{
		std::cout << "Error writing output file for _array_NGinh__spikespace." << endl;
	}
	ofstream outfile__array_NGinh_clock_1_dt;
	outfile__array_NGinh_clock_1_dt.open(results_dir + "_array_NGinh_clock_1_dt_1125660832", ios::binary | ios::out);
	if(outfile__array_NGinh_clock_1_dt.is_open())
	{
		outfile__array_NGinh_clock_1_dt.write(reinterpret_cast<char*>(_array_NGinh_clock_1_dt), 1*sizeof(_array_NGinh_clock_1_dt[0]));
		outfile__array_NGinh_clock_1_dt.close();
	} else
	{
		std::cout << "Error writing output file for _array_NGinh_clock_1_dt." << endl;
	}
	ofstream outfile__array_NGinh_clock_1_t;
	outfile__array_NGinh_clock_1_t.open(results_dir + "_array_NGinh_clock_1_t_1661774430", ios::binary | ios::out);
	if(outfile__array_NGinh_clock_1_t.is_open())
	{
		outfile__array_NGinh_clock_1_t.write(reinterpret_cast<char*>(_array_NGinh_clock_1_t), 1*sizeof(_array_NGinh_clock_1_t[0]));
		outfile__array_NGinh_clock_1_t.close();
	} else
	{
		std::cout << "Error writing output file for _array_NGinh_clock_1_t." << endl;
	}
	ofstream outfile__array_NGinh_clock_1_timestep;
	outfile__array_NGinh_clock_1_timestep.open(results_dir + "_array_NGinh_clock_1_timestep_587794015", ios::binary | ios::out);
	if(outfile__array_NGinh_clock_1_timestep.is_open())
	{
		outfile__array_NGinh_clock_1_timestep.write(reinterpret_cast<char*>(_array_NGinh_clock_1_timestep), 1*sizeof(_array_NGinh_clock_1_timestep[0]));
		outfile__array_NGinh_clock_1_timestep.close();
	} else
	{
		std::cout << "Error writing output file for _array_NGinh_clock_1_timestep." << endl;
	}
	ofstream outfile__array_NGinh_h;
	outfile__array_NGinh_h.open(results_dir + "_array_NGinh_h_2187434257", ios::binary | ios::out);
	if(outfile__array_NGinh_h.is_open())
	{
		outfile__array_NGinh_h.write(reinterpret_cast<char*>(_array_NGinh_h), 1*sizeof(_array_NGinh_h[0]));
		outfile__array_NGinh_h.close();
	} else
	{
		std::cout << "Error writing output file for _array_NGinh_h." << endl;
	}
	ofstream outfile__array_NGinh_i;
	outfile__array_NGinh_i.open(results_dir + "_array_NGinh_i_4117145991", ios::binary | ios::out);
	if(outfile__array_NGinh_i.is_open())
	{
		outfile__array_NGinh_i.write(reinterpret_cast<char*>(_array_NGinh_i), 1*sizeof(_array_NGinh_i[0]));
		outfile__array_NGinh_i.close();
	} else
	{
		std::cout << "Error writing output file for _array_NGinh_i." << endl;
	}
	ofstream outfile__array_NGinh_I_max;
	outfile__array_NGinh_I_max.open(results_dir + "_array_NGinh_I_max_1233860264", ios::binary | ios::out);
	if(outfile__array_NGinh_I_max.is_open())
	{
		outfile__array_NGinh_I_max.write(reinterpret_cast<char*>(_array_NGinh_I_max), 1*sizeof(_array_NGinh_I_max[0]));
		outfile__array_NGinh_I_max.close();
	} else
	{
		std::cout << "Error writing output file for _array_NGinh_I_max." << endl;
	}
	ofstream outfile__array_NGinh_lastspike;
	outfile__array_NGinh_lastspike.open(results_dir + "_array_NGinh_lastspike_2196339681", ios::binary | ios::out);
	if(outfile__array_NGinh_lastspike.is_open())
	{
		outfile__array_NGinh_lastspike.write(reinterpret_cast<char*>(_array_NGinh_lastspike), 1*sizeof(_array_NGinh_lastspike[0]));
		outfile__array_NGinh_lastspike.close();
	} else
	{
		std::cout << "Error writing output file for _array_NGinh_lastspike." << endl;
	}
	ofstream outfile__array_NGinh_m;
	outfile__array_NGinh_m.open(results_dir + "_array_NGinh_m_4060835230", ios::binary | ios::out);
	if(outfile__array_NGinh_m.is_open())
	{
		outfile__array_NGinh_m.write(reinterpret_cast<char*>(_array_NGinh_m), 1*sizeof(_array_NGinh_m[0]));
		outfile__array_NGinh_m.close();
	} else
	{
		std::cout << "Error writing output file for _array_NGinh_m." << endl;
	}
	ofstream outfile__array_NGinh_mask_noise;
	outfile__array_NGinh_mask_noise.open(results_dir + "_array_NGinh_mask_noise_3553104925", ios::binary | ios::out);
	if(outfile__array_NGinh_mask_noise.is_open())
	{
		outfile__array_NGinh_mask_noise.write(reinterpret_cast<char*>(_array_NGinh_mask_noise), 1*sizeof(_array_NGinh_mask_noise[0]));
		outfile__array_NGinh_mask_noise.close();
	} else
	{
		std::cout << "Error writing output file for _array_NGinh_mask_noise." << endl;
	}
	ofstream outfile__array_NGinh_n;
	outfile__array_NGinh_n.open(results_dir + "_array_NGinh_n_1795308580", ios::binary | ios::out);
	if(outfile__array_NGinh_n.is_open())
	{
		outfile__array_NGinh_n.write(reinterpret_cast<char*>(_array_NGinh_n), 1*sizeof(_array_NGinh_n[0]));
		outfile__array_NGinh_n.close();
	} else
	{
		std::cout << "Error writing output file for _array_NGinh_n." << endl;
	}
	ofstream outfile__array_NGinh_not_refractory;
	outfile__array_NGinh_not_refractory.open(results_dir + "_array_NGinh_not_refractory_4269115144", ios::binary | ios::out);
	if(outfile__array_NGinh_not_refractory.is_open())
	{
		outfile__array_NGinh_not_refractory.write(reinterpret_cast<char*>(_array_NGinh_not_refractory), 1*sizeof(_array_NGinh_not_refractory[0]));
		outfile__array_NGinh_not_refractory.close();
	} else
	{
		std::cout << "Error writing output file for _array_NGinh_not_refractory." << endl;
	}
	ofstream outfile__array_NGinh_p;
	outfile__array_NGinh_p.open(results_dir + "_array_NGinh_p_2433548615", ios::binary | ios::out);
	if(outfile__array_NGinh_p.is_open())
	{
		outfile__array_NGinh_p.write(reinterpret_cast<char*>(_array_NGinh_p), 1*sizeof(_array_NGinh_p[0]));
		outfile__array_NGinh_p.close();
	} else
	{
		std::cout << "Error writing output file for _array_NGinh_p." << endl;
	}
	ofstream outfile__array_NGinh_subgroup__sub_idx;
	outfile__array_NGinh_subgroup__sub_idx.open(results_dir + "_array_NGinh_subgroup__sub_idx_1003402364", ios::binary | ios::out);
	if(outfile__array_NGinh_subgroup__sub_idx.is_open())
	{
		outfile__array_NGinh_subgroup__sub_idx.write(reinterpret_cast<char*>(_array_NGinh_subgroup__sub_idx), 1*sizeof(_array_NGinh_subgroup__sub_idx[0]));
		outfile__array_NGinh_subgroup__sub_idx.close();
	} else
	{
		std::cout << "Error writing output file for _array_NGinh_subgroup__sub_idx." << endl;
	}
	ofstream outfile__array_NGinh_subgroup__sub_idx_1;
	outfile__array_NGinh_subgroup__sub_idx_1.open(results_dir + "_array_NGinh_subgroup__sub_idx_1_1654794751", ios::binary | ios::out);
	if(outfile__array_NGinh_subgroup__sub_idx_1.is_open())
	{
		outfile__array_NGinh_subgroup__sub_idx_1.write(reinterpret_cast<char*>(_array_NGinh_subgroup__sub_idx_1), 1*sizeof(_array_NGinh_subgroup__sub_idx_1[0]));
		outfile__array_NGinh_subgroup__sub_idx_1.close();
	} else
	{
		std::cout << "Error writing output file for _array_NGinh_subgroup__sub_idx_1." << endl;
	}
	ofstream outfile__array_NGinh_subgroup__sub_idx_2;
	outfile__array_NGinh_subgroup__sub_idx_2.open(results_dir + "_array_NGinh_subgroup__sub_idx_2_4222318661", ios::binary | ios::out);
	if(outfile__array_NGinh_subgroup__sub_idx_2.is_open())
	{
		outfile__array_NGinh_subgroup__sub_idx_2.write(reinterpret_cast<char*>(_array_NGinh_subgroup__sub_idx_2), 1*sizeof(_array_NGinh_subgroup__sub_idx_2[0]));
		outfile__array_NGinh_subgroup__sub_idx_2.close();
	} else
	{
		std::cout << "Error writing output file for _array_NGinh_subgroup__sub_idx_2." << endl;
	}
	ofstream outfile__array_NGinh_v;
	outfile__array_NGinh_v.open(results_dir + "_array_NGinh_v_2020516978", ios::binary | ios::out);
	if(outfile__array_NGinh_v.is_open())
	{
		outfile__array_NGinh_v.write(reinterpret_cast<char*>(_array_NGinh_v), 1*sizeof(_array_NGinh_v[0]));
		outfile__array_NGinh_v.close();
	} else
	{
		std::cout << "Error writing output file for _array_NGinh_v." << endl;
	}
	ofstream outfile__array_NGinh_x;
	outfile__array_NGinh_x.open(results_dir + "_array_NGinh_x_2681637237", ios::binary | ios::out);
	if(outfile__array_NGinh_x.is_open())
	{
		outfile__array_NGinh_x.write(reinterpret_cast<char*>(_array_NGinh_x), 1*sizeof(_array_NGinh_x[0]));
		outfile__array_NGinh_x.close();
	} else
	{
		std::cout << "Error writing output file for _array_NGinh_x." << endl;
	}
	ofstream outfile__array_NGinh_y;
	outfile__array_NGinh_y.open(results_dir + "_array_NGinh_y_3906058723", ios::binary | ios::out);
	if(outfile__array_NGinh_y.is_open())
	{
		outfile__array_NGinh_y.write(reinterpret_cast<char*>(_array_NGinh_y), 1*sizeof(_array_NGinh_y[0]));
		outfile__array_NGinh_y.close();
	} else
	{
		std::cout << "Error writing output file for _array_NGinh_y." << endl;
	}
	ofstream outfile__array_NGinh_z;
	outfile__array_NGinh_z.open(results_dir + "_array_NGinh_z_1910040665", ios::binary | ios::out);
	if(outfile__array_NGinh_z.is_open())
	{
		outfile__array_NGinh_z.write(reinterpret_cast<char*>(_array_NGinh_z), 1*sizeof(_array_NGinh_z[0]));
		outfile__array_NGinh_z.close();
	} else
	{
		std::cout << "Error writing output file for _array_NGinh_z." << endl;
	}
	ofstream outfile__array_spikemonitor_2__source_idx;
	outfile__array_spikemonitor_2__source_idx.open(results_dir + "_array_spikemonitor_2__source_idx_1793786228", ios::binary | ios::out);
	if(outfile__array_spikemonitor_2__source_idx.is_open())
	{
		outfile__array_spikemonitor_2__source_idx.write(reinterpret_cast<char*>(_array_spikemonitor_2__source_idx), 1*sizeof(_array_spikemonitor_2__source_idx[0]));
		outfile__array_spikemonitor_2__source_idx.close();
	} else
	{
		std::cout << "Error writing output file for _array_spikemonitor_2__source_idx." << endl;
	}
	ofstream outfile__array_spikemonitor_2_count;
	outfile__array_spikemonitor_2_count.open(results_dir + "_array_spikemonitor_2_count_3621222387", ios::binary | ios::out);
	if(outfile__array_spikemonitor_2_count.is_open())
	{
		outfile__array_spikemonitor_2_count.write(reinterpret_cast<char*>(_array_spikemonitor_2_count), 1*sizeof(_array_spikemonitor_2_count[0]));
		outfile__array_spikemonitor_2_count.close();
	} else
	{
		std::cout << "Error writing output file for _array_spikemonitor_2_count." << endl;
	}
	ofstream outfile__array_spikemonitor_2_N;
	outfile__array_spikemonitor_2_N.open(results_dir + "_array_spikemonitor_2_N_2352936276", ios::binary | ios::out);
	if(outfile__array_spikemonitor_2_N.is_open())
	{
		outfile__array_spikemonitor_2_N.write(reinterpret_cast<char*>(_array_spikemonitor_2_N), 1*sizeof(_array_spikemonitor_2_N[0]));
		outfile__array_spikemonitor_2_N.close();
	} else
	{
		std::cout << "Error writing output file for _array_spikemonitor_2_N." << endl;
	}
	ofstream outfile__array_spikemonitor_3__source_idx;
	outfile__array_spikemonitor_3__source_idx.open(results_dir + "_array_spikemonitor_3__source_idx_3078478065", ios::binary | ios::out);
	if(outfile__array_spikemonitor_3__source_idx.is_open())
	{
		outfile__array_spikemonitor_3__source_idx.write(reinterpret_cast<char*>(_array_spikemonitor_3__source_idx), 1*sizeof(_array_spikemonitor_3__source_idx[0]));
		outfile__array_spikemonitor_3__source_idx.close();
	} else
	{
		std::cout << "Error writing output file for _array_spikemonitor_3__source_idx." << endl;
	}
	ofstream outfile__array_spikemonitor_3_count;
	outfile__array_spikemonitor_3_count.open(results_dir + "_array_spikemonitor_3_count_1906342983", ios::binary | ios::out);
	if(outfile__array_spikemonitor_3_count.is_open())
	{
		outfile__array_spikemonitor_3_count.write(reinterpret_cast<char*>(_array_spikemonitor_3_count), 1*sizeof(_array_spikemonitor_3_count[0]));
		outfile__array_spikemonitor_3_count.close();
	} else
	{
		std::cout << "Error writing output file for _array_spikemonitor_3_count." << endl;
	}
	ofstream outfile__array_spikemonitor_3_N;
	outfile__array_spikemonitor_3_N.open(results_dir + "_array_spikemonitor_3_N_2382143331", ios::binary | ios::out);
	if(outfile__array_spikemonitor_3_N.is_open())
	{
		outfile__array_spikemonitor_3_N.write(reinterpret_cast<char*>(_array_spikemonitor_3_N), 1*sizeof(_array_spikemonitor_3_N[0]));
		outfile__array_spikemonitor_3_N.close();
	} else
	{
		std::cout << "Error writing output file for _array_spikemonitor_3_N." << endl;
	}
	ofstream outfile__array_svmon_exc_add__indices;
	outfile__array_svmon_exc_add__indices.open(results_dir + "_array_svmon_exc_add__indices_284829085", ios::binary | ios::out);
	if(outfile__array_svmon_exc_add__indices.is_open())
	{
		outfile__array_svmon_exc_add__indices.write(reinterpret_cast<char*>(_array_svmon_exc_add__indices), 1*sizeof(_array_svmon_exc_add__indices[0]));
		outfile__array_svmon_exc_add__indices.close();
	} else
	{
		std::cout << "Error writing output file for _array_svmon_exc_add__indices." << endl;
	}
	ofstream outfile__array_svmon_exc_add_clock_1_dt;
	outfile__array_svmon_exc_add_clock_1_dt.open(results_dir + "_array_svmon_exc_add_clock_1_dt_4122997639", ios::binary | ios::out);
	if(outfile__array_svmon_exc_add_clock_1_dt.is_open())
	{
		outfile__array_svmon_exc_add_clock_1_dt.write(reinterpret_cast<char*>(_array_svmon_exc_add_clock_1_dt), 1*sizeof(_array_svmon_exc_add_clock_1_dt[0]));
		outfile__array_svmon_exc_add_clock_1_dt.close();
	} else
	{
		std::cout << "Error writing output file for _array_svmon_exc_add_clock_1_dt." << endl;
	}
	ofstream outfile__array_svmon_exc_add_clock_1_t;
	outfile__array_svmon_exc_add_clock_1_t.open(results_dir + "_array_svmon_exc_add_clock_1_t_2733289569", ios::binary | ios::out);
	if(outfile__array_svmon_exc_add_clock_1_t.is_open())
	{
		outfile__array_svmon_exc_add_clock_1_t.write(reinterpret_cast<char*>(_array_svmon_exc_add_clock_1_t), 1*sizeof(_array_svmon_exc_add_clock_1_t[0]));
		outfile__array_svmon_exc_add_clock_1_t.close();
	} else
	{
		std::cout << "Error writing output file for _array_svmon_exc_add_clock_1_t." << endl;
	}
	ofstream outfile__array_svmon_exc_add_clock_1_timestep;
	outfile__array_svmon_exc_add_clock_1_timestep.open(results_dir + "_array_svmon_exc_add_clock_1_timestep_3054127365", ios::binary | ios::out);
	if(outfile__array_svmon_exc_add_clock_1_timestep.is_open())
	{
		outfile__array_svmon_exc_add_clock_1_timestep.write(reinterpret_cast<char*>(_array_svmon_exc_add_clock_1_timestep), 1*sizeof(_array_svmon_exc_add_clock_1_timestep[0]));
		outfile__array_svmon_exc_add_clock_1_timestep.close();
	} else
	{
		std::cout << "Error writing output file for _array_svmon_exc_add_clock_1_timestep." << endl;
	}
	ofstream outfile__array_svmon_exc_add_N;
	outfile__array_svmon_exc_add_N.open(results_dir + "_array_svmon_exc_add_N_684165610", ios::binary | ios::out);
	if(outfile__array_svmon_exc_add_N.is_open())
	{
		outfile__array_svmon_exc_add_N.write(reinterpret_cast<char*>(_array_svmon_exc_add_N), 1*sizeof(_array_svmon_exc_add_N[0]));
		outfile__array_svmon_exc_add_N.close();
	} else
	{
		std::cout << "Error writing output file for _array_svmon_exc_add_N." << endl;
	}
	ofstream outfile__array_svmon_inh_add__indices;
	outfile__array_svmon_inh_add__indices.open(results_dir + "_array_svmon_inh_add__indices_2658698757", ios::binary | ios::out);
	if(outfile__array_svmon_inh_add__indices.is_open())
	{
		outfile__array_svmon_inh_add__indices.write(reinterpret_cast<char*>(_array_svmon_inh_add__indices), 1*sizeof(_array_svmon_inh_add__indices[0]));
		outfile__array_svmon_inh_add__indices.close();
	} else
	{
		std::cout << "Error writing output file for _array_svmon_inh_add__indices." << endl;
	}
	ofstream outfile__array_svmon_inh_add_clock_1_dt;
	outfile__array_svmon_inh_add_clock_1_dt.open(results_dir + "_array_svmon_inh_add_clock_1_dt_66621253", ios::binary | ios::out);
	if(outfile__array_svmon_inh_add_clock_1_dt.is_open())
	{
		outfile__array_svmon_inh_add_clock_1_dt.write(reinterpret_cast<char*>(_array_svmon_inh_add_clock_1_dt), 1*sizeof(_array_svmon_inh_add_clock_1_dt[0]));
		outfile__array_svmon_inh_add_clock_1_dt.close();
	} else
	{
		std::cout << "Error writing output file for _array_svmon_inh_add_clock_1_dt." << endl;
	}
	ofstream outfile__array_svmon_inh_add_clock_1_t;
	outfile__array_svmon_inh_add_clock_1_t.open(results_dir + "_array_svmon_inh_add_clock_1_t_1555048884", ios::binary | ios::out);
	if(outfile__array_svmon_inh_add_clock_1_t.is_open())
	{
		outfile__array_svmon_inh_add_clock_1_t.write(reinterpret_cast<char*>(_array_svmon_inh_add_clock_1_t), 1*sizeof(_array_svmon_inh_add_clock_1_t[0]));
		outfile__array_svmon_inh_add_clock_1_t.close();
	} else
	{
		std::cout << "Error writing output file for _array_svmon_inh_add_clock_1_t." << endl;
	}
	ofstream outfile__array_svmon_inh_add_clock_1_timestep;
	outfile__array_svmon_inh_add_clock_1_timestep.open(results_dir + "_array_svmon_inh_add_clock_1_timestep_1635947666", ios::binary | ios::out);
	if(outfile__array_svmon_inh_add_clock_1_timestep.is_open())
	{
		outfile__array_svmon_inh_add_clock_1_timestep.write(reinterpret_cast<char*>(_array_svmon_inh_add_clock_1_timestep), 1*sizeof(_array_svmon_inh_add_clock_1_timestep[0]));
		outfile__array_svmon_inh_add_clock_1_timestep.close();
	} else
	{
		std::cout << "Error writing output file for _array_svmon_inh_add_clock_1_timestep." << endl;
	}
	ofstream outfile__array_svmon_inh_add_N;
	outfile__array_svmon_inh_add_N.open(results_dir + "_array_svmon_inh_add_N_1981814440", ios::binary | ios::out);
	if(outfile__array_svmon_inh_add_N.is_open())
	{
		outfile__array_svmon_inh_add_N.write(reinterpret_cast<char*>(_array_svmon_inh_add_N), 1*sizeof(_array_svmon_inh_add_N[0]));
		outfile__array_svmon_inh_add_N.close();
	} else
	{
		std::cout << "Error writing output file for _array_svmon_inh_add_N." << endl;
	}

	ofstream outfile__dynamic_array_spikemonitor_2_i;
	outfile__dynamic_array_spikemonitor_2_i.open(results_dir + "_dynamic_array_spikemonitor_2_i_2642822512", ios::binary | ios::out);
	if(outfile__dynamic_array_spikemonitor_2_i.is_open())
	{
        if (! _dynamic_array_spikemonitor_2_i.empty() )
        {
			outfile__dynamic_array_spikemonitor_2_i.write(reinterpret_cast<char*>(&_dynamic_array_spikemonitor_2_i[0]), _dynamic_array_spikemonitor_2_i.size()*sizeof(_dynamic_array_spikemonitor_2_i[0]));
		    outfile__dynamic_array_spikemonitor_2_i.close();
		}
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_spikemonitor_2_i." << endl;
	}
	ofstream outfile__dynamic_array_spikemonitor_2_t;
	outfile__dynamic_array_spikemonitor_2_t.open(results_dir + "_dynamic_array_spikemonitor_2_t_4269812137", ios::binary | ios::out);
	if(outfile__dynamic_array_spikemonitor_2_t.is_open())
	{
        if (! _dynamic_array_spikemonitor_2_t.empty() )
        {
			outfile__dynamic_array_spikemonitor_2_t.write(reinterpret_cast<char*>(&_dynamic_array_spikemonitor_2_t[0]), _dynamic_array_spikemonitor_2_t.size()*sizeof(_dynamic_array_spikemonitor_2_t[0]));
		    outfile__dynamic_array_spikemonitor_2_t.close();
		}
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_spikemonitor_2_t." << endl;
	}
	ofstream outfile__dynamic_array_spikemonitor_3_i;
	outfile__dynamic_array_spikemonitor_3_i.open(results_dir + "_dynamic_array_spikemonitor_3_i_2621714247", ios::binary | ios::out);
	if(outfile__dynamic_array_spikemonitor_3_i.is_open())
	{
        if (! _dynamic_array_spikemonitor_3_i.empty() )
        {
			outfile__dynamic_array_spikemonitor_3_i.write(reinterpret_cast<char*>(&_dynamic_array_spikemonitor_3_i[0]), _dynamic_array_spikemonitor_3_i.size()*sizeof(_dynamic_array_spikemonitor_3_i[0]));
		    outfile__dynamic_array_spikemonitor_3_i.close();
		}
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_spikemonitor_3_i." << endl;
	}
	ofstream outfile__dynamic_array_spikemonitor_3_t;
	outfile__dynamic_array_spikemonitor_3_t.open(results_dir + "_dynamic_array_spikemonitor_3_t_4282532766", ios::binary | ios::out);
	if(outfile__dynamic_array_spikemonitor_3_t.is_open())
	{
        if (! _dynamic_array_spikemonitor_3_t.empty() )
        {
			outfile__dynamic_array_spikemonitor_3_t.write(reinterpret_cast<char*>(&_dynamic_array_spikemonitor_3_t[0]), _dynamic_array_spikemonitor_3_t.size()*sizeof(_dynamic_array_spikemonitor_3_t[0]));
		    outfile__dynamic_array_spikemonitor_3_t.close();
		}
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_spikemonitor_3_t." << endl;
	}
	ofstream outfile__dynamic_array_svmon_exc_add_t;
	outfile__dynamic_array_svmon_exc_add_t.open(results_dir + "_dynamic_array_svmon_exc_add_t_1466968228", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_exc_add_t.is_open())
	{
        if (! _dynamic_array_svmon_exc_add_t.empty() )
        {
			outfile__dynamic_array_svmon_exc_add_t.write(reinterpret_cast<char*>(&_dynamic_array_svmon_exc_add_t[0]), _dynamic_array_svmon_exc_add_t.size()*sizeof(_dynamic_array_svmon_exc_add_t[0]));
		    outfile__dynamic_array_svmon_exc_add_t.close();
		}
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_exc_add_t." << endl;
	}
	ofstream outfile__dynamic_array_svmon_inh_add_t;
	outfile__dynamic_array_svmon_inh_add_t.open(results_dir + "_dynamic_array_svmon_inh_add_t_160937958", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_inh_add_t.is_open())
	{
        if (! _dynamic_array_svmon_inh_add_t.empty() )
        {
			outfile__dynamic_array_svmon_inh_add_t.write(reinterpret_cast<char*>(&_dynamic_array_svmon_inh_add_t[0]), _dynamic_array_svmon_inh_add_t.size()*sizeof(_dynamic_array_svmon_inh_add_t[0]));
		    outfile__dynamic_array_svmon_inh_add_t.close();
		}
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_inh_add_t." << endl;
	}

	ofstream outfile__dynamic_array_svmon_exc_add_a_m;
	outfile__dynamic_array_svmon_exc_add_a_m.open(results_dir + "_dynamic_array_svmon_exc_add_a_m_2865408381", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_exc_add_a_m.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_exc_add_a_m.n; n++)
        {
            if (! _dynamic_array_svmon_exc_add_a_m(n).empty())
            {
                outfile__dynamic_array_svmon_exc_add_a_m.write(reinterpret_cast<char*>(&_dynamic_array_svmon_exc_add_a_m(n, 0)), _dynamic_array_svmon_exc_add_a_m.m*sizeof(_dynamic_array_svmon_exc_add_a_m(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_exc_add_a_m.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_exc_add_a_m." << endl;
	}
	ofstream outfile__dynamic_array_svmon_exc_add_b_m;
	outfile__dynamic_array_svmon_exc_add_b_m.open(results_dir + "_dynamic_array_svmon_exc_add_b_m_2827753252", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_exc_add_b_m.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_exc_add_b_m.n; n++)
        {
            if (! _dynamic_array_svmon_exc_add_b_m(n).empty())
            {
                outfile__dynamic_array_svmon_exc_add_b_m.write(reinterpret_cast<char*>(&_dynamic_array_svmon_exc_add_b_m(n, 0)), _dynamic_array_svmon_exc_add_b_m.m*sizeof(_dynamic_array_svmon_exc_add_b_m(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_exc_add_b_m.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_exc_add_b_m." << endl;
	}
	ofstream outfile__dynamic_array_svmon_exc_add_h;
	outfile__dynamic_array_svmon_exc_add_h.open(results_dir + "_dynamic_array_svmon_exc_add_h_1131508971", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_exc_add_h.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_exc_add_h.n; n++)
        {
            if (! _dynamic_array_svmon_exc_add_h(n).empty())
            {
                outfile__dynamic_array_svmon_exc_add_h.write(reinterpret_cast<char*>(&_dynamic_array_svmon_exc_add_h(n, 0)), _dynamic_array_svmon_exc_add_h.m*sizeof(_dynamic_array_svmon_exc_add_h(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_exc_add_h.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_exc_add_h." << endl;
	}
	ofstream outfile__dynamic_array_svmon_exc_add_I_inj;
	outfile__dynamic_array_svmon_exc_add_I_inj.open(results_dir + "_dynamic_array_svmon_exc_add_I_inj_3571689140", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_exc_add_I_inj.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_exc_add_I_inj.n; n++)
        {
            if (! _dynamic_array_svmon_exc_add_I_inj(n).empty())
            {
                outfile__dynamic_array_svmon_exc_add_I_inj.write(reinterpret_cast<char*>(&_dynamic_array_svmon_exc_add_I_inj(n, 0)), _dynamic_array_svmon_exc_add_I_inj.m*sizeof(_dynamic_array_svmon_exc_add_I_inj(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_exc_add_I_inj.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_exc_add_I_inj." << endl;
	}
	ofstream outfile__dynamic_array_svmon_exc_add_I_Kd;
	outfile__dynamic_array_svmon_exc_add_I_Kd.open(results_dir + "_dynamic_array_svmon_exc_add_I_Kd_278628625", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_exc_add_I_Kd.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_exc_add_I_Kd.n; n++)
        {
            if (! _dynamic_array_svmon_exc_add_I_Kd(n).empty())
            {
                outfile__dynamic_array_svmon_exc_add_I_Kd.write(reinterpret_cast<char*>(&_dynamic_array_svmon_exc_add_I_Kd(n, 0)), _dynamic_array_svmon_exc_add_I_Kd.m*sizeof(_dynamic_array_svmon_exc_add_I_Kd(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_exc_add_I_Kd.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_exc_add_I_Kd." << endl;
	}
	ofstream outfile__dynamic_array_svmon_exc_add_I_L;
	outfile__dynamic_array_svmon_exc_add_I_L.open(results_dir + "_dynamic_array_svmon_exc_add_I_L_3506286203", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_exc_add_I_L.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_exc_add_I_L.n; n++)
        {
            if (! _dynamic_array_svmon_exc_add_I_L(n).empty())
            {
                outfile__dynamic_array_svmon_exc_add_I_L.write(reinterpret_cast<char*>(&_dynamic_array_svmon_exc_add_I_L(n, 0)), _dynamic_array_svmon_exc_add_I_L.m*sizeof(_dynamic_array_svmon_exc_add_I_L(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_exc_add_I_L.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_exc_add_I_L." << endl;
	}
	ofstream outfile__dynamic_array_svmon_exc_add_I_M;
	outfile__dynamic_array_svmon_exc_add_I_M.open(results_dir + "_dynamic_array_svmon_exc_add_I_M_2818219757", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_exc_add_I_M.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_exc_add_I_M.n; n++)
        {
            if (! _dynamic_array_svmon_exc_add_I_M(n).empty())
            {
                outfile__dynamic_array_svmon_exc_add_I_M.write(reinterpret_cast<char*>(&_dynamic_array_svmon_exc_add_I_M(n, 0)), _dynamic_array_svmon_exc_add_I_M.m*sizeof(_dynamic_array_svmon_exc_add_I_M(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_exc_add_I_M.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_exc_add_I_M." << endl;
	}
	ofstream outfile__dynamic_array_svmon_exc_add_I_max;
	outfile__dynamic_array_svmon_exc_add_I_max.open(results_dir + "_dynamic_array_svmon_exc_add_I_max_2815124463", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_exc_add_I_max.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_exc_add_I_max.n; n++)
        {
            if (! _dynamic_array_svmon_exc_add_I_max(n).empty())
            {
                outfile__dynamic_array_svmon_exc_add_I_max.write(reinterpret_cast<char*>(&_dynamic_array_svmon_exc_add_I_max(n, 0)), _dynamic_array_svmon_exc_add_I_max.m*sizeof(_dynamic_array_svmon_exc_add_I_max(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_exc_add_I_max.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_exc_add_I_max." << endl;
	}
	ofstream outfile__dynamic_array_svmon_exc_add_I_Na;
	outfile__dynamic_array_svmon_exc_add_I_Na.open(results_dir + "_dynamic_array_svmon_exc_add_I_Na_495356379", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_exc_add_I_Na.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_exc_add_I_Na.n; n++)
        {
            if (! _dynamic_array_svmon_exc_add_I_Na(n).empty())
            {
                outfile__dynamic_array_svmon_exc_add_I_Na.write(reinterpret_cast<char*>(&_dynamic_array_svmon_exc_add_I_Na(n, 0)), _dynamic_array_svmon_exc_add_I_Na.m*sizeof(_dynamic_array_svmon_exc_add_I_Na(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_exc_add_I_Na.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_exc_add_I_Na." << endl;
	}
	ofstream outfile__dynamic_array_svmon_exc_add_m;
	outfile__dynamic_array_svmon_exc_add_m.open(results_dir + "_dynamic_array_svmon_exc_add_m_857440356", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_exc_add_m.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_exc_add_m.n; n++)
        {
            if (! _dynamic_array_svmon_exc_add_m(n).empty())
            {
                outfile__dynamic_array_svmon_exc_add_m.write(reinterpret_cast<char*>(&_dynamic_array_svmon_exc_add_m(n, 0)), _dynamic_array_svmon_exc_add_m.m*sizeof(_dynamic_array_svmon_exc_add_m(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_exc_add_m.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_exc_add_m." << endl;
	}
	ofstream outfile__dynamic_array_svmon_exc_add_n;
	outfile__dynamic_array_svmon_exc_add_n.open(results_dir + "_dynamic_array_svmon_exc_add_n_2853360094", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_exc_add_n.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_exc_add_n.n; n++)
        {
            if (! _dynamic_array_svmon_exc_add_n(n).empty())
            {
                outfile__dynamic_array_svmon_exc_add_n.write(reinterpret_cast<char*>(&_dynamic_array_svmon_exc_add_n(n, 0)), _dynamic_array_svmon_exc_add_n.m*sizeof(_dynamic_array_svmon_exc_add_n(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_exc_add_n.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_exc_add_n." << endl;
	}
	ofstream outfile__dynamic_array_svmon_exc_add_p;
	outfile__dynamic_array_svmon_exc_add_p.open(results_dir + "_dynamic_array_svmon_exc_add_p_1344138429", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_exc_add_p.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_exc_add_p.n; n++)
        {
            if (! _dynamic_array_svmon_exc_add_p(n).empty())
            {
                outfile__dynamic_array_svmon_exc_add_p.write(reinterpret_cast<char*>(&_dynamic_array_svmon_exc_add_p(n, 0)), _dynamic_array_svmon_exc_add_p.m*sizeof(_dynamic_array_svmon_exc_add_p(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_exc_add_p.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_exc_add_p." << endl;
	}
	ofstream outfile__dynamic_array_svmon_exc_add_p_inf;
	outfile__dynamic_array_svmon_exc_add_p_inf.open(results_dir + "_dynamic_array_svmon_exc_add_p_inf_1902391400", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_exc_add_p_inf.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_exc_add_p_inf.n; n++)
        {
            if (! _dynamic_array_svmon_exc_add_p_inf(n).empty())
            {
                outfile__dynamic_array_svmon_exc_add_p_inf.write(reinterpret_cast<char*>(&_dynamic_array_svmon_exc_add_p_inf(n, 0)), _dynamic_array_svmon_exc_add_p_inf.m*sizeof(_dynamic_array_svmon_exc_add_p_inf(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_exc_add_p_inf.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_exc_add_p_inf." << endl;
	}
	ofstream outfile__dynamic_array_svmon_exc_add_tau_h;
	outfile__dynamic_array_svmon_exc_add_tau_h.open(results_dir + "_dynamic_array_svmon_exc_add_tau_h_2697554331", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_exc_add_tau_h.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_exc_add_tau_h.n; n++)
        {
            if (! _dynamic_array_svmon_exc_add_tau_h(n).empty())
            {
                outfile__dynamic_array_svmon_exc_add_tau_h.write(reinterpret_cast<char*>(&_dynamic_array_svmon_exc_add_tau_h(n, 0)), _dynamic_array_svmon_exc_add_tau_h.m*sizeof(_dynamic_array_svmon_exc_add_tau_h(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_exc_add_tau_h.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_exc_add_tau_h." << endl;
	}
	ofstream outfile__dynamic_array_svmon_exc_add_tau_m;
	outfile__dynamic_array_svmon_exc_add_tau_m.open(results_dir + "_dynamic_array_svmon_exc_add_tau_m_3500383508", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_exc_add_tau_m.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_exc_add_tau_m.n; n++)
        {
            if (! _dynamic_array_svmon_exc_add_tau_m(n).empty())
            {
                outfile__dynamic_array_svmon_exc_add_tau_m.write(reinterpret_cast<char*>(&_dynamic_array_svmon_exc_add_tau_m(n, 0)), _dynamic_array_svmon_exc_add_tau_m.m*sizeof(_dynamic_array_svmon_exc_add_tau_m(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_exc_add_tau_m.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_exc_add_tau_m." << endl;
	}
	ofstream outfile__dynamic_array_svmon_exc_add_tau_n;
	outfile__dynamic_array_svmon_exc_add_tau_n.open(results_dir + "_dynamic_array_svmon_exc_add_tau_n_1235930286", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_exc_add_tau_n.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_exc_add_tau_n.n; n++)
        {
            if (! _dynamic_array_svmon_exc_add_tau_n(n).empty())
            {
                outfile__dynamic_array_svmon_exc_add_tau_n.write(reinterpret_cast<char*>(&_dynamic_array_svmon_exc_add_tau_n(n, 0)), _dynamic_array_svmon_exc_add_tau_n.m*sizeof(_dynamic_array_svmon_exc_add_tau_n(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_exc_add_tau_n.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_exc_add_tau_n." << endl;
	}
	ofstream outfile__dynamic_array_svmon_exc_add_tau_p;
	outfile__dynamic_array_svmon_exc_add_tau_p.open(results_dir + "_dynamic_array_svmon_exc_add_tau_p_3013997005", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_exc_add_tau_p.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_exc_add_tau_p.n; n++)
        {
            if (! _dynamic_array_svmon_exc_add_tau_p(n).empty())
            {
                outfile__dynamic_array_svmon_exc_add_tau_p.write(reinterpret_cast<char*>(&_dynamic_array_svmon_exc_add_tau_p(n, 0)), _dynamic_array_svmon_exc_add_tau_p.m*sizeof(_dynamic_array_svmon_exc_add_tau_p(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_exc_add_tau_p.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_exc_add_tau_p." << endl;
	}
	ofstream outfile__dynamic_array_svmon_exc_add_v;
	outfile__dynamic_array_svmon_exc_add_v.open(results_dir + "_dynamic_array_svmon_exc_add_v_3112061320", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_exc_add_v.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_exc_add_v.n; n++)
        {
            if (! _dynamic_array_svmon_exc_add_v(n).empty())
            {
                outfile__dynamic_array_svmon_exc_add_v.write(reinterpret_cast<char*>(&_dynamic_array_svmon_exc_add_v(n, 0)), _dynamic_array_svmon_exc_add_v.m*sizeof(_dynamic_array_svmon_exc_add_v(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_exc_add_v.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_exc_add_v." << endl;
	}
	ofstream outfile__dynamic_array_svmon_inh_add_a_m;
	outfile__dynamic_array_svmon_inh_add_a_m.open(results_dir + "_dynamic_array_svmon_inh_add_a_m_137712840", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_inh_add_a_m.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_inh_add_a_m.n; n++)
        {
            if (! _dynamic_array_svmon_inh_add_a_m(n).empty())
            {
                outfile__dynamic_array_svmon_inh_add_a_m.write(reinterpret_cast<char*>(&_dynamic_array_svmon_inh_add_a_m(n, 0)), _dynamic_array_svmon_inh_add_a_m.m*sizeof(_dynamic_array_svmon_inh_add_a_m(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_inh_add_a_m.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_inh_add_a_m." << endl;
	}
	ofstream outfile__dynamic_array_svmon_inh_add_b_m;
	outfile__dynamic_array_svmon_inh_add_b_m.open(results_dir + "_dynamic_array_svmon_inh_add_b_m_175368849", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_inh_add_b_m.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_inh_add_b_m.n; n++)
        {
            if (! _dynamic_array_svmon_inh_add_b_m(n).empty())
            {
                outfile__dynamic_array_svmon_inh_add_b_m.write(reinterpret_cast<char*>(&_dynamic_array_svmon_inh_add_b_m(n, 0)), _dynamic_array_svmon_inh_add_b_m.m*sizeof(_dynamic_array_svmon_inh_add_b_m(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_inh_add_b_m.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_inh_add_b_m." << endl;
	}
	ofstream outfile__dynamic_array_svmon_inh_add_h;
	outfile__dynamic_array_svmon_inh_add_h.open(results_dir + "_dynamic_array_svmon_inh_add_h_496429993", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_inh_add_h.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_inh_add_h.n; n++)
        {
            if (! _dynamic_array_svmon_inh_add_h(n).empty())
            {
                outfile__dynamic_array_svmon_inh_add_h.write(reinterpret_cast<char*>(&_dynamic_array_svmon_inh_add_h(n, 0)), _dynamic_array_svmon_inh_add_h.m*sizeof(_dynamic_array_svmon_inh_add_h(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_inh_add_h.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_inh_add_h." << endl;
	}
	ofstream outfile__dynamic_array_svmon_inh_add_I_inj;
	outfile__dynamic_array_svmon_inh_add_I_inj.open(results_dir + "_dynamic_array_svmon_inh_add_I_inj_2388512791", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_inh_add_I_inj.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_inh_add_I_inj.n; n++)
        {
            if (! _dynamic_array_svmon_inh_add_I_inj(n).empty())
            {
                outfile__dynamic_array_svmon_inh_add_I_inj.write(reinterpret_cast<char*>(&_dynamic_array_svmon_inh_add_I_inj(n, 0)), _dynamic_array_svmon_inh_add_I_inj.m*sizeof(_dynamic_array_svmon_inh_add_I_inj(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_inh_add_I_inj.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_inh_add_I_inj." << endl;
	}
	ofstream outfile__dynamic_array_svmon_inh_add_I_Kd;
	outfile__dynamic_array_svmon_inh_add_I_Kd.open(results_dir + "_dynamic_array_svmon_inh_add_I_Kd_2872193519", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_inh_add_I_Kd.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_inh_add_I_Kd.n; n++)
        {
            if (! _dynamic_array_svmon_inh_add_I_Kd(n).empty())
            {
                outfile__dynamic_array_svmon_inh_add_I_Kd.write(reinterpret_cast<char*>(&_dynamic_array_svmon_inh_add_I_Kd(n, 0)), _dynamic_array_svmon_inh_add_I_Kd.m*sizeof(_dynamic_array_svmon_inh_add_I_Kd(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_inh_add_I_Kd.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_inh_add_I_Kd." << endl;
	}
	ofstream outfile__dynamic_array_svmon_inh_add_I_L;
	outfile__dynamic_array_svmon_inh_add_I_L.open(results_dir + "_dynamic_array_svmon_inh_add_I_L_1912755150", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_inh_add_I_L.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_inh_add_I_L.n; n++)
        {
            if (! _dynamic_array_svmon_inh_add_I_L(n).empty())
            {
                outfile__dynamic_array_svmon_inh_add_I_L.write(reinterpret_cast<char*>(&_dynamic_array_svmon_inh_add_I_L(n, 0)), _dynamic_array_svmon_inh_add_I_L.m*sizeof(_dynamic_array_svmon_inh_add_I_L(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_inh_add_I_L.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_inh_add_I_L." << endl;
	}
	ofstream outfile__dynamic_array_svmon_inh_add_I_M;
	outfile__dynamic_array_svmon_inh_add_I_M.open(results_dir + "_dynamic_array_svmon_inh_add_I_M_84239192", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_inh_add_I_M.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_inh_add_I_M.n; n++)
        {
            if (! _dynamic_array_svmon_inh_add_I_M(n).empty())
            {
                outfile__dynamic_array_svmon_inh_add_I_M.write(reinterpret_cast<char*>(&_dynamic_array_svmon_inh_add_I_M(n, 0)), _dynamic_array_svmon_inh_add_I_M.m*sizeof(_dynamic_array_svmon_inh_add_I_M(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_inh_add_I_M.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_inh_add_I_M." << endl;
	}
	ofstream outfile__dynamic_array_svmon_inh_add_I_max;
	outfile__dynamic_array_svmon_inh_add_I_max.open(results_dir + "_dynamic_array_svmon_inh_add_I_max_4252308812", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_inh_add_I_max.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_inh_add_I_max.n; n++)
        {
            if (! _dynamic_array_svmon_inh_add_I_max(n).empty())
            {
                outfile__dynamic_array_svmon_inh_add_I_max.write(reinterpret_cast<char*>(&_dynamic_array_svmon_inh_add_I_max(n, 0)), _dynamic_array_svmon_inh_add_I_max.m*sizeof(_dynamic_array_svmon_inh_add_I_max(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_inh_add_I_max.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_inh_add_I_max." << endl;
	}
	ofstream outfile__dynamic_array_svmon_inh_add_I_Na;
	outfile__dynamic_array_svmon_inh_add_I_Na.open(results_dir + "_dynamic_array_svmon_inh_add_I_Na_2788110629", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_inh_add_I_Na.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_inh_add_I_Na.n; n++)
        {
            if (! _dynamic_array_svmon_inh_add_I_Na(n).empty())
            {
                outfile__dynamic_array_svmon_inh_add_I_Na.write(reinterpret_cast<char*>(&_dynamic_array_svmon_inh_add_I_Na(n, 0)), _dynamic_array_svmon_inh_add_I_Na.m*sizeof(_dynamic_array_svmon_inh_add_I_Na(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_inh_add_I_Na.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_inh_add_I_Na." << endl;
	}
	ofstream outfile__dynamic_array_svmon_inh_add_m;
	outfile__dynamic_array_svmon_inh_add_m.open(results_dir + "_dynamic_array_svmon_inh_add_m_1845239590", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_inh_add_m.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_inh_add_m.n; n++)
        {
            if (! _dynamic_array_svmon_inh_add_m(n).empty())
            {
                outfile__dynamic_array_svmon_inh_add_m.write(reinterpret_cast<char*>(&_dynamic_array_svmon_inh_add_m(n, 0)), _dynamic_array_svmon_inh_add_m.m*sizeof(_dynamic_array_svmon_inh_add_m(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_inh_add_m.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_inh_add_m." << endl;
	}
	ofstream outfile__dynamic_array_svmon_inh_add_n;
	outfile__dynamic_array_svmon_inh_add_n.open(results_dir + "_dynamic_array_svmon_inh_add_n_4109717148", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_inh_add_n.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_inh_add_n.n; n++)
        {
            if (! _dynamic_array_svmon_inh_add_n(n).empty())
            {
                outfile__dynamic_array_svmon_inh_add_n.write(reinterpret_cast<char*>(&_dynamic_array_svmon_inh_add_n(n, 0)), _dynamic_array_svmon_inh_add_n.m*sizeof(_dynamic_array_svmon_inh_add_n(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_inh_add_n.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_inh_add_n." << endl;
	}
	ofstream outfile__dynamic_array_svmon_inh_add_p;
	outfile__dynamic_array_svmon_inh_add_p.open(results_dir + "_dynamic_array_svmon_inh_add_p_251294719", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_inh_add_p.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_inh_add_p.n; n++)
        {
            if (! _dynamic_array_svmon_inh_add_p(n).empty())
            {
                outfile__dynamic_array_svmon_inh_add_p.write(reinterpret_cast<char*>(&_dynamic_array_svmon_inh_add_p(n, 0)), _dynamic_array_svmon_inh_add_p.m*sizeof(_dynamic_array_svmon_inh_add_p(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_inh_add_p.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_inh_add_p." << endl;
	}
	ofstream outfile__dynamic_array_svmon_inh_add_p_inf;
	outfile__dynamic_array_svmon_inh_add_p_inf.open(results_dir + "_dynamic_array_svmon_inh_add_p_inf_735725259", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_inh_add_p_inf.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_inh_add_p_inf.n; n++)
        {
            if (! _dynamic_array_svmon_inh_add_p_inf(n).empty())
            {
                outfile__dynamic_array_svmon_inh_add_p_inf.write(reinterpret_cast<char*>(&_dynamic_array_svmon_inh_add_p_inf(n, 0)), _dynamic_array_svmon_inh_add_p_inf.m*sizeof(_dynamic_array_svmon_inh_add_p_inf(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_inh_add_p_inf.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_inh_add_p_inf." << endl;
	}
	ofstream outfile__dynamic_array_svmon_inh_add_tau_h;
	outfile__dynamic_array_svmon_inh_add_tau_h.open(results_dir + "_dynamic_array_svmon_inh_add_tau_h_4202110776", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_inh_add_tau_h.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_inh_add_tau_h.n; n++)
        {
            if (! _dynamic_array_svmon_inh_add_tau_h(n).empty())
            {
                outfile__dynamic_array_svmon_inh_add_tau_h.write(reinterpret_cast<char*>(&_dynamic_array_svmon_inh_add_tau_h(n, 0)), _dynamic_array_svmon_inh_add_tau_h.m*sizeof(_dynamic_array_svmon_inh_add_tau_h(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_inh_add_tau_h.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_inh_add_tau_h." << endl;
	}
	ofstream outfile__dynamic_array_svmon_inh_add_tau_m;
	outfile__dynamic_array_svmon_inh_add_tau_m.open(results_dir + "_dynamic_array_svmon_inh_add_tau_m_2317216695", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_inh_add_tau_m.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_inh_add_tau_m.n; n++)
        {
            if (! _dynamic_array_svmon_inh_add_tau_m(n).empty())
            {
                outfile__dynamic_array_svmon_inh_add_tau_m.write(reinterpret_cast<char*>(&_dynamic_array_svmon_inh_add_tau_m(n, 0)), _dynamic_array_svmon_inh_add_tau_m.m*sizeof(_dynamic_array_svmon_inh_add_tau_m(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_inh_add_tau_m.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_inh_add_tau_m." << endl;
	}
	ofstream outfile__dynamic_array_svmon_inh_add_tau_n;
	outfile__dynamic_array_svmon_inh_add_tau_n.open(results_dir + "_dynamic_array_svmon_inh_add_tau_n_320125453", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_inh_add_tau_n.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_inh_add_tau_n.n; n++)
        {
            if (! _dynamic_array_svmon_inh_add_tau_n(n).empty())
            {
                outfile__dynamic_array_svmon_inh_add_tau_n.write(reinterpret_cast<char*>(&_dynamic_array_svmon_inh_add_tau_n(n, 0)), _dynamic_array_svmon_inh_add_tau_n.m*sizeof(_dynamic_array_svmon_inh_add_tau_n(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_inh_add_tau_n.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_inh_add_tau_n." << endl;
	}
	ofstream outfile__dynamic_array_svmon_inh_add_tau_p;
	outfile__dynamic_array_svmon_inh_add_tau_p.open(results_dir + "_dynamic_array_svmon_inh_add_tau_p_3910895470", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_inh_add_tau_p.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_inh_add_tau_p.n; n++)
        {
            if (! _dynamic_array_svmon_inh_add_tau_p(n).empty())
            {
                outfile__dynamic_array_svmon_inh_add_tau_p.write(reinterpret_cast<char*>(&_dynamic_array_svmon_inh_add_tau_p(n, 0)), _dynamic_array_svmon_inh_add_tau_p.m*sizeof(_dynamic_array_svmon_inh_add_tau_p(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_inh_add_tau_p.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_inh_add_tau_p." << endl;
	}
	ofstream outfile__dynamic_array_svmon_inh_add_v;
	outfile__dynamic_array_svmon_inh_add_v.open(results_dir + "_dynamic_array_svmon_inh_add_v_3885618890", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_inh_add_v.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_inh_add_v.n; n++)
        {
            if (! _dynamic_array_svmon_inh_add_v(n).empty())
            {
                outfile__dynamic_array_svmon_inh_add_v.write(reinterpret_cast<char*>(&_dynamic_array_svmon_inh_add_v(n, 0)), _dynamic_array_svmon_inh_add_v.m*sizeof(_dynamic_array_svmon_inh_add_v(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_inh_add_v.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_inh_add_v." << endl;
	}
	// Write last run info to disk
	ofstream outfile_last_run_info;
	outfile_last_run_info.open(results_dir + "last_run_info.txt", ios::out);
	if(outfile_last_run_info.is_open())
	{
		outfile_last_run_info << (Network::_last_run_time) << " " << (Network::_last_run_completed_fraction) << std::endl;
		outfile_last_run_info.close();
	} else
	{
	    std::cout << "Error writing last run info to file." << std::endl;
	}
}

void _dealloc_arrays()
{
	using namespace brian;


	// static arrays
	if(_static_array__index__array_NGexc_x!=0)
	{
		delete [] _static_array__index__array_NGexc_x;
		_static_array__index__array_NGexc_x = 0;
	}
	if(_static_array__index__array_NGexc_y!=0)
	{
		delete [] _static_array__index__array_NGexc_y;
		_static_array__index__array_NGexc_y = 0;
	}
	if(_static_array__index__array_NGexc_z!=0)
	{
		delete [] _static_array__index__array_NGexc_z;
		_static_array__index__array_NGexc_z = 0;
	}
	if(_static_array__index__array_NGinh_x!=0)
	{
		delete [] _static_array__index__array_NGinh_x;
		_static_array__index__array_NGinh_x = 0;
	}
	if(_static_array__index__array_NGinh_y!=0)
	{
		delete [] _static_array__index__array_NGinh_y;
		_static_array__index__array_NGinh_y = 0;
	}
	if(_static_array__index__array_NGinh_z!=0)
	{
		delete [] _static_array__index__array_NGinh_z;
		_static_array__index__array_NGinh_z = 0;
	}
	if(_static_array__value__array_NGexc_x!=0)
	{
		delete [] _static_array__value__array_NGexc_x;
		_static_array__value__array_NGexc_x = 0;
	}
	if(_static_array__value__array_NGexc_y!=0)
	{
		delete [] _static_array__value__array_NGexc_y;
		_static_array__value__array_NGexc_y = 0;
	}
	if(_static_array__value__array_NGexc_z!=0)
	{
		delete [] _static_array__value__array_NGexc_z;
		_static_array__value__array_NGexc_z = 0;
	}
	if(_static_array__value__array_NGinh_x!=0)
	{
		delete [] _static_array__value__array_NGinh_x;
		_static_array__value__array_NGinh_x = 0;
	}
	if(_static_array__value__array_NGinh_y!=0)
	{
		delete [] _static_array__value__array_NGinh_y;
		_static_array__value__array_NGinh_y = 0;
	}
	if(_static_array__value__array_NGinh_z!=0)
	{
		delete [] _static_array__value__array_NGinh_z;
		_static_array__value__array_NGinh_z = 0;
	}
	if(_timedarray_3_values!=0)
	{
		delete [] _timedarray_3_values;
		_timedarray_3_values = 0;
	}
	if(_timedarray_4_values!=0)
	{
		delete [] _timedarray_4_values;
		_timedarray_4_values = 0;
	}
}

