

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
    if (name == "hh_1._spikespace") {
        var_size = 2;
        data_size = 2*sizeof(int32_t);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<int32_t>(name, _array_hh_1__spikespace, var_size, (int32_t)atoi(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_hh_1__spikespace, data_size, s_value);
        }
        return;
    }
    if (name == "hh_1.h") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_hh_1_h, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_hh_1_h, data_size, s_value);
        }
        return;
    }
    if (name == "hh_1.m") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_hh_1_m, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_hh_1_m, data_size, s_value);
        }
        return;
    }
    if (name == "hh_1.n") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_hh_1_n, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_hh_1_n, data_size, s_value);
        }
        return;
    }
    if (name == "hh_1.n_Cl_E") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_hh_1_n_Cl_E, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_hh_1_n_Cl_E, data_size, s_value);
        }
        return;
    }
    if (name == "hh_1.n_Cl_N") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_hh_1_n_Cl_N, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_hh_1_n_Cl_N, data_size, s_value);
        }
        return;
    }
    if (name == "hh_1.n_K_E") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_hh_1_n_K_E, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_hh_1_n_K_E, data_size, s_value);
        }
        return;
    }
    if (name == "hh_1.n_K_N") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_hh_1_n_K_N, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_hh_1_n_K_N, data_size, s_value);
        }
        return;
    }
    if (name == "hh_1.n_Na_E") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_hh_1_n_Na_E, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_hh_1_n_Na_E, data_size, s_value);
        }
        return;
    }
    if (name == "hh_1.n_Na_N") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_hh_1_n_Na_N, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_hh_1_n_Na_N, data_size, s_value);
        }
        return;
    }
    if (name == "hh_1.v") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_hh_1_v, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_hh_1_v, data_size, s_value);
        }
        return;
    }
    // dynamic arrays (1d)
    if (name == "_timedarray_1.values") {
        var_size = 4;
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _timedarray_1_values, var_size, (double)atof(s_value.c_str()));


        } else {
            // set from file
            set_variable_from_file(name, _timedarray_1_values, data_size, s_value);
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
int32_t * _array_hh_1__spikespace;
const int _num__array_hh_1__spikespace = 2;
double * _array_hh_1_clock_dt;
const int _num__array_hh_1_clock_dt = 1;
double * _array_hh_1_clock_t;
const int _num__array_hh_1_clock_t = 1;
int64_t * _array_hh_1_clock_timestep;
const int _num__array_hh_1_clock_timestep = 1;
double * _array_hh_1_h;
const int _num__array_hh_1_h = 1;
int32_t * _array_hh_1_i;
const int _num__array_hh_1_i = 1;
double * _array_hh_1_m;
const int _num__array_hh_1_m = 1;
double * _array_hh_1_n;
const int _num__array_hh_1_n = 1;
double * _array_hh_1_n_Cl_E;
const int _num__array_hh_1_n_Cl_E = 1;
double * _array_hh_1_n_Cl_N;
const int _num__array_hh_1_n_Cl_N = 1;
double * _array_hh_1_n_K_E;
const int _num__array_hh_1_n_K_E = 1;
double * _array_hh_1_n_K_N;
const int _num__array_hh_1_n_K_N = 1;
double * _array_hh_1_n_Na_E;
const int _num__array_hh_1_n_Na_E = 1;
double * _array_hh_1_n_Na_N;
const int _num__array_hh_1_n_Na_N = 1;
double * _array_hh_1_v;
const int _num__array_hh_1_v = 1;
int32_t * _array_svmon__indices;
const int _num__array_svmon__indices = 1;
double * _array_svmon_C_Cl_N;
const int _num__array_svmon_C_Cl_N = (0, 1);
double * _array_svmon_C_K_N;
const int _num__array_svmon_C_K_N = (0, 1);
double * _array_svmon_C_Na_N;
const int _num__array_svmon_C_Na_N = (0, 1);
double * _array_svmon_clock_1_dt;
const int _num__array_svmon_clock_1_dt = 1;
double * _array_svmon_clock_1_t;
const int _num__array_svmon_clock_1_t = 1;
int64_t * _array_svmon_clock_1_timestep;
const int _num__array_svmon_clock_1_timestep = 1;
double * _array_svmon_E_Cl;
const int _num__array_svmon_E_Cl = (0, 1);
double * _array_svmon_E_K;
const int _num__array_svmon_E_K = (0, 1);
double * _array_svmon_E_Na;
const int _num__array_svmon_E_Na = (0, 1);
double * _array_svmon_I_Cl_L;
const int _num__array_svmon_I_Cl_L = (0, 1);
double * _array_svmon_I_K;
const int _num__array_svmon_I_K = (0, 1);
double * _array_svmon_I_KCC;
const int _num__array_svmon_I_KCC = (0, 1);
double * _array_svmon_I_Na;
const int _num__array_svmon_I_Na = (0, 1);
double * _array_svmon_I_Na_L;
const int _num__array_svmon_I_Na_L = (0, 1);
double * _array_svmon_I_NKP;
const int _num__array_svmon_I_NKP = (0, 1);
int32_t * _array_svmon_N;
const int _num__array_svmon_N = 1;
double * _array_svmon_v;
const int _num__array_svmon_v = (0, 1);

//////////////// dynamic arrays 1d /////////
std::vector<double> _dynamic_array_svmon_t;

//////////////// dynamic arrays 2d /////////
DynamicArray2D<double> _dynamic_array_svmon_C_Cl_N;
DynamicArray2D<double> _dynamic_array_svmon_C_K_N;
DynamicArray2D<double> _dynamic_array_svmon_C_Na_N;
DynamicArray2D<double> _dynamic_array_svmon_E_Cl;
DynamicArray2D<double> _dynamic_array_svmon_E_K;
DynamicArray2D<double> _dynamic_array_svmon_E_Na;
DynamicArray2D<double> _dynamic_array_svmon_I_Cl_L;
DynamicArray2D<double> _dynamic_array_svmon_I_K;
DynamicArray2D<double> _dynamic_array_svmon_I_KCC;
DynamicArray2D<double> _dynamic_array_svmon_I_Na;
DynamicArray2D<double> _dynamic_array_svmon_I_Na_L;
DynamicArray2D<double> _dynamic_array_svmon_I_NKP;
DynamicArray2D<double> _dynamic_array_svmon_v;

/////////////// static arrays /////////////
double * _timedarray_1_values;
const int _num__timedarray_1_values = 4;

//////////////// synapses /////////////////

//////////////// clocks ///////////////////
Clock hh_1_clock;  // attributes will be set in run.cpp
Clock svmon_clock_1;  // attributes will be set in run.cpp

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

	_array_hh_1__spikespace = new int32_t[2];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<2; i++) _array_hh_1__spikespace[i] = 0;

	_array_hh_1_clock_dt = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_hh_1_clock_dt[i] = 0;

	_array_hh_1_clock_t = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_hh_1_clock_t[i] = 0;

	_array_hh_1_clock_timestep = new int64_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_hh_1_clock_timestep[i] = 0;

	_array_hh_1_h = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_hh_1_h[i] = 0;

	_array_hh_1_i = new int32_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_hh_1_i[i] = 0;

	_array_hh_1_m = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_hh_1_m[i] = 0;

	_array_hh_1_n = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_hh_1_n[i] = 0;

	_array_hh_1_n_Cl_E = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_hh_1_n_Cl_E[i] = 0;

	_array_hh_1_n_Cl_N = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_hh_1_n_Cl_N[i] = 0;

	_array_hh_1_n_K_E = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_hh_1_n_K_E[i] = 0;

	_array_hh_1_n_K_N = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_hh_1_n_K_N[i] = 0;

	_array_hh_1_n_Na_E = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_hh_1_n_Na_E[i] = 0;

	_array_hh_1_n_Na_N = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_hh_1_n_Na_N[i] = 0;

	_array_hh_1_v = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_hh_1_v[i] = 0;

	_array_svmon__indices = new int32_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_svmon__indices[i] = 0;

	_array_svmon_clock_1_dt = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_svmon_clock_1_dt[i] = 0;

	_array_svmon_clock_1_t = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_svmon_clock_1_t[i] = 0;

	_array_svmon_clock_1_timestep = new int64_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_svmon_clock_1_timestep[i] = 0;

	_array_svmon_N = new int32_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_svmon_N[i] = 0;


	// Arrays initialized to an "arange"
	_array_hh_1_i = new int32_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_hh_1_i[i] = 0 + i;


	// static arrays
	_timedarray_1_values = new double[4];

	// Random number generator states
	for (int i=0; i<2; i++)
	    _mersenne_twister_states.push_back(new rk_state());
}

void _load_arrays()
{
	using namespace brian;

	ifstream f_timedarray_1_values;
	f_timedarray_1_values.open("static_arrays/_timedarray_1_values", ios::in | ios::binary);
	if(f_timedarray_1_values.is_open())
	{
		f_timedarray_1_values.read(reinterpret_cast<char*>(_timedarray_1_values), 4*sizeof(double));
	} else
	{
		std::cout << "Error opening static array _timedarray_1_values." << endl;
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
	ofstream outfile__array_hh_1__spikespace;
	outfile__array_hh_1__spikespace.open(results_dir + "_array_hh_1__spikespace_1223258136", ios::binary | ios::out);
	if(outfile__array_hh_1__spikespace.is_open())
	{
		outfile__array_hh_1__spikespace.write(reinterpret_cast<char*>(_array_hh_1__spikespace), 2*sizeof(_array_hh_1__spikespace[0]));
		outfile__array_hh_1__spikespace.close();
	} else
	{
		std::cout << "Error writing output file for _array_hh_1__spikespace." << endl;
	}
	ofstream outfile__array_hh_1_clock_dt;
	outfile__array_hh_1_clock_dt.open(results_dir + "_array_hh_1_clock_dt_2750978926", ios::binary | ios::out);
	if(outfile__array_hh_1_clock_dt.is_open())
	{
		outfile__array_hh_1_clock_dt.write(reinterpret_cast<char*>(_array_hh_1_clock_dt), 1*sizeof(_array_hh_1_clock_dt[0]));
		outfile__array_hh_1_clock_dt.close();
	} else
	{
		std::cout << "Error writing output file for _array_hh_1_clock_dt." << endl;
	}
	ofstream outfile__array_hh_1_clock_t;
	outfile__array_hh_1_clock_t.open(results_dir + "_array_hh_1_clock_t_1447459412", ios::binary | ios::out);
	if(outfile__array_hh_1_clock_t.is_open())
	{
		outfile__array_hh_1_clock_t.write(reinterpret_cast<char*>(_array_hh_1_clock_t), 1*sizeof(_array_hh_1_clock_t[0]));
		outfile__array_hh_1_clock_t.close();
	} else
	{
		std::cout << "Error writing output file for _array_hh_1_clock_t." << endl;
	}
	ofstream outfile__array_hh_1_clock_timestep;
	outfile__array_hh_1_clock_timestep.open(results_dir + "_array_hh_1_clock_timestep_1089539397", ios::binary | ios::out);
	if(outfile__array_hh_1_clock_timestep.is_open())
	{
		outfile__array_hh_1_clock_timestep.write(reinterpret_cast<char*>(_array_hh_1_clock_timestep), 1*sizeof(_array_hh_1_clock_timestep[0]));
		outfile__array_hh_1_clock_timestep.close();
	} else
	{
		std::cout << "Error writing output file for _array_hh_1_clock_timestep." << endl;
	}
	ofstream outfile__array_hh_1_h;
	outfile__array_hh_1_h.open(results_dir + "_array_hh_1_h_1853437229", ios::binary | ios::out);
	if(outfile__array_hh_1_h.is_open())
	{
		outfile__array_hh_1_h.write(reinterpret_cast<char*>(_array_hh_1_h), 1*sizeof(_array_hh_1_h[0]));
		outfile__array_hh_1_h.close();
	} else
	{
		std::cout << "Error writing output file for _array_hh_1_h." << endl;
	}
	ofstream outfile__array_hh_1_i;
	outfile__array_hh_1_i.open(results_dir + "_array_hh_1_i_427689403", ios::binary | ios::out);
	if(outfile__array_hh_1_i.is_open())
	{
		outfile__array_hh_1_i.write(reinterpret_cast<char*>(_array_hh_1_i), 1*sizeof(_array_hh_1_i[0]));
		outfile__array_hh_1_i.close();
	} else
	{
		std::cout << "Error writing output file for _array_hh_1_i." << endl;
	}
	ofstream outfile__array_hh_1_m;
	outfile__array_hh_1_m.open(results_dir + "_array_hh_1_m_504611234", ios::binary | ios::out);
	if(outfile__array_hh_1_m.is_open())
	{
		outfile__array_hh_1_m.write(reinterpret_cast<char*>(_array_hh_1_m), 1*sizeof(_array_hh_1_m[0]));
		outfile__array_hh_1_m.close();
	} else
	{
		std::cout << "Error writing output file for _array_hh_1_m." << endl;
	}
	ofstream outfile__array_hh_1_n;
	outfile__array_hh_1_n.open(results_dir + "_array_hh_1_n_2266664984", ios::binary | ios::out);
	if(outfile__array_hh_1_n.is_open())
	{
		outfile__array_hh_1_n.write(reinterpret_cast<char*>(_array_hh_1_n), 1*sizeof(_array_hh_1_n[0]));
		outfile__array_hh_1_n.close();
	} else
	{
		std::cout << "Error writing output file for _array_hh_1_n." << endl;
	}
	ofstream outfile__array_hh_1_n_Cl_E;
	outfile__array_hh_1_n_Cl_E.open(results_dir + "_array_hh_1_n_Cl_E_121342937", ios::binary | ios::out);
	if(outfile__array_hh_1_n_Cl_E.is_open())
	{
		outfile__array_hh_1_n_Cl_E.write(reinterpret_cast<char*>(_array_hh_1_n_Cl_E), 1*sizeof(_array_hh_1_n_Cl_E[0]));
		outfile__array_hh_1_n_Cl_E.close();
	} else
	{
		std::cout << "Error writing output file for _array_hh_1_n_Cl_E." << endl;
	}
	ofstream outfile__array_hh_1_n_Cl_N;
	outfile__array_hh_1_n_Cl_N.open(results_dir + "_array_hh_1_n_Cl_N_2431210065", ios::binary | ios::out);
	if(outfile__array_hh_1_n_Cl_N.is_open())
	{
		outfile__array_hh_1_n_Cl_N.write(reinterpret_cast<char*>(_array_hh_1_n_Cl_N), 1*sizeof(_array_hh_1_n_Cl_N[0]));
		outfile__array_hh_1_n_Cl_N.close();
	} else
	{
		std::cout << "Error writing output file for _array_hh_1_n_Cl_N." << endl;
	}
	ofstream outfile__array_hh_1_n_K_E;
	outfile__array_hh_1_n_K_E.open(results_dir + "_array_hh_1_n_K_E_1581219221", ios::binary | ios::out);
	if(outfile__array_hh_1_n_K_E.is_open())
	{
		outfile__array_hh_1_n_K_E.write(reinterpret_cast<char*>(_array_hh_1_n_K_E), 1*sizeof(_array_hh_1_n_K_E[0]));
		outfile__array_hh_1_n_K_E.close();
	} else
	{
		std::cout << "Error writing output file for _array_hh_1_n_K_E." << endl;
	}
	ofstream outfile__array_hh_1_n_K_N;
	outfile__array_hh_1_n_K_N.open(results_dir + "_array_hh_1_n_K_N_3387794461", ios::binary | ios::out);
	if(outfile__array_hh_1_n_K_N.is_open())
	{
		outfile__array_hh_1_n_K_N.write(reinterpret_cast<char*>(_array_hh_1_n_K_N), 1*sizeof(_array_hh_1_n_K_N[0]));
		outfile__array_hh_1_n_K_N.close();
	} else
	{
		std::cout << "Error writing output file for _array_hh_1_n_K_N." << endl;
	}
	ofstream outfile__array_hh_1_n_Na_E;
	outfile__array_hh_1_n_Na_E.open(results_dir + "_array_hh_1_n_Na_E_4253663319", ios::binary | ios::out);
	if(outfile__array_hh_1_n_Na_E.is_open())
	{
		outfile__array_hh_1_n_Na_E.write(reinterpret_cast<char*>(_array_hh_1_n_Na_E), 1*sizeof(_array_hh_1_n_Na_E[0]));
		outfile__array_hh_1_n_Na_E.close();
	} else
	{
		std::cout << "Error writing output file for _array_hh_1_n_Na_E." << endl;
	}
	ofstream outfile__array_hh_1_n_Na_N;
	outfile__array_hh_1_n_Na_N.open(results_dir + "_array_hh_1_n_Na_N_1784355295", ios::binary | ios::out);
	if(outfile__array_hh_1_n_Na_N.is_open())
	{
		outfile__array_hh_1_n_Na_N.write(reinterpret_cast<char*>(_array_hh_1_n_Na_N), 1*sizeof(_array_hh_1_n_Na_N[0]));
		outfile__array_hh_1_n_Na_N.close();
	} else
	{
		std::cout << "Error writing output file for _array_hh_1_n_Na_N." << endl;
	}
	ofstream outfile__array_hh_1_v;
	outfile__array_hh_1_v.open(results_dir + "_array_hh_1_v_2490763342", ios::binary | ios::out);
	if(outfile__array_hh_1_v.is_open())
	{
		outfile__array_hh_1_v.write(reinterpret_cast<char*>(_array_hh_1_v), 1*sizeof(_array_hh_1_v[0]));
		outfile__array_hh_1_v.close();
	} else
	{
		std::cout << "Error writing output file for _array_hh_1_v." << endl;
	}
	ofstream outfile__array_svmon__indices;
	outfile__array_svmon__indices.open(results_dir + "_array_svmon__indices_2661819678", ios::binary | ios::out);
	if(outfile__array_svmon__indices.is_open())
	{
		outfile__array_svmon__indices.write(reinterpret_cast<char*>(_array_svmon__indices), 1*sizeof(_array_svmon__indices[0]));
		outfile__array_svmon__indices.close();
	} else
	{
		std::cout << "Error writing output file for _array_svmon__indices." << endl;
	}
	ofstream outfile__array_svmon_clock_1_dt;
	outfile__array_svmon_clock_1_dt.open(results_dir + "_array_svmon_clock_1_dt_3843143261", ios::binary | ios::out);
	if(outfile__array_svmon_clock_1_dt.is_open())
	{
		outfile__array_svmon_clock_1_dt.write(reinterpret_cast<char*>(_array_svmon_clock_1_dt), 1*sizeof(_array_svmon_clock_1_dt[0]));
		outfile__array_svmon_clock_1_dt.close();
	} else
	{
		std::cout << "Error writing output file for _array_svmon_clock_1_dt." << endl;
	}
	ofstream outfile__array_svmon_clock_1_t;
	outfile__array_svmon_clock_1_t.open(results_dir + "_array_svmon_clock_1_t_3604295931", ios::binary | ios::out);
	if(outfile__array_svmon_clock_1_t.is_open())
	{
		outfile__array_svmon_clock_1_t.write(reinterpret_cast<char*>(_array_svmon_clock_1_t), 1*sizeof(_array_svmon_clock_1_t[0]));
		outfile__array_svmon_clock_1_t.close();
	} else
	{
		std::cout << "Error writing output file for _array_svmon_clock_1_t." << endl;
	}
	ofstream outfile__array_svmon_clock_1_timestep;
	outfile__array_svmon_clock_1_timestep.open(results_dir + "_array_svmon_clock_1_timestep_607022131", ios::binary | ios::out);
	if(outfile__array_svmon_clock_1_timestep.is_open())
	{
		outfile__array_svmon_clock_1_timestep.write(reinterpret_cast<char*>(_array_svmon_clock_1_timestep), 1*sizeof(_array_svmon_clock_1_timestep[0]));
		outfile__array_svmon_clock_1_timestep.close();
	} else
	{
		std::cout << "Error writing output file for _array_svmon_clock_1_timestep." << endl;
	}
	ofstream outfile__array_svmon_N;
	outfile__array_svmon_N.open(results_dir + "_array_svmon_N_525316449", ios::binary | ios::out);
	if(outfile__array_svmon_N.is_open())
	{
		outfile__array_svmon_N.write(reinterpret_cast<char*>(_array_svmon_N), 1*sizeof(_array_svmon_N[0]));
		outfile__array_svmon_N.close();
	} else
	{
		std::cout << "Error writing output file for _array_svmon_N." << endl;
	}

	ofstream outfile__dynamic_array_svmon_t;
	outfile__dynamic_array_svmon_t.open(results_dir + "_dynamic_array_svmon_t_2014749096", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_t.is_open())
	{
        if (! _dynamic_array_svmon_t.empty() )
        {
			outfile__dynamic_array_svmon_t.write(reinterpret_cast<char*>(&_dynamic_array_svmon_t[0]), _dynamic_array_svmon_t.size()*sizeof(_dynamic_array_svmon_t[0]));
		    outfile__dynamic_array_svmon_t.close();
		}
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_t." << endl;
	}

	ofstream outfile__dynamic_array_svmon_C_Cl_N;
	outfile__dynamic_array_svmon_C_Cl_N.open(results_dir + "_dynamic_array_svmon_C_Cl_N_1274607234", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_C_Cl_N.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_C_Cl_N.n; n++)
        {
            if (! _dynamic_array_svmon_C_Cl_N(n).empty())
            {
                outfile__dynamic_array_svmon_C_Cl_N.write(reinterpret_cast<char*>(&_dynamic_array_svmon_C_Cl_N(n, 0)), _dynamic_array_svmon_C_Cl_N.m*sizeof(_dynamic_array_svmon_C_Cl_N(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_C_Cl_N.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_C_Cl_N." << endl;
	}
	ofstream outfile__dynamic_array_svmon_C_K_N;
	outfile__dynamic_array_svmon_C_K_N.open(results_dir + "_dynamic_array_svmon_C_K_N_1669374263", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_C_K_N.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_C_K_N.n; n++)
        {
            if (! _dynamic_array_svmon_C_K_N(n).empty())
            {
                outfile__dynamic_array_svmon_C_K_N.write(reinterpret_cast<char*>(&_dynamic_array_svmon_C_K_N(n, 0)), _dynamic_array_svmon_C_K_N.m*sizeof(_dynamic_array_svmon_C_K_N(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_C_K_N.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_C_K_N." << endl;
	}
	ofstream outfile__dynamic_array_svmon_C_Na_N;
	outfile__dynamic_array_svmon_C_Na_N.open(results_dir + "_dynamic_array_svmon_C_Na_N_2974465292", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_C_Na_N.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_C_Na_N.n; n++)
        {
            if (! _dynamic_array_svmon_C_Na_N(n).empty())
            {
                outfile__dynamic_array_svmon_C_Na_N.write(reinterpret_cast<char*>(&_dynamic_array_svmon_C_Na_N(n, 0)), _dynamic_array_svmon_C_Na_N.m*sizeof(_dynamic_array_svmon_C_Na_N(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_C_Na_N.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_C_Na_N." << endl;
	}
	ofstream outfile__dynamic_array_svmon_E_Cl;
	outfile__dynamic_array_svmon_E_Cl.open(results_dir + "_dynamic_array_svmon_E_Cl_3687131831", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_E_Cl.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_E_Cl.n; n++)
        {
            if (! _dynamic_array_svmon_E_Cl(n).empty())
            {
                outfile__dynamic_array_svmon_E_Cl.write(reinterpret_cast<char*>(&_dynamic_array_svmon_E_Cl(n, 0)), _dynamic_array_svmon_E_Cl.m*sizeof(_dynamic_array_svmon_E_Cl(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_E_Cl.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_E_Cl." << endl;
	}
	ofstream outfile__dynamic_array_svmon_E_K;
	outfile__dynamic_array_svmon_E_K.open(results_dir + "_dynamic_array_svmon_E_K_2136119634", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_E_K.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_E_K.n; n++)
        {
            if (! _dynamic_array_svmon_E_K(n).empty())
            {
                outfile__dynamic_array_svmon_E_K.write(reinterpret_cast<char*>(&_dynamic_array_svmon_E_K(n, 0)), _dynamic_array_svmon_E_K.m*sizeof(_dynamic_array_svmon_E_K(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_E_K.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_E_K." << endl;
	}
	ofstream outfile__dynamic_array_svmon_E_Na;
	outfile__dynamic_array_svmon_E_Na.open(results_dir + "_dynamic_array_svmon_E_Na_282732615", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_E_Na.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_E_Na.n; n++)
        {
            if (! _dynamic_array_svmon_E_Na(n).empty())
            {
                outfile__dynamic_array_svmon_E_Na.write(reinterpret_cast<char*>(&_dynamic_array_svmon_E_Na(n, 0)), _dynamic_array_svmon_E_Na.m*sizeof(_dynamic_array_svmon_E_Na(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_E_Na.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_E_Na." << endl;
	}
	ofstream outfile__dynamic_array_svmon_I_Cl_L;
	outfile__dynamic_array_svmon_I_Cl_L.open(results_dir + "_dynamic_array_svmon_I_Cl_L_74298568", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_I_Cl_L.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_I_Cl_L.n; n++)
        {
            if (! _dynamic_array_svmon_I_Cl_L(n).empty())
            {
                outfile__dynamic_array_svmon_I_Cl_L.write(reinterpret_cast<char*>(&_dynamic_array_svmon_I_Cl_L(n, 0)), _dynamic_array_svmon_I_Cl_L.m*sizeof(_dynamic_array_svmon_I_Cl_L(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_I_Cl_L.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_I_Cl_L." << endl;
	}
	ofstream outfile__dynamic_array_svmon_I_K;
	outfile__dynamic_array_svmon_I_K.open(results_dir + "_dynamic_array_svmon_I_K_1984454710", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_I_K.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_I_K.n; n++)
        {
            if (! _dynamic_array_svmon_I_K(n).empty())
            {
                outfile__dynamic_array_svmon_I_K.write(reinterpret_cast<char*>(&_dynamic_array_svmon_I_K(n, 0)), _dynamic_array_svmon_I_K.m*sizeof(_dynamic_array_svmon_I_K(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_I_K.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_I_K." << endl;
	}
	ofstream outfile__dynamic_array_svmon_I_KCC;
	outfile__dynamic_array_svmon_I_KCC.open(results_dir + "_dynamic_array_svmon_I_KCC_2985728118", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_I_KCC.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_I_KCC.n; n++)
        {
            if (! _dynamic_array_svmon_I_KCC(n).empty())
            {
                outfile__dynamic_array_svmon_I_KCC.write(reinterpret_cast<char*>(&_dynamic_array_svmon_I_KCC(n, 0)), _dynamic_array_svmon_I_KCC.m*sizeof(_dynamic_array_svmon_I_KCC(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_I_KCC.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_I_KCC." << endl;
	}
	ofstream outfile__dynamic_array_svmon_I_Na;
	outfile__dynamic_array_svmon_I_Na.open(results_dir + "_dynamic_array_svmon_I_Na_1510774783", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_I_Na.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_I_Na.n; n++)
        {
            if (! _dynamic_array_svmon_I_Na(n).empty())
            {
                outfile__dynamic_array_svmon_I_Na.write(reinterpret_cast<char*>(&_dynamic_array_svmon_I_Na(n, 0)), _dynamic_array_svmon_I_Na.m*sizeof(_dynamic_array_svmon_I_Na(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_I_Na.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_I_Na." << endl;
	}
	ofstream outfile__dynamic_array_svmon_I_Na_L;
	outfile__dynamic_array_svmon_I_Na_L.open(results_dir + "_dynamic_array_svmon_I_Na_L_4276092742", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_I_Na_L.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_I_Na_L.n; n++)
        {
            if (! _dynamic_array_svmon_I_Na_L(n).empty())
            {
                outfile__dynamic_array_svmon_I_Na_L.write(reinterpret_cast<char*>(&_dynamic_array_svmon_I_Na_L(n, 0)), _dynamic_array_svmon_I_Na_L.m*sizeof(_dynamic_array_svmon_I_Na_L(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_I_Na_L.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_I_Na_L." << endl;
	}
	ofstream outfile__dynamic_array_svmon_I_NKP;
	outfile__dynamic_array_svmon_I_NKP.open(results_dir + "_dynamic_array_svmon_I_NKP_4217016651", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_I_NKP.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_I_NKP.n; n++)
        {
            if (! _dynamic_array_svmon_I_NKP(n).empty())
            {
                outfile__dynamic_array_svmon_I_NKP.write(reinterpret_cast<char*>(&_dynamic_array_svmon_I_NKP(n, 0)), _dynamic_array_svmon_I_NKP.m*sizeof(_dynamic_array_svmon_I_NKP(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_I_NKP.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_I_NKP." << endl;
	}
	ofstream outfile__dynamic_array_svmon_v;
	outfile__dynamic_array_svmon_v.open(results_dir + "_dynamic_array_svmon_v_2518204548", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_v.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_v.n; n++)
        {
            if (! _dynamic_array_svmon_v(n).empty())
            {
                outfile__dynamic_array_svmon_v.write(reinterpret_cast<char*>(&_dynamic_array_svmon_v(n, 0)), _dynamic_array_svmon_v.m*sizeof(_dynamic_array_svmon_v(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_v.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_v." << endl;
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
	if(_timedarray_1_values!=0)
	{
		delete [] _timedarray_1_values;
		_timedarray_1_values = 0;
	}
}

