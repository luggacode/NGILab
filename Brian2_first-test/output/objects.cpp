

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
Network network;

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
    if (name == "hh._spikespace") {
        var_size = 2;
        data_size = 2*sizeof(int32_t);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<int32_t>(name, _array_hh__spikespace, var_size, (int32_t)atoi(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_hh__spikespace, data_size, s_value);
        }
        return;
    }
    if (name == "hh.h") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_hh_h, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_hh_h, data_size, s_value);
        }
        return;
    }
    if (name == "hh.m") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_hh_m, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_hh_m, data_size, s_value);
        }
        return;
    }
    if (name == "hh.n") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_hh_n, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_hh_n, data_size, s_value);
        }
        return;
    }
    if (name == "hh.n_Cl_E") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_hh_n_Cl_E, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_hh_n_Cl_E, data_size, s_value);
        }
        return;
    }
    if (name == "hh.n_Cl_N") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_hh_n_Cl_N, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_hh_n_Cl_N, data_size, s_value);
        }
        return;
    }
    if (name == "hh.n_K_E") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_hh_n_K_E, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_hh_n_K_E, data_size, s_value);
        }
        return;
    }
    if (name == "hh.n_K_N") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_hh_n_K_N, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_hh_n_K_N, data_size, s_value);
        }
        return;
    }
    if (name == "hh.n_Na_E") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_hh_n_Na_E, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_hh_n_Na_E, data_size, s_value);
        }
        return;
    }
    if (name == "hh.n_Na_N") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_hh_n_Na_N, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_hh_n_Na_N, data_size, s_value);
        }
        return;
    }
    if (name == "hh.v") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_hh_v, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_hh_v, data_size, s_value);
        }
        return;
    }
    // dynamic arrays (1d)
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
int32_t * _array_hh__spikespace;
const int _num__array_hh__spikespace = 2;
double * _array_hh_clock_dt;
const int _num__array_hh_clock_dt = 1;
double * _array_hh_clock_t;
const int _num__array_hh_clock_t = 1;
int64_t * _array_hh_clock_timestep;
const int _num__array_hh_clock_timestep = 1;
double * _array_hh_h;
const int _num__array_hh_h = 1;
int32_t * _array_hh_i;
const int _num__array_hh_i = 1;
double * _array_hh_m;
const int _num__array_hh_m = 1;
double * _array_hh_n;
const int _num__array_hh_n = 1;
double * _array_hh_n_Cl_E;
const int _num__array_hh_n_Cl_E = 1;
double * _array_hh_n_Cl_N;
const int _num__array_hh_n_Cl_N = 1;
double * _array_hh_n_K_E;
const int _num__array_hh_n_K_E = 1;
double * _array_hh_n_K_N;
const int _num__array_hh_n_K_N = 1;
double * _array_hh_n_Na_E;
const int _num__array_hh_n_Na_E = 1;
double * _array_hh_n_Na_N;
const int _num__array_hh_n_Na_N = 1;
double * _array_hh_v;
const int _num__array_hh_v = 1;
int32_t * _array_svmon__indices;
const int _num__array_svmon__indices = 1;
double * _array_svmon_C_Cl_N;
const int _num__array_svmon_C_Cl_N = (0, 1);
double * _array_svmon_C_K_N;
const int _num__array_svmon_C_K_N = (0, 1);
double * _array_svmon_C_Na_N;
const int _num__array_svmon_C_Na_N = (0, 1);
double * _array_svmon_clock_dt;
const int _num__array_svmon_clock_dt = 1;
double * _array_svmon_clock_t;
const int _num__array_svmon_clock_t = 1;
int64_t * _array_svmon_clock_timestep;
const int _num__array_svmon_clock_timestep = 1;
double * _array_svmon_E_Cl;
const int _num__array_svmon_E_Cl = (0, 1);
double * _array_svmon_I_Cl;
const int _num__array_svmon_I_Cl = (0, 1);
double * _array_svmon_I_K;
const int _num__array_svmon_I_K = (0, 1);
double * _array_svmon_I_Na;
const int _num__array_svmon_I_Na = (0, 1);
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
DynamicArray2D<double> _dynamic_array_svmon_I_Cl;
DynamicArray2D<double> _dynamic_array_svmon_I_K;
DynamicArray2D<double> _dynamic_array_svmon_I_Na;
DynamicArray2D<double> _dynamic_array_svmon_v;

/////////////// static arrays /////////////

//////////////// synapses /////////////////

//////////////// clocks ///////////////////
Clock hh_clock;  // attributes will be set in run.cpp
Clock svmon_clock;  // attributes will be set in run.cpp

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

	_array_hh__spikespace = new int32_t[2];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<2; i++) _array_hh__spikespace[i] = 0;

	_array_hh_clock_dt = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_hh_clock_dt[i] = 0;

	_array_hh_clock_t = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_hh_clock_t[i] = 0;

	_array_hh_clock_timestep = new int64_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_hh_clock_timestep[i] = 0;

	_array_hh_h = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_hh_h[i] = 0;

	_array_hh_i = new int32_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_hh_i[i] = 0;

	_array_hh_m = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_hh_m[i] = 0;

	_array_hh_n = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_hh_n[i] = 0;

	_array_hh_n_Cl_E = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_hh_n_Cl_E[i] = 0;

	_array_hh_n_Cl_N = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_hh_n_Cl_N[i] = 0;

	_array_hh_n_K_E = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_hh_n_K_E[i] = 0;

	_array_hh_n_K_N = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_hh_n_K_N[i] = 0;

	_array_hh_n_Na_E = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_hh_n_Na_E[i] = 0;

	_array_hh_n_Na_N = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_hh_n_Na_N[i] = 0;

	_array_hh_v = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_hh_v[i] = 0;

	_array_svmon__indices = new int32_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_svmon__indices[i] = 0;

	_array_svmon_clock_dt = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_svmon_clock_dt[i] = 0;

	_array_svmon_clock_t = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_svmon_clock_t[i] = 0;

	_array_svmon_clock_timestep = new int64_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_svmon_clock_timestep[i] = 0;

	_array_svmon_N = new int32_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_svmon_N[i] = 0;


	// Arrays initialized to an "arange"
	_array_hh_i = new int32_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_hh_i[i] = 0 + i;


	// static arrays

	// Random number generator states
	for (int i=0; i<2; i++)
	    _mersenne_twister_states.push_back(new rk_state());
}

void _load_arrays()
{
	using namespace brian;

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
	ofstream outfile__array_hh__spikespace;
	outfile__array_hh__spikespace.open(results_dir + "_array_hh__spikespace_1182925430", ios::binary | ios::out);
	if(outfile__array_hh__spikespace.is_open())
	{
		outfile__array_hh__spikespace.write(reinterpret_cast<char*>(_array_hh__spikespace), 2*sizeof(_array_hh__spikespace[0]));
		outfile__array_hh__spikespace.close();
	} else
	{
		std::cout << "Error writing output file for _array_hh__spikespace." << endl;
	}
	ofstream outfile__array_hh_clock_dt;
	outfile__array_hh_clock_dt.open(results_dir + "_array_hh_clock_dt_150566293", ios::binary | ios::out);
	if(outfile__array_hh_clock_dt.is_open())
	{
		outfile__array_hh_clock_dt.write(reinterpret_cast<char*>(_array_hh_clock_dt), 1*sizeof(_array_hh_clock_dt[0]));
		outfile__array_hh_clock_dt.close();
	} else
	{
		std::cout << "Error writing output file for _array_hh_clock_dt." << endl;
	}
	ofstream outfile__array_hh_clock_t;
	outfile__array_hh_clock_t.open(results_dir + "_array_hh_clock_t_2257967227", ios::binary | ios::out);
	if(outfile__array_hh_clock_t.is_open())
	{
		outfile__array_hh_clock_t.write(reinterpret_cast<char*>(_array_hh_clock_t), 1*sizeof(_array_hh_clock_t[0]));
		outfile__array_hh_clock_t.close();
	} else
	{
		std::cout << "Error writing output file for _array_hh_clock_t." << endl;
	}
	ofstream outfile__array_hh_clock_timestep;
	outfile__array_hh_clock_timestep.open(results_dir + "_array_hh_clock_timestep_3318827205", ios::binary | ios::out);
	if(outfile__array_hh_clock_timestep.is_open())
	{
		outfile__array_hh_clock_timestep.write(reinterpret_cast<char*>(_array_hh_clock_timestep), 1*sizeof(_array_hh_clock_timestep[0]));
		outfile__array_hh_clock_timestep.close();
	} else
	{
		std::cout << "Error writing output file for _array_hh_clock_timestep." << endl;
	}
	ofstream outfile__array_hh_h;
	outfile__array_hh_h.open(results_dir + "_array_hh_h_4084671545", ios::binary | ios::out);
	if(outfile__array_hh_h.is_open())
	{
		outfile__array_hh_h.write(reinterpret_cast<char*>(_array_hh_h), 1*sizeof(_array_hh_h[0]));
		outfile__array_hh_h.close();
	} else
	{
		std::cout << "Error writing output file for _array_hh_h." << endl;
	}
	ofstream outfile__array_hh_i;
	outfile__array_hh_i.open(results_dir + "_array_hh_i_2221937839", ios::binary | ios::out);
	if(outfile__array_hh_i.is_open())
	{
		outfile__array_hh_i.write(reinterpret_cast<char*>(_array_hh_i), 1*sizeof(_array_hh_i[0]));
		outfile__array_hh_i.close();
	} else
	{
		std::cout << "Error writing output file for _array_hh_i." << endl;
	}
	ofstream outfile__array_hh_m;
	outfile__array_hh_m.open(results_dir + "_array_hh_m_2199769270", ios::binary | ios::out);
	if(outfile__array_hh_m.is_open())
	{
		outfile__array_hh_m.write(reinterpret_cast<char*>(_array_hh_m), 1*sizeof(_array_hh_m[0]));
		outfile__array_hh_m.close();
	} else
	{
		std::cout << "Error writing output file for _array_hh_m." << endl;
	}
	ofstream outfile__array_hh_n;
	outfile__array_hh_n.open(results_dir + "_array_hh_n_437551372", ios::binary | ios::out);
	if(outfile__array_hh_n.is_open())
	{
		outfile__array_hh_n.write(reinterpret_cast<char*>(_array_hh_n), 1*sizeof(_array_hh_n[0]));
		outfile__array_hh_n.close();
	} else
	{
		std::cout << "Error writing output file for _array_hh_n." << endl;
	}
	ofstream outfile__array_hh_n_Cl_E;
	outfile__array_hh_n_Cl_E.open(results_dir + "_array_hh_n_Cl_E_3029291836", ios::binary | ios::out);
	if(outfile__array_hh_n_Cl_E.is_open())
	{
		outfile__array_hh_n_Cl_E.write(reinterpret_cast<char*>(_array_hh_n_Cl_E), 1*sizeof(_array_hh_n_Cl_E[0]));
		outfile__array_hh_n_Cl_E.close();
	} else
	{
		std::cout << "Error writing output file for _array_hh_n_Cl_E." << endl;
	}
	ofstream outfile__array_hh_n_Cl_N;
	outfile__array_hh_n_Cl_N.open(results_dir + "_array_hh_n_Cl_N_593332916", ios::binary | ios::out);
	if(outfile__array_hh_n_Cl_N.is_open())
	{
		outfile__array_hh_n_Cl_N.write(reinterpret_cast<char*>(_array_hh_n_Cl_N), 1*sizeof(_array_hh_n_Cl_N[0]));
		outfile__array_hh_n_Cl_N.close();
	} else
	{
		std::cout << "Error writing output file for _array_hh_n_Cl_N." << endl;
	}
	ofstream outfile__array_hh_n_K_E;
	outfile__array_hh_n_K_E.open(results_dir + "_array_hh_n_K_E_2359146093", ios::binary | ios::out);
	if(outfile__array_hh_n_K_E.is_open())
	{
		outfile__array_hh_n_K_E.write(reinterpret_cast<char*>(_array_hh_n_K_E), 1*sizeof(_array_hh_n_K_E[0]));
		outfile__array_hh_n_K_E.close();
	} else
	{
		std::cout << "Error writing output file for _array_hh_n_K_E." << endl;
	}
	ofstream outfile__array_hh_n_K_N;
	outfile__array_hh_n_K_N.open(results_dir + "_array_hh_n_K_N_458190821", ios::binary | ios::out);
	if(outfile__array_hh_n_K_N.is_open())
	{
		outfile__array_hh_n_K_N.write(reinterpret_cast<char*>(_array_hh_n_K_N), 1*sizeof(_array_hh_n_K_N[0]));
		outfile__array_hh_n_K_N.close();
	} else
	{
		std::cout << "Error writing output file for _array_hh_n_K_N." << endl;
	}
	ofstream outfile__array_hh_n_Na_E;
	outfile__array_hh_n_Na_E.open(results_dir + "_array_hh_n_Na_E_1312626866", ios::binary | ios::out);
	if(outfile__array_hh_n_Na_E.is_open())
	{
		outfile__array_hh_n_Na_E.write(reinterpret_cast<char*>(_array_hh_n_Na_E), 1*sizeof(_array_hh_n_Na_E[0]));
		outfile__array_hh_n_Na_E.close();
	} else
	{
		std::cout << "Error writing output file for _array_hh_n_Na_E." << endl;
	}
	ofstream outfile__array_hh_n_Na_N;
	outfile__array_hh_n_Na_N.open(results_dir + "_array_hh_n_Na_N_3656368442", ios::binary | ios::out);
	if(outfile__array_hh_n_Na_N.is_open())
	{
		outfile__array_hh_n_Na_N.write(reinterpret_cast<char*>(_array_hh_n_Na_N), 1*sizeof(_array_hh_n_Na_N[0]));
		outfile__array_hh_n_Na_N.close();
	} else
	{
		std::cout << "Error writing output file for _array_hh_n_Na_N." << endl;
	}
	ofstream outfile__array_hh_v;
	outfile__array_hh_v.open(results_dir + "_array_hh_v_158865754", ios::binary | ios::out);
	if(outfile__array_hh_v.is_open())
	{
		outfile__array_hh_v.write(reinterpret_cast<char*>(_array_hh_v), 1*sizeof(_array_hh_v[0]));
		outfile__array_hh_v.close();
	} else
	{
		std::cout << "Error writing output file for _array_hh_v." << endl;
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
	ofstream outfile__array_svmon_clock_dt;
	outfile__array_svmon_clock_dt.open(results_dir + "_array_svmon_clock_dt_1737044899", ios::binary | ios::out);
	if(outfile__array_svmon_clock_dt.is_open())
	{
		outfile__array_svmon_clock_dt.write(reinterpret_cast<char*>(_array_svmon_clock_dt), 1*sizeof(_array_svmon_clock_dt[0]));
		outfile__array_svmon_clock_dt.close();
	} else
	{
		std::cout << "Error writing output file for _array_svmon_clock_dt." << endl;
	}
	ofstream outfile__array_svmon_clock_t;
	outfile__array_svmon_clock_t.open(results_dir + "_array_svmon_clock_t_1185157037", ios::binary | ios::out);
	if(outfile__array_svmon_clock_t.is_open())
	{
		outfile__array_svmon_clock_t.write(reinterpret_cast<char*>(_array_svmon_clock_t), 1*sizeof(_array_svmon_clock_t[0]));
		outfile__array_svmon_clock_t.close();
	} else
	{
		std::cout << "Error writing output file for _array_svmon_clock_t." << endl;
	}
	ofstream outfile__array_svmon_clock_timestep;
	outfile__array_svmon_clock_timestep.open(results_dir + "_array_svmon_clock_timestep_3908381981", ios::binary | ios::out);
	if(outfile__array_svmon_clock_timestep.is_open())
	{
		outfile__array_svmon_clock_timestep.write(reinterpret_cast<char*>(_array_svmon_clock_timestep), 1*sizeof(_array_svmon_clock_timestep[0]));
		outfile__array_svmon_clock_timestep.close();
	} else
	{
		std::cout << "Error writing output file for _array_svmon_clock_timestep." << endl;
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
	ofstream outfile__dynamic_array_svmon_I_Cl;
	outfile__dynamic_array_svmon_I_Cl.open(results_dir + "_dynamic_array_svmon_I_Cl_2433979663", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_I_Cl.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_I_Cl.n; n++)
        {
            if (! _dynamic_array_svmon_I_Cl(n).empty())
            {
                outfile__dynamic_array_svmon_I_Cl.write(reinterpret_cast<char*>(&_dynamic_array_svmon_I_Cl(n, 0)), _dynamic_array_svmon_I_Cl.m*sizeof(_dynamic_array_svmon_I_Cl(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_I_Cl.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_I_Cl." << endl;
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
}

