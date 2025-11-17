

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
Network network_3;

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
    if (name == "hh_4._spikespace") {
        var_size = 2;
        data_size = 2*sizeof(int32_t);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<int32_t>(name, _array_hh_4__spikespace, var_size, (int32_t)atoi(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_hh_4__spikespace, data_size, s_value);
        }
        return;
    }
    if (name == "hh_4.h") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_hh_4_h, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_hh_4_h, data_size, s_value);
        }
        return;
    }
    if (name == "hh_4.m") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_hh_4_m, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_hh_4_m, data_size, s_value);
        }
        return;
    }
    if (name == "hh_4.n") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_hh_4_n, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_hh_4_n, data_size, s_value);
        }
        return;
    }
    if (name == "hh_4.v") {
        var_size = 1;
        data_size = 1*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_hh_4_v, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_hh_4_v, data_size, s_value);
        }
        return;
    }
    // dynamic arrays (1d)
    if (name == "_timedarray_5.values") {
        var_size = 4;
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _timedarray_5_values, var_size, (double)atof(s_value.c_str()));


        } else {
            // set from file
            set_variable_from_file(name, _timedarray_5_values, data_size, s_value);
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
int32_t * _array_hh_4__spikespace;
const int _num__array_hh_4__spikespace = 2;
double * _array_hh_4_clock_dt;
const int _num__array_hh_4_clock_dt = 1;
double * _array_hh_4_clock_t;
const int _num__array_hh_4_clock_t = 1;
int64_t * _array_hh_4_clock_timestep;
const int _num__array_hh_4_clock_timestep = 1;
double * _array_hh_4_h;
const int _num__array_hh_4_h = 1;
int32_t * _array_hh_4_i;
const int _num__array_hh_4_i = 1;
double * _array_hh_4_m;
const int _num__array_hh_4_m = 1;
double * _array_hh_4_n;
const int _num__array_hh_4_n = 1;
double * _array_hh_4_v;
const int _num__array_hh_4_v = 1;
int32_t * _array_svmon__indices;
const int _num__array_svmon__indices = 1;
double * _array_svmon_clock_3_dt;
const int _num__array_svmon_clock_3_dt = 1;
double * _array_svmon_clock_3_t;
const int _num__array_svmon_clock_3_t = 1;
int64_t * _array_svmon_clock_3_timestep;
const int _num__array_svmon_clock_3_timestep = 1;
double * _array_svmon_h;
const int _num__array_svmon_h = (0, 1);
double * _array_svmon_m;
const int _num__array_svmon_m = (0, 1);
int32_t * _array_svmon_N;
const int _num__array_svmon_N = 1;
double * _array_svmon_n;
const int _num__array_svmon_n = (0, 1);
double * _array_svmon_v;
const int _num__array_svmon_v = (0, 1);

//////////////// dynamic arrays 1d /////////
std::vector<double> _dynamic_array_svmon_t;

//////////////// dynamic arrays 2d /////////
DynamicArray2D<double> _dynamic_array_svmon_h;
DynamicArray2D<double> _dynamic_array_svmon_m;
DynamicArray2D<double> _dynamic_array_svmon_n;
DynamicArray2D<double> _dynamic_array_svmon_v;

/////////////// static arrays /////////////
double * _timedarray_5_values;
const int _num__timedarray_5_values = 4;

//////////////// synapses /////////////////

//////////////// clocks ///////////////////
Clock hh_4_clock;  // attributes will be set in run.cpp
Clock svmon_clock_3;  // attributes will be set in run.cpp

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

	_array_hh_4__spikespace = new int32_t[2];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<2; i++) _array_hh_4__spikespace[i] = 0;

	_array_hh_4_clock_dt = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_hh_4_clock_dt[i] = 0;

	_array_hh_4_clock_t = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_hh_4_clock_t[i] = 0;

	_array_hh_4_clock_timestep = new int64_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_hh_4_clock_timestep[i] = 0;

	_array_hh_4_h = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_hh_4_h[i] = 0;

	_array_hh_4_i = new int32_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_hh_4_i[i] = 0;

	_array_hh_4_m = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_hh_4_m[i] = 0;

	_array_hh_4_n = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_hh_4_n[i] = 0;

	_array_hh_4_v = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_hh_4_v[i] = 0;

	_array_svmon__indices = new int32_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_svmon__indices[i] = 0;

	_array_svmon_clock_3_dt = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_svmon_clock_3_dt[i] = 0;

	_array_svmon_clock_3_t = new double[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_svmon_clock_3_t[i] = 0;

	_array_svmon_clock_3_timestep = new int64_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_svmon_clock_3_timestep[i] = 0;

	_array_svmon_N = new int32_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_svmon_N[i] = 0;


	// Arrays initialized to an "arange"
	_array_hh_4_i = new int32_t[1];
    #pragma omp parallel for schedule(static)
	for(int i=0; i<1; i++) _array_hh_4_i[i] = 0 + i;


	// static arrays
	_timedarray_5_values = new double[4];

	// Random number generator states
	for (int i=0; i<1; i++)
	    _mersenne_twister_states.push_back(new rk_state());
}

void _load_arrays()
{
	using namespace brian;

	ifstream f_timedarray_5_values;
	f_timedarray_5_values.open("static_arrays/_timedarray_5_values", ios::in | ios::binary);
	if(f_timedarray_5_values.is_open())
	{
		f_timedarray_5_values.read(reinterpret_cast<char*>(_timedarray_5_values), 4*sizeof(double));
	} else
	{
		std::cout << "Error opening static array _timedarray_5_values." << endl;
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
	ofstream outfile__array_hh_4__spikespace;
	outfile__array_hh_4__spikespace.open(results_dir + "_array_hh_4__spikespace_1439089419", ios::binary | ios::out);
	if(outfile__array_hh_4__spikespace.is_open())
	{
		outfile__array_hh_4__spikespace.write(reinterpret_cast<char*>(_array_hh_4__spikespace), 2*sizeof(_array_hh_4__spikespace[0]));
		outfile__array_hh_4__spikespace.close();
	} else
	{
		std::cout << "Error writing output file for _array_hh_4__spikespace." << endl;
	}
	ofstream outfile__array_hh_4_clock_dt;
	outfile__array_hh_4_clock_dt.open(results_dir + "_array_hh_4_clock_dt_1171345450", ios::binary | ios::out);
	if(outfile__array_hh_4_clock_dt.is_open())
	{
		outfile__array_hh_4_clock_dt.write(reinterpret_cast<char*>(_array_hh_4_clock_dt), 1*sizeof(_array_hh_4_clock_dt[0]));
		outfile__array_hh_4_clock_dt.close();
	} else
	{
		std::cout << "Error writing output file for _array_hh_4_clock_dt." << endl;
	}
	ofstream outfile__array_hh_4_clock_t;
	outfile__array_hh_4_clock_t.open(results_dir + "_array_hh_4_clock_t_483475227", ios::binary | ios::out);
	if(outfile__array_hh_4_clock_t.is_open())
	{
		outfile__array_hh_4_clock_t.write(reinterpret_cast<char*>(_array_hh_4_clock_t), 1*sizeof(_array_hh_4_clock_t[0]));
		outfile__array_hh_4_clock_t.close();
	} else
	{
		std::cout << "Error writing output file for _array_hh_4_clock_t." << endl;
	}
	ofstream outfile__array_hh_4_clock_timestep;
	outfile__array_hh_4_clock_timestep.open(results_dir + "_array_hh_4_clock_timestep_984141907", ios::binary | ios::out);
	if(outfile__array_hh_4_clock_timestep.is_open())
	{
		outfile__array_hh_4_clock_timestep.write(reinterpret_cast<char*>(_array_hh_4_clock_timestep), 1*sizeof(_array_hh_4_clock_timestep[0]));
		outfile__array_hh_4_clock_timestep.close();
	} else
	{
		std::cout << "Error writing output file for _array_hh_4_clock_timestep." << endl;
	}
	ofstream outfile__array_hh_4_h;
	outfile__array_hh_4_h.open(results_dir + "_array_hh_4_h_1756559302", ios::binary | ios::out);
	if(outfile__array_hh_4_h.is_open())
	{
		outfile__array_hh_4_h.write(reinterpret_cast<char*>(_array_hh_4_h), 1*sizeof(_array_hh_4_h[0]));
		outfile__array_hh_4_h.close();
	} else
	{
		std::cout << "Error writing output file for _array_hh_4_h." << endl;
	}
	ofstream outfile__array_hh_4_i;
	outfile__array_hh_4_i.open(results_dir + "_array_hh_4_i_532006736", ios::binary | ios::out);
	if(outfile__array_hh_4_i.is_open())
	{
		outfile__array_hh_4_i.write(reinterpret_cast<char*>(_array_hh_4_i), 1*sizeof(_array_hh_4_i[0]));
		outfile__array_hh_4_i.close();
	} else
	{
		std::cout << "Error writing output file for _array_hh_4_i." << endl;
	}
	ofstream outfile__array_hh_4_m;
	outfile__array_hh_4_m.open(results_dir + "_array_hh_4_m_416809801", ios::binary | ios::out);
	if(outfile__array_hh_4_m.is_open())
	{
		outfile__array_hh_4_m.write(reinterpret_cast<char*>(_array_hh_4_m), 1*sizeof(_array_hh_4_m[0]));
		outfile__array_hh_4_m.close();
	} else
	{
		std::cout << "Error writing output file for _array_hh_4_m." << endl;
	}
	ofstream outfile__array_hh_4_n;
	outfile__array_hh_4_n.open(results_dir + "_array_hh_4_n_2177979123", ios::binary | ios::out);
	if(outfile__array_hh_4_n.is_open())
	{
		outfile__array_hh_4_n.write(reinterpret_cast<char*>(_array_hh_4_n), 1*sizeof(_array_hh_4_n[0]));
		outfile__array_hh_4_n.close();
	} else
	{
		std::cout << "Error writing output file for _array_hh_4_n." << endl;
	}
	ofstream outfile__array_hh_4_v;
	outfile__array_hh_4_v.open(results_dir + "_array_hh_4_v_2461911717", ios::binary | ios::out);
	if(outfile__array_hh_4_v.is_open())
	{
		outfile__array_hh_4_v.write(reinterpret_cast<char*>(_array_hh_4_v), 1*sizeof(_array_hh_4_v[0]));
		outfile__array_hh_4_v.close();
	} else
	{
		std::cout << "Error writing output file for _array_hh_4_v." << endl;
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
	ofstream outfile__array_svmon_clock_3_dt;
	outfile__array_svmon_clock_3_dt.open(results_dir + "_array_svmon_clock_3_dt_1327005398", ios::binary | ios::out);
	if(outfile__array_svmon_clock_3_dt.is_open())
	{
		outfile__array_svmon_clock_3_dt.write(reinterpret_cast<char*>(_array_svmon_clock_3_dt), 1*sizeof(_array_svmon_clock_3_dt[0]));
		outfile__array_svmon_clock_3_dt.close();
	} else
	{
		std::cout << "Error writing output file for _array_svmon_clock_3_dt." << endl;
	}
	ofstream outfile__array_svmon_clock_3_t;
	outfile__array_svmon_clock_3_t.open(results_dir + "_array_svmon_clock_3_t_3578913941", ios::binary | ios::out);
	if(outfile__array_svmon_clock_3_t.is_open())
	{
		outfile__array_svmon_clock_3_t.write(reinterpret_cast<char*>(_array_svmon_clock_3_t), 1*sizeof(_array_svmon_clock_3_t[0]));
		outfile__array_svmon_clock_3_t.close();
	} else
	{
		std::cout << "Error writing output file for _array_svmon_clock_3_t." << endl;
	}
	ofstream outfile__array_svmon_clock_3_timestep;
	outfile__array_svmon_clock_3_timestep.open(results_dir + "_array_svmon_clock_3_timestep_551271438", ios::binary | ios::out);
	if(outfile__array_svmon_clock_3_timestep.is_open())
	{
		outfile__array_svmon_clock_3_timestep.write(reinterpret_cast<char*>(_array_svmon_clock_3_timestep), 1*sizeof(_array_svmon_clock_3_timestep[0]));
		outfile__array_svmon_clock_3_timestep.close();
	} else
	{
		std::cout << "Error writing output file for _array_svmon_clock_3_timestep." << endl;
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

	ofstream outfile__dynamic_array_svmon_h;
	outfile__dynamic_array_svmon_h.open(results_dir + "_dynamic_array_svmon_h_1813511655", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_h.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_h.n; n++)
        {
            if (! _dynamic_array_svmon_h(n).empty())
            {
                outfile__dynamic_array_svmon_h.write(reinterpret_cast<char*>(&_dynamic_array_svmon_h(n, 0)), _dynamic_array_svmon_h.m*sizeof(_dynamic_array_svmon_h(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_h.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_h." << endl;
	}
	ofstream outfile__dynamic_array_svmon_m;
	outfile__dynamic_array_svmon_m.open(results_dir + "_dynamic_array_svmon_m_477956456", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_m.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_m.n; n++)
        {
            if (! _dynamic_array_svmon_m(n).empty())
            {
                outfile__dynamic_array_svmon_m.write(reinterpret_cast<char*>(&_dynamic_array_svmon_m(n, 0)), _dynamic_array_svmon_m.m*sizeof(_dynamic_array_svmon_m(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_m.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_m." << endl;
	}
	ofstream outfile__dynamic_array_svmon_n;
	outfile__dynamic_array_svmon_n.open(results_dir + "_dynamic_array_svmon_n_2238994642", ios::binary | ios::out);
	if(outfile__dynamic_array_svmon_n.is_open())
	{
        for (int n=0; n<_dynamic_array_svmon_n.n; n++)
        {
            if (! _dynamic_array_svmon_n(n).empty())
            {
                outfile__dynamic_array_svmon_n.write(reinterpret_cast<char*>(&_dynamic_array_svmon_n(n, 0)), _dynamic_array_svmon_n.m*sizeof(_dynamic_array_svmon_n(0, 0)));
            }
        }
        outfile__dynamic_array_svmon_n.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_svmon_n." << endl;
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
	if(_timedarray_5_values!=0)
	{
		delete [] _timedarray_5_values;
		_timedarray_5_values = 0;
	}
}

