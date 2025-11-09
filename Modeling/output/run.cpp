#include<stdlib.h>
#include "objects.h"
#include<ctime>
#include "randomkit.h"

#include "code_objects/hh_1_spike_resetter_codeobject.h"
#include "code_objects/hh_1_spike_thresholder_codeobject.h"
#include "code_objects/hh_1_stateupdater_codeobject.h"
#include "code_objects/svmon_codeobject_1.h"


void brian_start()
{
	_init_arrays();
	_load_arrays();
	// Initialize clocks (link timestep and dt to the respective arrays)
    brian::hh_1_clock.timestep = brian::_array_hh_1_clock_timestep;
    brian::hh_1_clock.dt = brian::_array_hh_1_clock_dt;
    brian::hh_1_clock.t = brian::_array_hh_1_clock_t;
    brian::svmon_clock_1.timestep = brian::_array_svmon_clock_1_timestep;
    brian::svmon_clock_1.dt = brian::_array_svmon_clock_1_dt;
    brian::svmon_clock_1.t = brian::_array_svmon_clock_1_t;
    for (int i=0; i<2; i++)
	    rk_randomseed(brian::_mersenne_twister_states[i]);  // Note that this seed can be potentially replaced in main.cpp
}

void brian_end()
{
	_write_arrays();
	_dealloc_arrays();
}


