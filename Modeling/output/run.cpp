#include<stdlib.h>
#include "objects.h"
#include<ctime>
#include "randomkit.h"

#include "code_objects/NGexc_spike_resetter_codeobject_1.h"
#include "code_objects/NGexc_spike_thresholder_codeobject_1.h"
#include "code_objects/NGexc_stateupdater_codeobject_1.h"
#include "code_objects/NGinh_spike_resetter_codeobject_1.h"
#include "code_objects/NGinh_spike_thresholder_codeobject_1.h"
#include "code_objects/NGinh_stateupdater_codeobject_1.h"
#include "code_objects/spikemonitor_2_codeobject.h"
#include "code_objects/spikemonitor_3_codeobject.h"
#include "code_objects/svmon_exc_add_codeobject_1.h"
#include "code_objects/svmon_inh_add_codeobject_1.h"


void brian_start()
{
	_init_arrays();
	_load_arrays();
	// Initialize clocks (link timestep and dt to the respective arrays)
    brian::NGexc_clock_1.timestep = brian::_array_NGexc_clock_1_timestep;
    brian::NGexc_clock_1.dt = brian::_array_NGexc_clock_1_dt;
    brian::NGexc_clock_1.t = brian::_array_NGexc_clock_1_t;
    brian::NGinh_clock_1.timestep = brian::_array_NGinh_clock_1_timestep;
    brian::NGinh_clock_1.dt = brian::_array_NGinh_clock_1_dt;
    brian::NGinh_clock_1.t = brian::_array_NGinh_clock_1_t;
    brian::svmon_exc_add_clock_1.timestep = brian::_array_svmon_exc_add_clock_1_timestep;
    brian::svmon_exc_add_clock_1.dt = brian::_array_svmon_exc_add_clock_1_dt;
    brian::svmon_exc_add_clock_1.t = brian::_array_svmon_exc_add_clock_1_t;
    brian::svmon_inh_add_clock_1.timestep = brian::_array_svmon_inh_add_clock_1_timestep;
    brian::svmon_inh_add_clock_1.dt = brian::_array_svmon_inh_add_clock_1_dt;
    brian::svmon_inh_add_clock_1.t = brian::_array_svmon_inh_add_clock_1_t;
    for (int i=0; i<1; i++)
	    rk_randomseed(brian::_mersenne_twister_states[i]);  // Note that this seed can be potentially replaced in main.cpp
}

void brian_end()
{
	_write_arrays();
	_dealloc_arrays();
}


