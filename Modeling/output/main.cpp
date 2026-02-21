#include <stdlib.h>
#include "objects.h"
#include <ctime>
#include <time.h>
#include <omp.h>
#include "run.h"
#include "brianlib/common_math.h"
#include "randomkit.h"

#include "code_objects/NGexc_spike_resetter_codeobject_1.h"
#include "code_objects/NGexc_spike_thresholder_codeobject_1.h"
#include "code_objects/after_run_NGexc_spike_thresholder_codeobject_1.h"
#include "code_objects/NGexc_stateupdater_codeobject_1.h"
#include "code_objects/NGinh_spike_resetter_codeobject_1.h"
#include "code_objects/NGinh_spike_thresholder_codeobject_1.h"
#include "code_objects/after_run_NGinh_spike_thresholder_codeobject_1.h"
#include "code_objects/NGinh_stateupdater_codeobject_1.h"
#include "code_objects/spikemonitor_2_codeobject.h"
#include "code_objects/spikemonitor_3_codeobject.h"
#include "code_objects/svmon_exc_add_codeobject_1.h"
#include "code_objects/svmon_inh_add_codeobject_1.h"


#include <iostream>
#include <fstream>
#include <string>


        std::string _format_time(float time_in_s)
        {
            float divisors[] = {24*60*60, 60*60, 60, 1};
            char letters[] = {'d', 'h', 'm', 's'};
            float remaining = time_in_s;
            std::string text = "";
            int time_to_represent;
            for (int i =0; i < sizeof(divisors)/sizeof(float); i++)
            {
                time_to_represent = int(remaining / divisors[i]);
                remaining -= time_to_represent * divisors[i];
                if (time_to_represent > 0 || text.length())
                {
                    if(text.length() > 0)
                    {
                        text += " ";
                    }
                    text += (std::to_string(time_to_represent)+letters[i]);
                }
            }
            //less than one second
            if(text.length() == 0)
            {
                text = "< 1s";
            }
            return text;
        }
        void report_progress(const double elapsed, const double completed, const double start, const double duration)
        {
            if (completed == 0.0)
            {
                std::cout << "Starting simulation at t=" << start << " s for duration " << duration << " s";
            } else
            {
                std::cout << completed*duration << " s (" << (int)(completed*100.) << "%) simulated in " << _format_time(elapsed);
                if (completed < 1.0)
                {
                    const int remaining = (int)((1-completed)/completed*elapsed+0.5);
                    std::cout << ", estimated " << _format_time(remaining) << " remaining.";
                }
            }

            std::cout << std::endl << std::flush;
        }
        


void set_from_command_line(const std::vector<std::string> args)
{
    for (const auto& arg : args) {
		// Split into two parts
		size_t equal_sign = arg.find("=");
		auto name = arg.substr(0, equal_sign);
		auto value = arg.substr(equal_sign + 1, arg.length());
		brian::set_variable_by_name(name, value);
	}
}
int main(int argc, char **argv)
{
	std::vector<std::string> args(argv + 1, argv + argc);
	if (args.size() >=2 && args[0] == "--results_dir")
	{
		brian::results_dir = args[1];
		#ifdef DEBUG
		std::cout << "Setting results dir to '" << brian::results_dir << "'" << std::endl;
		#endif
		args.erase(args.begin(), args.begin()+2);
	}
        

	brian_start();
        

	{
		using namespace brian;

		omp_set_dynamic(0);
omp_set_num_threads(1);
                
        _array_defaultclock_dt[0] = 0.0001;
        _array_defaultclock_dt[0] = 0.0001;
        _array_NGexc_clock_1_dt[0] = 1e-07;
        _array_NGexc_clock_1_dt[0] = 1e-07;
        _array_NGexc_lastspike[0] = - 10000.0;
        _array_NGexc_not_refractory[0] = true;
        _array_NGexc_v[0] = - 0.0725;
        _array_NGexc_m[0] = 0.0003919714799822308;
        _array_NGexc_h[0] = 0.9999991878248916;
        _array_NGexc_n[0] = 0.001983630308690279;
        _array_NGexc_p[0] = 0.022977369910025615;
        _array_NGexc_I_max[0] = 0.015;
        _array_NGexc_mask_noise[0] = true;
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__static_array__index__array_NGexc_x; i++)
                        {
                            _array_NGexc_x[_static_array__index__array_NGexc_x[i]] = _static_array__value__array_NGexc_x[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__static_array__index__array_NGexc_y; i++)
                        {
                            _array_NGexc_y[_static_array__index__array_NGexc_y[i]] = _static_array__value__array_NGexc_y[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__static_array__index__array_NGexc_z; i++)
                        {
                            _array_NGexc_z[_static_array__index__array_NGexc_z[i]] = _static_array__value__array_NGexc_z[i];
                        }
                        
        _array_NGinh_clock_1_dt[0] = 1e-07;
        _array_NGinh_clock_1_dt[0] = 1e-07;
        _array_NGinh_lastspike[0] = - 10000.0;
        _array_NGinh_not_refractory[0] = true;
        _array_NGinh_v[0] = - 0.0648;
        _array_NGinh_m[0] = 0.0022991085011206743;
        _array_NGinh_h[0] = 0.9999941893858119;
        _array_NGinh_n[0] = 0.008462011693223123;
        _array_NGinh_p[0] = 0.04833762871990523;
        _array_NGinh_I_max[0] = 0.015;
        _array_NGinh_mask_noise[0] = true;
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__static_array__index__array_NGinh_x; i++)
                        {
                            _array_NGinh_x[_static_array__index__array_NGinh_x[i]] = _static_array__value__array_NGinh_x[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__static_array__index__array_NGinh_y; i++)
                        {
                            _array_NGinh_y[_static_array__index__array_NGinh_y[i]] = _static_array__value__array_NGinh_y[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__static_array__index__array_NGinh_z; i++)
                        {
                            _array_NGinh_z[_static_array__index__array_NGinh_z[i]] = _static_array__value__array_NGinh_z[i];
                        }
                        
        _array_svmon_exc_add_clock_1_dt[0] = 1e-05;
        _array_svmon_exc_add_clock_1_dt[0] = 1e-05;
        _array_svmon_exc_add__indices[0] = false;
        _array_svmon_inh_add_clock_1_dt[0] = 1e-05;
        _array_svmon_inh_add_clock_1_dt[0] = 1e-05;
        _array_svmon_inh_add__indices[0] = false;
        _array_NGexc_clock_1_timestep[0] = false;
        _array_NGexc_clock_1_t[0] = false;
        _array_svmon_inh_add_clock_1_timestep[0] = false;
        _array_svmon_inh_add_clock_1_t[0] = false;
        _array_svmon_exc_add_clock_1_timestep[0] = false;
        _array_svmon_exc_add_clock_1_t[0] = false;
        _array_NGinh_clock_1_timestep[0] = false;
        _array_NGinh_clock_1_t[0] = false;
        network_1.clear();
        network_1.add(&svmon_exc_add_clock_1, _run_svmon_exc_add_codeobject_1);
        network_1.add(&svmon_inh_add_clock_1, _run_svmon_inh_add_codeobject_1);
        network_1.add(&NGexc_clock_1, _run_NGexc_stateupdater_codeobject_1);
        network_1.add(&NGinh_clock_1, _run_NGinh_stateupdater_codeobject_1);
        network_1.add(&NGexc_clock_1, _run_NGexc_spike_thresholder_codeobject_1);
        network_1.add(&NGinh_clock_1, _run_NGinh_spike_thresholder_codeobject_1);
        network_1.add(&NGexc_clock_1, _run_spikemonitor_2_codeobject);
        network_1.add(&NGinh_clock_1, _run_spikemonitor_3_codeobject);
        network_1.add(&NGexc_clock_1, _run_NGexc_spike_resetter_codeobject_1);
        network_1.add(&NGinh_clock_1, _run_NGinh_spike_resetter_codeobject_1);
        set_from_command_line(args);
        network_1.run(2.0, report_progress, 10.0);
        _after_run_NGexc_spike_thresholder_codeobject_1();
        _after_run_NGinh_spike_thresholder_codeobject_1();
        #ifdef DEBUG
        _debugmsg_spikemonitor_2_codeobject();
        #endif
        
        #ifdef DEBUG
        _debugmsg_spikemonitor_3_codeobject();
        #endif

	}
        

	brian_end();
        

	return 0;
}