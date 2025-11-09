#include <stdlib.h>
#include "objects.h"
#include <ctime>
#include <time.h>
#include <omp.h>
#include "run.h"
#include "brianlib/common_math.h"
#include "randomkit.h"

#include "code_objects/hh_1_spike_resetter_codeobject.h"
#include "code_objects/hh_1_spike_thresholder_codeobject.h"
#include "code_objects/after_run_hh_1_spike_thresholder_codeobject.h"
#include "code_objects/hh_1_stateupdater_codeobject.h"
#include "code_objects/svmon_codeobject_1.h"


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
omp_set_num_threads(2);
                
        _array_defaultclock_dt[0] = 0.0001;
        _array_defaultclock_dt[0] = 0.0001;
        _array_hh_1_clock_dt[0] = 5e-08;
        _array_hh_1_clock_dt[0] = 5e-08;
        _array_hh_1_m[0] = 0.007070226997565917;
        _array_hh_1_h[0] = 0.6607563687658172;
        _array_hh_1_n[0] = 0.00010681238059392482;
        _array_hh_1_v[0] = - 0.07;
        _array_svmon_clock_1_dt[0] = 0.0002;
        _array_svmon_clock_1_dt[0] = 0.0002;
        _array_svmon__indices[0] = false;
        _array_svmon_clock_1_timestep[0] = false;
        _array_svmon_clock_1_t[0] = false;
        _array_hh_1_clock_timestep[0] = false;
        _array_hh_1_clock_t[0] = false;
        network_1.clear();
        network_1.add(&svmon_clock_1, _run_svmon_codeobject_1);
        network_1.add(&hh_1_clock, _run_hh_1_stateupdater_codeobject);
        network_1.add(&hh_1_clock, _run_hh_1_spike_thresholder_codeobject);
        network_1.add(&hh_1_clock, _run_hh_1_spike_resetter_codeobject);
        set_from_command_line(args);
        network_1.run(0.2, report_progress, 10.0);
        _after_run_hh_1_spike_thresholder_codeobject();

	}
        

	brian_end();
        

	return 0;
}