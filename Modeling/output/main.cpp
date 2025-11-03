#include <stdlib.h>
#include "objects.h"
#include <ctime>
#include <time.h>
#include <omp.h>
#include "run.h"
#include "brianlib/common_math.h"
#include "randomkit.h"

#include "code_objects/hh_spike_thresholder_codeobject.h"
#include "code_objects/after_run_hh_spike_thresholder_codeobject.h"
#include "code_objects/hh_stateupdater_codeobject.h"
#include "code_objects/svmon_codeobject.h"


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
        _array_hh_clock_dt[0] = 1e-07;
        _array_hh_clock_dt[0] = 1e-07;
        _array_hh_m[0] = 0.01;
        _array_hh_h[0] = 0.99;
        _array_hh_n[0] = 0.01;
        _array_hh_v[0] = - 0.08707803552883991;
        _array_svmon_clock_dt[0] = 0.0001;
        _array_svmon_clock_dt[0] = 0.0001;
        _array_svmon__indices[0] = false;
        _array_hh_clock_timestep[0] = false;
        _array_hh_clock_t[0] = false;
        _array_svmon_clock_timestep[0] = false;
        _array_svmon_clock_t[0] = false;
        network.clear();
        network.add(&svmon_clock, _run_svmon_codeobject);
        network.add(&hh_clock, _run_hh_stateupdater_codeobject);
        network.add(&hh_clock, _run_hh_spike_thresholder_codeobject);
        set_from_command_line(args);
        network.run(1.0, report_progress, 10.0);
        _after_run_hh_spike_thresholder_codeobject();

	}
        

	brian_end();
        

	return 0;
}