PROGRAM = main

SRCS = brianlib/randomkit/randomkit.c code_objects/NGexc_spike_resetter_codeobject_1.cpp code_objects/NGexc_spike_thresholder_codeobject_1.cpp code_objects/NGexc_stateupdater_codeobject_1.cpp code_objects/NGinh_spike_resetter_codeobject_1.cpp code_objects/NGinh_spike_thresholder_codeobject_1.cpp code_objects/NGinh_stateupdater_codeobject_1.cpp code_objects/after_run_NGexc_spike_thresholder_codeobject_1.cpp code_objects/after_run_NGinh_spike_thresholder_codeobject_1.cpp code_objects/spikemonitor_2_codeobject.cpp code_objects/spikemonitor_3_codeobject.cpp code_objects/svmon_exc_add_codeobject_1.cpp code_objects/svmon_inh_add_codeobject_1.cpp main.cpp network.cpp objects.cpp run.cpp synapses_classes.cpp
H_SRCS = brianlib/clocks.h brianlib/common_math.h brianlib/dynamic_array.h code_objects/NGexc_spike_resetter_codeobject_1.h code_objects/NGexc_spike_thresholder_codeobject_1.h code_objects/NGexc_stateupdater_codeobject_1.h code_objects/NGinh_spike_resetter_codeobject_1.h code_objects/NGinh_spike_thresholder_codeobject_1.h code_objects/NGinh_stateupdater_codeobject_1.h code_objects/after_run_NGexc_spike_thresholder_codeobject_1.h code_objects/after_run_NGinh_spike_thresholder_codeobject_1.h code_objects/spikemonitor_2_codeobject.h code_objects/spikemonitor_3_codeobject.h code_objects/svmon_exc_add_codeobject_1.h code_objects/svmon_inh_add_codeobject_1.h network.h objects.h run.h synapses_classes.h
OBJS = ${SRCS:.cpp=.o}
OBJS := ${OBJS:.c=.o}
OPTIMISATIONS = -Ibrianlib/randomkit -I/opt/anaconda3/envs/Brian2_NGILab/include -w -O3 -ffast-math -fno-finite-math-only -march=native -std=c++11
CXXFLAGS = -c -Wno-write-strings $(OPTIMISATIONS) -I. -fopenmp 
LFLAGS = -fopenmp -Lbrianlib/randomkit -L/opt/anaconda3/envs/Brian2_NGILab/lib -Wl,-rpath,/opt/anaconda3/envs/Brian2_NGILab/lib 
DEPS = make.deps

all: $(PROGRAM)

.PHONY: all clean

$(PROGRAM): $(OBJS) $(DEPS) makefile
	$(CXX) $(OBJS) -o $(PROGRAM) $(LFLAGS)

clean:
	rm $(OBJS) $(PROGRAM) $(DEPS)

make.deps: $(SRCS) $(H_SRCS)
	$(CXX) $(CXXFLAGS) -MM $(SRCS) > make.deps
	
ifneq ($(wildcard $(DEPS)), )
include $(DEPS)
endif

%.o : %.cpp makefile
	$(CXX) $(CXXFLAGS) $< -o $@