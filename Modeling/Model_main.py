from brian2 import *
defaultclock.dt = 0.01*ms
codegen.target='cython'
start_scope()

device.reinit()
device.activate()

#-----------------------------------------------------------------------------------------------------------------------
# Default Imports
#-----------------------------------------------------------------------------------------------------------------------
import numpy as np
import math
import scipy.constants as spc
import importlib
from brian2tools import *
from importnb import Notebook
from numpy import random

#-----------------------------------------------------------------------------------------------------------------------
# Brian2 import: we use Brian CPP-standalone code generation for fast parallelized simulations
#-----------------------------------------------------------------------------------------------------------------------
from brian2 import *
code_dir = './codegen'
prefs.GSL.directory = '/opt/anaconda3/envs/Brian2_NGILab/include/'   ## The directory where the GSL library headings are found
set_device('cpp_standalone',directory=code_dir,build_on_run=False)
prefs.devices.cpp_standalone.openmp_threads = 1 ## The number of threads used in the parallelization (machine-dependent)
prefs.logging.file_log = False
prefs.logging.delete_log_on_exit = True


import matplotlib.pyplot as plt
import Equations_calc
import parameters
importlib.reload(Equations_calc)
importlib.reload(parameters)
# %matplotlib widget
from Equations_calc import return_HH_equations, return_plotting_list
with Notebook():
    import Calculator as calc
# from Equations_2 import calc_leakage_conductance, I_Na_inf, f_NaK, sigma, calculate_I_NKP_max
from parameters import return_initial_parameters, return_synapse_pars


#-----------------------------------------------------------------------------------------------------------------------
## Utilities
#-----------------------------------------------------------------------------------------------------------------------
def varargin(pars, **kwargs):
    """
    varargin-like option for user-defined parameters in any function/module
    Use:
    pars = varargin(pars,**kwargs)

    Input:
    - pars     : the dictionary of parameters of the calling function
    - **kwargs : a dictionary of user-defined parameters

    Output:
    - pars     : modified dictionary of parameters to be used inside the calling
                 (parent) function

    Maurizio De Pitta', The University of Chicago, August 27th, 2014.
    """
    for key, val in kwargs.items():
        if key in pars:
            pars[key] = val
    return pars

def visualise_connectivity(S):
    Ns = len(S.source)
    Nt = len(S.target)
    figure(figsize=(10, 4))
    subplot(121)
    plot(zeros(Ns), arange(Ns), 'ok', ms=10)
    plot(ones(Nt), arange(Nt), 'ok', ms=10)
    for i, j in zip(S.i, S.j):
        plot([0, 1], [i, j], '-k')
    xticks([0, 1], ['Source', 'Target'])
    ylabel('Neuron index')
    xlim(-0.1, 1.1)
    ylim(-1, max(Ns, Nt))
    subplot(122)
    plot(S.i, S.j, 'ok')
    xlim(-1, Ns)
    ylim(-1, Nt)
    xlabel('Source neuron index')
    ylabel('Target neuron index')

#-----------------------------------------------------------------------------------------------------------------------
## Build User-defined convenience functions to be also called by equations in neuron models by Brian  
#-----------------------------------------------------------------------------------------------------------------------
def Hill(x,K,n):
    return x**n/(x**n+K**n)
Hill = Function(Hill,arg_units=[mmolar,mmolar,1], return_unit=1,auto_vectorise=False)
Hill_cpp = '''
    #include <math.h>
    double Hill(double x,double K,double n)
    {
        return pow(x,n)/(pow(x,n)+pow(K,n));
    };
    '''
Hill.implementations.add_implementation('cpp',Hill_cpp,compiler_kwds={'headers': ['"math.h"']})

def ThermalVoltage(T):
    return 8.314472000000000e+00*(T+273.15)/9.648534297750000e+04
ThermalPotential = Function(ThermalVoltage,arg_units=[1], return_unit=1,auto_vectorise=False)
ThermalVoltage_cpp = '''
    #include <gsl/gsl_const_mksa.h>
    double ThermalPotential(const double T)
    {
        const double R = GSL_CONST_MKSA_MOLAR_GAS;
        const double F = GSL_CONST_MKSA_FARADAY;
        return R*(T+273.15)/F;
    }
    '''
ThermalPotential.implementations.add_implementation('cpp',ThermalVoltage_cpp,
                                                  compiler_kwds={'headers': ['"gsl_const_mksa.h"'],
                                                                 'include_dirs': ['/opt/anaconda3/envs/Brian2_NGILab/include/gsl/']})

def NernstPotential(x_e,x_i,z,T):
    """
    Nernst potential in volts (w/out units)

    Input parameters (w/out units):
    - x_e : float   Intracellular concentration
    - x_i : float   Extracellyular concentration
    - z   : int     Ion valence
    - T   : float   Temperature in ^C

    Return:
    - E_x : Nernst Reverse Potential in volt (W/OUT units)
    """
    V_T = ThermalVoltage(T)
    return V_T/z*np.log(x_e/x_i)
NernstPotential = Function(NernstPotential,arg_units=[mmolar,mmolar,1,1], return_unit=1,auto_vectorise=False)
NernstPotential_cpp = '''
    #include <gsl/gsl_const_mksa.h>
    const double R = GSL_CONST_MKSA_MOLAR_GAS; 
    const double F = GSL_CONST_MKSA_FARADAY;
    double ThermalVoltage(const double T)
    {
        return R*(T+273.15)/F;
    };
    double NernstPotential(double x_e,double x_i,double z,const double T)
    {
        return ThermalVoltage(T)*log(x_e/x_i)/z;
    };
    '''
NernstPotential.implementations.add_implementation('cpp',NernstPotential_cpp,
                                                   dependencies={'log': DEFAULT_FUNCTIONS['log']},
                                                   compiler_kwds={'headers': ['"gsl_const_mksa.h"'],
                                                                  'include_dirs': ['/opt/anaconda3/envs/Brian2_NGILab/include/gsl/']})

def Heaviside_Function(t_step, t_end, t):
    if (t < t_step):
        return 0
    elif (t_step <= t <= t_end):
        return 1
    else:
        return 0
Heaviside_Function = Function(Heaviside_Function,arg_units=[second,second,second], return_unit=1,auto_vectorise=False)
Heaviside_Function_cpp = '''
    double Heaviside_Function(double t_step, double t_end, double t)
    {
    if (t < t_step)
        return 0;
    else if (t <= t_end)
        return 1;
    else
        return 0;
    };
    '''
Heaviside_Function.implementations.add_implementation('cpp',Heaviside_Function_cpp,
                                                   dependencies={'log': DEFAULT_FUNCTIONS['log']},
                                                    )

#-----------------------------------------------------------------------------------------------------------------------
## PARAMETER preparation 
#-----------------------------------------------------------------------------------------------------------------------
def lpc5_parameters(model='hh-neuron',**kwargs):
    ## Get parameters from parameters.py
    pars = return_initial_parameters()

    ## Conductance sodium leakage channel
    pars['g_Na_L'] = calc.calc_leakage_conductance()

    ## Adding I_NKP_max
    I_Na_inf_calc = calc.I_Na_inf(calc.calc_leakage_conductance(), pars['resting_potential']*mvolt)
    f_NaK_calc = calc.f_NaK(pars['resting_potential']*mvolt, ThermalVoltage(pars['T_exp'])*1000*mvolt)
    Hill_K = Hill(pars['C0_K_E'], pars['zeta_K'], 1)
    Hill_Na = Hill(pars['C0_Na_N'], pars['zeta_Na'], 1.5)
    pars['I_NKP_max'] = calc.calculate_I_NKP_max(I_Na_inf_calc, f_NaK_calc, Hill_Na, Hill_K) 
    return pars

#-----------------------------------------------------------------------------------------------------------------------
## SIMULATION preparation 
#-----------------------------------------------------------------------------------------------------------------------
def neuron_simulator():
    indices = array([0, 0, 0, 0])
    times = array([30, 60, 90, 120])*ms
    pre_neuron = SpikeGeneratorGroup(1, indices, times)
    return pre_neuron

def lpc5_neuron(N,params,model='hh-ecs',name='hh*',dt=None, I_max = 0 * nA/cm**2):
    def DESCRIPTION():
        """
    Method that generate the "neuron" model with given parameters.

    Input arguments:
    - N      : int      Number of neurons to simulate
    - params : dict     Dictionary of neuron parameters
    - model  : {'hh-neuron'} | 'hh-ecs'  Model equations
    - name   : String   The name assigned to the NeuronGroup
    - dt     : float*second minimum time step for numerical integration

    Return:
    - neuron : NeuronGroup
    """
    coords = np.array([[0, 0, 0]])
    for i in range(N):
        additional = np.array([[random(), random(), random()]])
        coords = np.append(coords, additional, axis=0)

    # get model equations
    eqs = return_HH_equations(model)

    # setup neuron based on parameters and equations and with refractory time
    neuron = NeuronGroup(N,eqs,
                         threshold='v>v_thr',
                         refractory=3*ms,
                         reset='',
                         namespace=params,
                         name=name,
                         method='euler',
                         dt=dt)

    # Initialize state variables
    neuron.m = calc.gating_variable_m(params['resting_potential']*mvolt)
    neuron.h = calc.gating_variable_h(params['resting_potential']*mvolt)
    neuron.n = calc.gating_variable_n(params['resting_potential']*mvolt)
    neuron.mask = 0
    neuron.mask_noise = 1
    if model == 'hh-ecs' or model == 'hh-ecs-synapse':
        neuron.n_K_N = 0 * mole
        neuron.n_Na_N = 0 * mole
        neuron.n_Cl_N = 0 * mole
    neuron.v = params['resting_potential']*mvolt
    #neuron.v = -90 * mvolt
    for i in range(N):
        neuron[i].x = coords[i, 0]
        neuron[i].y = coords[i, 1]
        neuron[i].z = coords[i, 2]

    return neuron

def exc_synapse(input_signal, neurons_2):     
    # setup synapse parameters
    Synapse_pars = return_synapse_pars()
    
    # setup synapse equations
    #  A                : 1
    synapse_eqs = ''' 
    g_syn_e_post = g_s_exc * r : siemens/meter ** 2 (summed)
    dr/dt = (-r + x_syn * (1-r)* A)/tau_r_x   : 1
    dx_syn/dt = (-x_syn)/tau_x      : 1
    '''
    synapse_action = '''
    x_syn += J*(1-x_syn)
    '''
    # setup synapses
    S = Synapses(input_signal, neurons_2, model=synapse_eqs, on_pre=synapse_action, namespace=Synapse_pars)
    S.connect(p = 0.2)
    # Initialize state variables 
    S.r = 0
    S.x_syn = 0
    return S

def inh_synapse(input_signal, neurons_2):     
    # setup synapse parameters
    Synapse_pars = return_synapse_pars()
    
    # setup synapse equations
    synapse_eqs = ''' 
    g_syn_i_post = g_s_inh * r : siemens/meter ** 2 (summed)
    dr/dt = (-r + x_syn * (1-r) * A)/tau_r_x   : 1
    dx_syn/dt = (-x_syn)/tau_x      : 1 
    '''
    synapse_action = '''
    x_syn += J*(1-x_syn)
    '''
    # setup synapses
    S = Synapses(input_signal, neurons_2, model=synapse_eqs, on_pre=synapse_action, namespace=Synapse_pars)
    S.connect(p = 0.2)
    # Initialize state variables 
    S.r = 0
    S.x_syn = 0
    return S

def lpc5_simulation(duration=1.0,model='hh-ecs', show_monitor = False, save_plots = False, **kwargs):
    def DESCRIPTION_2():
        """
    This is the actual simulation of the Neuron. It is provided as a standalone method for modularity.

    Input parameters:
    - duration   : float*second   Simulated Time
    - model      : {'hh-neuron'} | 'hh-ecs'   Use only 'hh-neuron' for now. Other cases will be added later.
    :param kwargs:
    :return:
    """
    
    # Reinitialization
    device.delete(force=True)
    
    ## BUILD NEURONS AND SYNAPSES
    # Initialize neurons
    n_e = 10
    n_i = 0
    params = lpc5_parameters(model=model,T_ramp=duration*second,**kwargs)
    neurons = lpc5_neuron(n_e + n_i,params,model=model,name='hh*',dt=0.5*us)
    # neurons = lpc5_neuron(n_e + n_i,params,model=model,name='hh*',dt=0.5*us)
    # stimulating_neuron = neurons[0]
    # monitored_neurons = neurons[1:]
    exc_neurons = neurons #[:n_e]
    # inh_neurons = neurons[n_e:]
    neurons.mask = 1

    ## parameter testing:
    decreased_conductances_g_Na = []
    increased_conductances_g_Na = []
    decreased_conductances_g_K = []
    increased_conductances_g_K = []
    for i in range(int(n_e/2)):
        decreased_conductances_g_Na.append(1 ** (int(n_e)/2 - i))
        increased_conductances_g_Na.append(1 ** (i + 1))
        decreased_conductances_g_K.append(0.8 ** (int(n_e)/2 - i))
        increased_conductances_g_K.append(1.2 ** (i + 1))

    test_conductances_g_Na =  decreased_conductances_g_Na + increased_conductances_g_Na
    print(test_conductances_g_Na)
    test_conductances_g_K =  decreased_conductances_g_K + increased_conductances_g_K

    neurons.g_Na_mod = test_conductances_g_Na
    neurons.g_K_mod = test_conductances_g_K
    # for element in test_conductances:
    #     index = test_conductances.index(element)
    #     neurons.g_Na[index] = element

    # Initialize Synapse
    
    
    # exc_synapses = exc_synapse(exc_neurons, neurons)
    # inh_synapses = inh_synapse(inh_neurons, neurons)

    #exc_synapses.delay = 'sqrt((x_pre-x_post)**2 + (y_pre-y_post)**2 + (z_pre-z_post)**2) * 40 * ms' 
    # A_values = [0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8]
    # for element in A_values:
    #     index = A_values.index(element)
    #     exc_synapses.A[index] = element


    # inh_synapses = inh_synapse(inh_neurons, neurons)

    ## SET MONITORS
    #synapse monitor
    if model=='hh-neuron':
        # neuron_variables = ['v', 'm', 'n', 'h']
        neuron_variables = ['v', 'I_noise']
    elif model=='hh-ecs':
        neuron_variables = ['v', 'E_K', 'E_Cl', 'E_Na', 'I_Na', 'I_Na_L', 'I_K', 'I_Cl_L', 'I_KCC', 'I_NKP', 'C_Cl_N','C_Na_N', 'C_K_N', 'C_Cl_E','C_Na_E', 'C_K_E', 'I_inj', 'Check_1', 'Check_2', 'Check_3', 'n_Na_N', 'n_K_N', 'n_Cl_N', 'g_Na_mod']
    elif model=='hh-ecs-synapse':
        neuron_variables = ['v', 'E_K', 'E_Cl', 'E_Na', 'I_Na', 'I_Na_L', 'I_K', 'I_Cl_L', 'I_KCC', 'I_NKP', 'I_syn', 'C_Cl_N','C_Na_N', 'C_K_N', 'C_Cl_E','C_Na_E', 'C_K_E', 'I_inj', 'Check_1', 'Check_2', 'Check_3', 'n_Na_N', 'n_K_N', 'n_Cl_N'] #E_syn
    synapse_variables = ['r', 'x_syn']

    # exc_synapse_monitor = StateMonitor(exc_synapses, variables=synapse_variables,record = [0], dt=0.1*ms, name='exc_synapse_monitor')
    # inh_synapse_monitor = StateMonitor(inh_synapses, variables=synapse_variables,record = [0], dt=0.1*ms, name='inh_synapse_monitor')
    # delay_mon = StateMonitor(exc_synapses.pre, 'delay', record = [0], dt=0.01*ms, name = 'delay_mon')
    #neuron monitor
    sv_mon = StateMonitor(neurons, variables=neuron_variables, record=np.arange(10,dtype=int), dt=0.1*ms, name='svmon_add')
    spike_mon_exc = SpikeMonitor(exc_neurons)
    # spike_mon_inh = SpikeMonitor(inh_neurons)

    ## SIMULATION
    # Gather all objects required for the simulation
    # network = Network([neurons, sv_mon, spike_mon_exc, spike_mon_inh, exc_synapses, exc_synapse_monitor]) #, inh_synapses, inh_synapse_monitor
    network = Network([neurons, sv_mon, spike_mon_exc]) #, spike_mon_inh, exc_synapses, inh_synapses
    
    
    # Run the simulator
    network.run(duration=duration*second,report='text')
    device.build(directory=code_dir, clean=True)
    # print(network.get_states())

    # Set show_monitor -> true if you want to print the monitored variables
    # if show_monitor:
    #     for variable in neuron_variables:
    #         print(variable)
    #         print(getattr(sv_mon, variable)[0][:5])
    #     print(getattr(sv_mon, 't')[:5])

    ## DATA VISUALIZATION

    if model=='hh-neuron':
        plotting_list = return_plotting_list('hh-neuron')
    elif model=='hh-ecs':
        plotting_list = return_plotting_list('hh-ecs')
    elif model=='hh-ecs-synapse':
        plotting_list = return_plotting_list('hh-ecs')

    #fig, ax = plt.subplots(1,1)
    # ax.plot(sv_mon.t,sv_mon.v[0].T)
    # indexes = spike_mon.i==0
    # ax.vlines(spike_mon.t[indexes],spike_mon.i[indexes],spike_mon.i[indexes]+2,linewidth=1,colors='r')
    
    #ax.plot(spike_mon_exc.t,spike_mon_exc.i,"|",lw=0.8, color='k')
    # ax.plot(spike_mon_inh.t,spike_mon_inh.i + 100,"|",lw=0.8, color='r')

    # plot_raster(spike_mon.i, spike_mon.t)
    
    # print(plotting_list)
    fig, ax = plt.subplots(plotting_list[-1]['plot_number']+2, 1,sharex=True, figsize=(8, len(plotting_list) * 3.5))
    # ax[-1].set_xlabel("time (s)")
    for element in plotting_list:
        data = getattr(sv_mon, element['variable'])[:]/element['unit']
        for idx, neuron in enumerate(data):
            # print(neuron)
            ax[element['plot_number']].plot(sv_mon.t/second, neuron, label = 'neuron ' + str(idx))
            #ax[element['plot_number']].plot(sv_mon.t/second, getattr(sv_mon, element['variable'])[:].T/element['unit'], label = element['variable'])    
        ax[element['plot_number']].set_ylabel(element['axis'])
        ax[element['plot_number']].legend()
        # ax[plotting_list[-1]['plot_number']+1].plot(delay_mon.t/second, getattr(delay_mon, 'delay')[:].T/1, label = 'delay')
        #ax[plotting_list[-1]['plot_number']+1].legend()
        # ax[plotting_list[-1]['plot_number']+1].plot(exc_synapse_monitor.t/second, getattr(exc_synapse_monitor, 'x_syn')[:].T/1, label = 'x_syn')
        # ax[plotting_list[-1]['plot_number']+1].legend()
    ax[-1].plot(neurons.g_K_mod, spike_mon_exc.count/(duration*second), label = 'spike frequency')
    plt.legend()
    if save_plots:
        plt.savefig("output/graphs/Test_legend")
    
    
    # visualise_connectivity(exc_synapses)
    # visualise_connectivity(inh_synapses)

    
    device.delete(force=True)
    

if __name__=="__main__":
    #-------------------------------------------------------------------------------------------------------------------
    # Generate Neuron Simulation
    #-------------------------------------------------------------------------------------------------------------------
    device.reinit()
    device.activate()
    set_device('cpp_standalone',directory=code_dir,build_on_run=False)
    prefs.devices.cpp_standalone.openmp_threads = 1 ## The number of threads used in the parallelization (machine-dependent)

    # Multiple HH-neurons connected by synapses
    lpc5_simulation(duration=10,model='hh-ecs', I_dc=10**4.2*nA/cm**2, show_monitor = False, save_plots = True)
    # plt.show()

    # Multiple neurons HH-ecs connected by synapses
    # lpc5_simulation(duration=0.1, model='hh-ecs-synapse', I_dc=10**4.2*nA/cm**2, show_monitor=False, save_plots = True)

    # #-------------------------------------------------------------------------------------------------------------------
    ## Show figures
    # #-------------------------------------------------------------------------------------------------------------------
    # plt.show()