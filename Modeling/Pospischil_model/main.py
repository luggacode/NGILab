#-----------------------------------------------------------------------------------------------------------------------
# Default Imports
#-----------------------------------------------------------------------------------------------------------------------
import numpy as np
import json
from pathlib import Path
import importlib
from brian2tools import *
from importnb import Notebook
from numpy import random
import matplotlib.pyplot as plt
import plotly.io as pio
import plotly.graph_objects as go
from plotly.subplots import make_subplots
pio.renderers.default = "browser"

#-----------------------------------------------------------------------------------------------------------------------
# File imports
#-----------------------------------------------------------------------------------------------------------------------
from equations import return_HH_equations
from parameters import return_initial_parameters, return_synapse_pars
from plotting import return_plotting_list
from data_storing import store_data
import calculator as calc

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
def lpc5_parameters(model, **kwargs):
    ## Get parameters from parameters.py
    pars = return_initial_parameters(model)
    if model == 'hh-ecs_exc' or model == 'hh-ecs_exc':
        ## Conductance sodium leakage channel
        pars['g_Na_L'] = calc.calc_leakage_conductance()

        ## Adding I_NKP_max
        I_Na_inf_calc = calc.I_Na_inf(calc.calc_leakage_conductance(), pars['resting_potential']*mvolt)
        f_NaK_calc = calc.f_NaK(pars['resting_potential']*mvolt, ThermalVoltage(pars['T_exp'])*1000*mvolt)
        Hill_K = Hill(pars['C0_K_E'], pars['zeta_K'], 1)
        Hill_Na = Hill(pars['C0_Na_N'], pars['zeta_Na'], 1.5)
        pars['I_NKP_max'] = calc.calculate_I_NKP_max(I_Na_inf_calc, f_NaK_calc, Hill_Na, Hill_K) 
    print("Parameter generation successful------------------------------------------------------------------------------")
    return pars

#-----------------------------------------------------------------------------------------------------------------------
## Neuron Setup
#-----------------------------------------------------------------------------------------------------------------------

def lpc5_neuron(N,params,model='hh-ecs',name='hh*',dt=None, I_max = 0 * nA/cm**2):
    def DESCRIPTION():
        """
    Method that generate the "neuron" model with given parameters.

    Input arguments:
    - N      : int      Number of neurons to simulate
    - params : dict     Dictionary of neuron parameters
    - model  : 'hh-neuron' | 'hh-ecs' | 'hh-ecs-synapse' Model equations
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
    print("Equation generation successful-------------------------------------------------------------------------------")

    # setup neuron based on parameters and equations and with refractory time
    # threshold='v>v_thr',
    neuron = NeuronGroup(N,eqs, 
                         threshold = 'v>v_thr',
                         refractory='v >= v_thr',
                         reset='',
                         namespace=params,
                         name=name,
                         method='euler',
                         dt=dt)
    print("NeuronGroup generation successful----------------------------------------------------------------------------")

    # Initialize state variables
    neuron.m = calc.gating_variable_m(params['resting_potential'] * mvolt, params['V_T'])
    neuron.h = calc.gating_variable_h(params['resting_potential'] * mvolt, params['V_T'])
    neuron.n = calc.gating_variable_n(params['resting_potential'] * mvolt, params['V_T'])
    neuron.p = calc.gating_variable_p(params['resting_potential'] * mvolt)
    neuron.I_max = params['I_max']
    # neuron.mask = 1
    neuron.mask_noise = 1
    # neuron.I_inj_mod = 1
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
    print("Gating variables initialized---------------------------------------------------------------------------------")
    return neuron

#-----------------------------------------------------------------------------------------------------------------------
## Synapses Setup
#-----------------------------------------------------------------------------------------------------------------------
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
    S = Synapses(input_signal, neurons_2, model=synapse_eqs, on_pre=synapse_action, namespace=Synapse_pars, method='euler')
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
    S = Synapses(input_signal, neurons_2, model=synapse_eqs, on_pre=synapse_action, namespace=Synapse_pars, method='euler')
    S.connect(p = 0.2)
    # Initialize state variables 
    S.r = 0
    S.x_syn = 0
    return S

#-----------------------------------------------------------------------------------------------------------------------
## Simulation 
#-----------------------------------------------------------------------------------------------------------------------
def lpc5_simulation(duration=10.0,model='hh-ecs', show_monitor = False, save_plots = False, data_name = '', **kwargs):
    def DESCRIPTION_2():
        """
    This is the actual simulation of the Neuron. It is provided as a standalone method for modularity.

    Input parameters:
    - duration   : float*second   Simulated Time
    - model      : {'hh-neuron'} | 'hh-ecs'   Use only 'hh-neuron' for now. Other cases will be added later.
    :param kwargs:
    :return:
    """
    
    # Reinitialization to allow multiple runs
    device.delete(force=True)
    print(model)
    ## BUILD NEURONS AND SYNAPSES
    # INITIALIZE NEURONS------------------------------------------------------------------------------------------------------------------------------------------------
    exc_params = lpc5_parameters(model=model + '_exc', T_ramp=duration*second, **kwargs)
    inh_params = lpc5_parameters(model=model + '_inh', T_ramp=duration*second, **kwargs)
    if model == 'hh-neuron' or model == 'hh-ecs':
        n_e = 1
        n_i = 1
        if n_e > 0 and n_i > 0:
            exc_neurons = lpc5_neuron(n_e, exc_params, model=model, name='hh*', dt=0.1*us)
            inh_neurons = lpc5_neuron(n_i, inh_params, model=model, name='hh*', dt=0.1*us)
        elif not n_i > 0:
            exc_neurons = lpc5_neuron(n_e, exc_params, model=model, name='hh*', dt=0.1*us)
        elif not n_e > 0:
            inh_neurons = lpc5_neuron(n_i, inh_params, model=model, name='hh*', dt=0.1*us)
    elif model == 'hh-neuron-synapse' or model == 'hh-ecs-synapse':
        n_e = 16
        n_i = 4
        neurons = lpc5_neuron(n_e + n_i, exc_params, model=model, name='hh*', dt=0.1*us)
        exc_neurons = neurons[:n_e]
        inh_neurons = neurons[n_e:]

    

     # PARAMETER PREPARATION
    if model == 'hh-neuron' or model == 'hh-ecs':
        if n_e > 0 and n_i > 0:
            equilibrium_potential = -70 * mvolt
            exc_neurons.v = equilibrium_potential
            exc_neurons.m = calc.gating_variable_m(equilibrium_potential, exc_params['V_T'])
            exc_neurons.h = calc.gating_variable_h(equilibrium_potential, exc_params['V_T'])
            exc_neurons.n = calc.gating_variable_n(equilibrium_potential, exc_params['V_T'])
            exc_neurons.p = calc.gating_variable_p(exc_params['V_T'])
            exc_neurons.g_Na_mod = 1 
            exc_neurons.g_Kd_mod = 1
            exc_neurons.g_M_mod = 1
            inh_neurons.v = equilibrium_potential
            inh_neurons.m = calc.gating_variable_m(equilibrium_potential, inh_params['V_T'])
            inh_neurons.h = calc.gating_variable_h(equilibrium_potential, inh_params['V_T'])
            inh_neurons.n = calc.gating_variable_n(equilibrium_potential, inh_params['V_T'])
            inh_neurons.p = calc.gating_variable_p(inh_params['V_T'])
            inh_neurons.g_Na_mod = 1 
            inh_neurons.g_Na_mod = 1 
            inh_neurons.g_Kd_mod = 1
            inh_neurons.g_M_mod = 1
        elif not n_i > 0:
            exc_neurons.g_Na_mod = 1 
            exc_neurons.g_Kd_mod = 1
            exc_neurons.g_M_mod = 1
        elif not n_e > 0:
            inh_neurons.g_Na_mod = 1 
            inh_neurons.g_Kd_mod = 1
            inh_neurons.g_M_mod = 1
    elif model == 'hh-neuron-synapse' or model == 'hh-ecs-synapse':
        neurons.g_Na_mod = 1
        neurons.g_Kd_mod = 1
        neurons.g_M_mod = 1
   
    # modifier = [1, 4, 5, 6, 7, 8.5, 10, 12.5, 15, 18, 25, 32.5, 40, 50,  60, 75, 90, 105, 112.5, 120] * 3
    # I_max_mods = [params['I_max']/(namp/cm**2)] * n
    # I_max_mods = np.multiply(modifier, I_max_mods) * (namp/cm**2)
    # neurons.I_max = I_max_mods
    print("Parameter initialization successful---------------------------------------------------------------------------")

    ## INITIALIZE SYNAPSES------------------------------------------------------------------------------------------------------------------------------------------------
    if model == 'hh-neuron-synapse' or model == 'hh-ecs-synapse':
        exc_synapses = exc_synapse(exc_neurons, neurons)
        inh_synapses = inh_synapse(inh_neurons, neurons)
        exc_synapses.delay = 'sqrt((x_pre-x_post)**2 + (y_pre-y_post)**2 + (z_pre-z_post)**2) * 40 * ms' 
        # A_values = [0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8]
        # for element in A_values:
        #     index = A_values.index(element)
        #     exc_synapses.A[index] = element

    ## SET MONITORS-------------------------------------------------------------------------------------------------------------------------------------------------------
    # Define monitored variables
    if model=='hh-neuron':
        neuron_variables = ['v', 'I_inj', 'I_Na', 'I_Kd', 'I_M', 'I_L', 'm', 'n', 'h', 'p', 'a_m', 'b_m', 'p_inf', 'tau_m', 'tau_h', 'tau_n', 'tau_p', 'I_max'] # I_noise
    elif model=='hh-ecs':
        neuron_variables = ['v', 'E_K', 'E_Cl', 'E_Na', 'I_Na', 'I_Na_L', 'I_Kd', 'I_M', 'I_Cl_L', 'I_KCC', 'I_NKP', 'C_Cl_N','C_Na_N', 'C_K_N', 'C_Cl_E','C_Na_E', 'C_K_E', 'I_inj', 'Check_1', 'Check_2', 'Check_3', 'n_Na_N', 'n_K_N', 'n_Cl_N', 'g_Na_mod', 'm', 'n', 'h', 'p', 'a_m', 'b_m']
    elif model == 'hh-neuron-synapse':
        neuron_variables = ['v', 'I_inj', 'I_Na', 'I_Kd', 'I_M', 'I_L', 'm', 'n', 'h', 'p', 'a_m', 'b_m', 'p_inf', 'tau_m', 'tau_h', 'tau_n', 'tau_p', 'I_max'] # I_noise
        synapse_variables = ['r', 'x_syn']
    elif model=='hh-ecs-synapse':
        neuron_variables = ['v', 'E_K', 'E_Cl', 'E_Na', 'I_Na', 'I_Na_L', 'I_K', 'I_Cl_L', 'I_KCC', 'I_NKP', 'I_syn', 'C_Cl_N','C_Na_N', 'C_K_N', 'C_Cl_E','C_Na_E', 'C_K_E', 'I_inj', 'Check_1', 'Check_2', 'Check_3', 'n_Na_N', 'n_K_N', 'n_Cl_N'] #E_syn
        synapse_variables = ['r', 'x_syn']
        
    # Define monitors
    if model == 'hh-neuron' or model =='hh-ecs':
        if n_e > 0 and n_i > 0:
            sv_mon_exc = StateMonitor(exc_neurons, variables=neuron_variables, record=np.arange(n_e,dtype=int), dt=0.01*ms, name='svmon_exc_add')
            sv_mon_inh = StateMonitor(inh_neurons, variables=neuron_variables, record=np.arange(n_i,dtype=int), dt=0.01*ms, name='svmon_inh_add')
        elif not n_i > 0:
            sv_mon_exc = StateMonitor(exc_neurons, variables=neuron_variables, record=np.arange(n_e,dtype=int), dt=0.01*ms, name='svmon_exc_add')
        elif not n_e > 0:
            sv_mon_inh = StateMonitor(inh_neurons, variables=neuron_variables, record=np.arange(n_i,dtype=int), dt=0.01*ms, name='svmon_inh_add')
    elif model == 'hh-neuron-synapse' or model == 'hh-ecs-synapse':
        sv_mon = StateMonitor(neurons, variables=neuron_variables, record=np.arange(n_e + n_i,dtype=int), dt=0.01*ms, name='svmon_add')
    if model == 'hh-neuron' or model == 'hh-ecs':
        if n_e > 0 and n_i > 0:
            spike_mon_exc = SpikeMonitor(exc_neurons)
            spike_mon_inh = SpikeMonitor(inh_neurons)
        elif not n_i > 0:
            spike_mon_exc = SpikeMonitor(exc_neurons)
        elif not n_e > 0:
            spike_mon_inh = SpikeMonitor(inh_neurons)
    elif model == 'hh-neuron-synapse' or model == 'hh-ecs-synapse':
        # exc_synapse_monitor = StateMonitor(exc_synapses, synapse_variables, record= True)
        # inh_synapse_monitor = StateMonitor(inh_synapses, synapse_variables, record= True)
        spike_mon_exc = SpikeMonitor(exc_neurons)
        spike_mon_inh = SpikeMonitor(inh_neurons)


    ## SIMULATION---------------------------------------------------------------------------------------------------------------------------------------------------------
    # Gather all objects required for the simulation
    if model == 'hh-neuron' or model == 'hh-ecs':
        if n_e > 0 and n_i > 0:
            network = Network([exc_neurons, inh_neurons, sv_mon_exc, sv_mon_inh, spike_mon_exc, spike_mon_inh]) 
        elif not n_i > 0:
            network = Network([exc_neurons, sv_mon_exc, spike_mon_exc]) 
        elif not n_e > 0:
            network = Network([inh_neurons, sv_mon_inh, spike_mon_inh]) #, spike_mon_inh, exc_synapses, inh_synapses
    elif model == 'hh-neuron-synapse' or model == 'hh-ecs-synapse':
        network = Network([neurons, sv_mon, spike_mon_exc, spike_mon_inh, exc_synapses, inh_synapses]) #exc_synapse_monitor, inh_synapse_monitor
    
    # Run the simulator
    network.run(duration=duration*second,report='text')
    device.build(directory=code_dir, clean=True)
    # print(network.get_states())

    # Set show_monitor -> true if you want to print the monitored variables
    if show_monitor:
        for variable in neuron_variables:
            print(variable)
            readout = getattr(sv_mon_exc, variable)[0][:5]
            print(readout)
            print(readout.unit)
        print(getattr(sv_mon_exc, 't')[:5])

    ## DATA VISUALIZATION---------------------------------------------------------------------------------------------------------------------------------------------------
    
    # Get plottting list for required variables 
    if model=='hh-neuron' or model == 'hh-neuron-synapse':
        plotting_list = return_plotting_list('hh-neuron')
    elif model=='hh-ecs':
        plotting_list = return_plotting_list('hh-ecs')
    elif model=='hh-ecs-synapse':
        plotting_list = return_plotting_list('hh-ecs')

    # indexes = spike_mon.i==0
    # ax.vlines(spike_mon.t[indexes],spike_mon.i[indexes],spike_mon.i[indexes]+2,linewidth=1,colors='r')
    
    fig, ax = plt.subplots(1,2, figsize = (10, 5))
    ax[0].plot(sv_mon_exc.t/second, sv_mon_exc.v[:].T/mvolt)
    ax[1].plot(sv_mon_inh.t/second, sv_mon_inh.v[:].T/mvolt)
    plt.show

    if model == 'hh-neuron-synapse' or model == 'hh-ecs-synapse':
        # fig, ax = plt.subplots(1,1)
        # ax.plot(spike_mon_exc.t,spike_mon_exc.i,"|",lw=0.1, color='k', label = 'neuron index')
        # ax.plot(spike_mon_inh.t,spike_mon_inh.i + n_e,"|",lw=0.1, color='r', label = 'time (s)')
        # ax.set_xlabel('time (s)', fontsize = 11)
        # ax.set_ylabel('neuron idx', fontsize = 11)
        total_connectivity = np.zeros(n_e + n_i)
        print(exc_synapses.j)
        for target in exc_synapses.j:
            total_connectivity[target] += 1
        for target in inh_synapses.j:
            total_connectivity[target] -=1
        print(total_connectivity)
        # print(exc_synapses.i)
        # print(exc_synapses.j)
        # print(inh_synapses.i)
        # visualise_connectivity(exc_synapses)
        # visualise_connectivity(inh_synapses)
        # plt.show()  

    store_data('malgucken', [sv_mon_exc, sv_mon_inh, spike_mon_exc, spike_mon_inh], n_e, n_i)
    
    # if store_data:
    #     v_mV = np.asarray(sv_mon.v / mvolt)     # shape: (n_neurons, n_timepoints)
    #     t_ms = np.asarray(sv_mon.t / second)        # shape: (n_timepoints,)
    #     I_inj = np.asarray(sv_mon.I_inj / (namp/cm**2))
    #     # I_max = np.asarray(I_max_mods/(namp/cm**2))
    #     I_Na = np.asarray(sv_mon.I_Na / (namp/cm**2))
    #     I_Kd = np.asarray(sv_mon.I_Kd / (namp/cm**2))
    #     I_M = np.asarray(sv_mon.I_M / (namp/cm**2))
    #     I_L = np.asarray(sv_mon.I_L / (namp/cm**2))
    #     m = np.asarray(sv_mon.m)
    #     h = np.asarray(sv_mon.h)
    #     n = np.asarray(sv_mon.n)
    #     p = np.asarray(sv_mon.p)
    #     tau_m = np.asarray(sv_mon.tau_m / second)
    #     tau_h = np.asarray(sv_mon.tau_h / second)
    #     tau_n = np.asarray(sv_mon.tau_n / second)
    #     tau_p = np.asarray(sv_mon.tau_p / second)
    #     # spike_times = np.asarray(spike_mon.spike_trains()[0])
    #     # spike_frequency = np.asarray(spike_mon.count/duration)
    #     if model == 'hh-neuron-synapse' or model == 'hh-ecs-synapse':
    #         spike_mon_exc_t = np.asarray(spike_mon_exc.t / second)
    #         spike_mon_exc_i = np.asarray(spike_mon_exc.i)
    #         spike_mon_inh_t = np.asarray(spike_mon_inh.t / second)
    #         spike_mon_inh_i = np.asarray(spike_mon_inh.i)
    #         number_neurons = np.asarray(n_e + n_i)
    #         number_exc_neurons = np.asarray(n_e)
    #         number_inh_neurons = np.asarray(n_i)
    #     np.savez(data_name + ".npz",
    #      v=v_mV,
    #      t=t_ms,
    #     #  I_max = I_max,
    #      I_inj=I_inj,
    #      I_Na=I_Na,
    #      I_Kd=I_Kd,
    #      I_M=I_M,
    #      I_L=I_L,
    #      m = m,
    #      h = h,
    #      n = n,
    #      p = p,
    #      tau_m = tau_m,
    #      tau_h = tau_h,
    #      tau_n = tau_n,
    #      tau_p = tau_p,
    #     #  spike_times = spike_times,
    #     #  spike_frequency = spike_frequency,
    #      ## below just for synaptic connections
    #      spike_mon_exc_t = spike_mon_exc_t,
    #      spike_mon_exc_i = spike_mon_exc_i,
    #      spike_mon_inh_t = spike_mon_inh_t,
    #      spike_mon_inh_i = spike_mon_inh_i,
    #      number_neurons = number_neurons,
    #      number_exc_neurons = number_exc_neurons,
    #      number_inh_neurons = number_inh_neurons,
    #      total_connectivity = total_connectivity
    #      )

    if save_plots:
        plt.savefig("output/graphs/Test_legend")
        
    device.delete(force=True)
    
#-------------------------------------------------------------------------------------------------------------------
# Generate Neuron Simulation
#-------------------------------------------------------------------------------------------------------------------
if __name__=="__main__":
    device.reinit()
    device.activate()
    set_device('cpp_standalone',directory=code_dir,build_on_run=False)
    prefs.devices.cpp_standalone.openmp_threads = 1 ## The number of threads used in the parallelization (machine-dependent)
    # Neuron model after Pospischil implementation and testing
    lpc5_simulation(duration=1,
                    model='hh-neuron', 
                    show_monitor = False, 
                    save_plots = False, 
                    data_name = 'kek')
    plt.show()