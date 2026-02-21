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
from calculator import Calculator

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
def varargin(pars, pars_mods):
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
    for key, val in pars_mods.items():
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

def show_monitor_vars(monitor, variables):
     for variable in variables:
            print(variable)
            readout = getattr(monitor, variable)[0][:5]
            print(readout)
            print(readout.unit)

def calculate_network_connectivities(exc_to_exc_synapses, exc_to_inh_synapses, inh_to_exc_synapses, inh_to_inh_synapses, n_e, n_i):
    total_connectivity = np.zeros(n_e + n_i)
    for target in exc_to_exc_synapses.j:
        total_connectivity[target] += 1
    for target in exc_to_inh_synapses.j:
        total_connectivity[target + n_e] += 1
    for target in inh_to_exc_synapses.j:
        total_connectivity[target] -=1
    for target in inh_to_inh_synapses.j:
        total_connectivity[target + n_e] -= 1
    inhibitory_connectivity = np.zeros(n_e + n_i)
    for target in inh_to_exc_synapses.j:
        inhibitory_connectivity[target] -=1
    for target in inh_to_inh_synapses.j:
        inhibitory_connectivity[target + n_e] -= 1
    return total_connectivity, inhibitory_connectivity
    
def params_handover(model):
    return return_initial_parameters(model)   
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

# V_T = ThermalVoltage(36)*1000
# print(V_T)

V_T = -56.2

def alpha_m(V):
    return -0.32 * (V - V_T - 13)/(np.expm1(-(V - V_T - 13)/4))
def beta_m(V):
    return 0.28 * (V - V_T - 40)/(np.expm1((V - V_T - 40)/5))
def m_inf(V):
    return alpha_m(V)/(alpha_m(V) + beta_m(V))

def alpha_h(V):
    return 0.128 * np.exp((-V - V_T - 17)/18)
def beta_h(V):
    return 4/(1 + np.exp(-(V - V_T - 40)/5))
def h_inf(V):
    return alpha_h(V)/(alpha_h(V) + beta_h(V))

def alpha_n(V):
    return -0.032 * (V - V_T - 15)/(np.exp((-(V - V_T - 15)/5))- 1)
def beta_n(V):
    return 0.5 * np.exp(-(V -V_T - 10)/40)
def n_inf(V):
    return alpha_n(V)/(alpha_n(V) + beta_n(V))

def p_inf(V):
    return 1/(1 + np.exp(- (V + 35)/10))

def set_NeuronGroup_state_variables(NeuronGroup, N, model, params):
    NeuronGroup.v = params['resting_potential'] * mvolt
    NeuronGroup.m = m_inf(params['resting_potential'])
    NeuronGroup.h = h_inf(params['resting_potential'])
    NeuronGroup.n = n_inf(params['resting_potential'])
    NeuronGroup.p = p_inf(params['resting_potential'])
    NeuronGroup.I_max = params['I_max']
    NeuronGroup.mask_noise = 1

    coords = np.array([[0, 0, 0]])
    for i in range(N):
        additional = np.array([[random(), random(), random()]])
        coords = np.append(coords, additional, axis=0)
    for i in range(N):
        NeuronGroup[i].x = coords[i, 0]
        NeuronGroup[i].y = coords[i, 1]
        NeuronGroup[i].z = coords[i, 2]
    
    if model == 'hh-ecs' or model == 'hh-ecs-synapse':
        NeuronGroup.n_K_N = 0 * mole
        NeuronGroup.n_Na_N = 0 * mole
        NeuronGroup.n_Cl_N = 0 * mole

#-----------------------------------------------------------------------------------------------------------------------
## PARAMETER preparation 
#-----------------------------------------------------------------------------------------------------------------------
def lpc5_parameters(model, **kwargs):
    ## Get parameters from parameters.py
    pars = return_initial_parameters(model)
    if model == 'hh-ecs_exc' or model == 'hh-ecs_inh' or model == 'hh-ecs-synapse_exc' or model == 'hh-ecs-synapse_inh' :
        calc = Calculator(model)
        print(calc.params)
        ## Conductance sodium leakage channel
        pars['g_Na_L'] = calc.calc_leakage_conductance()
        ## Adding I_NKP_max
        I_Na_inf_calc = calc.I_Na_inf(pars['g_Na_L'])
        # I_Na_inf_calc = calc.I_Na_inf(calc.calc_leakage_conductance(), pars['resting_potential']*mvolt)
        f_NaK_calc = calc.f_NaK()
        # f_NaK_calc = calc.f_NaK(pars['resting_potential']*mvolt, ThermalVoltage(pars['T_exp'])*1000*mvolt)
        Hill_K = Hill(pars['C0_K_E'], pars['zeta_K'], 1)
        Hill_Na = Hill(pars['C0_Na_N'], pars['zeta_Na'], 1.5)
        pars['I_NKP_max'] = calc.calculate_I_NKP_max(I_Na_inf_calc, f_NaK_calc, Hill_Na, Hill_K)
    print("Parameter generation successful------------------------------------------------------------------------------")
    print(pars['g_Na_L'])
    print(pars['resting_potential'])
    return pars

#-----------------------------------------------------------------------------------------------------------------------
## Neuron Setup
#-----------------------------------------------------------------------------------------------------------------------

def lpc5_neuron(N, params, model='hh-ecs', name='hh*', dt=None):
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
    # Get model equations
    eqs = return_HH_equations(model)
    print("Equation generation successful-------------------------------------------------------------------------------")

    # Initialize NeuronGroup
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
    set_NeuronGroup_state_variables(neuron, N, model, params)
    print("State variables initialized----------------------------------------------------------------------------------")
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
    S.r = 0.1
    S.x_syn = 0
    return S

#-----------------------------------------------------------------------------------------------------------------------
## Simulation 
#-----------------------------------------------------------------------------------------------------------------------
def lpc5_simulation(model_dict, show_monitor = False, save_plots = False, data_name = None, exc_params_mods = None, inh_params_mods = None, **kwargs):
    def DESCRIPTION_2():
        """
    This is the actual simulation of the Neuron. It is provided as a standalone method for modularity.

    Input parameters:
    - duration   : float*second   Simulated Time
    - model      : {'hh-neuron'} | 'hh-ecs'   Use only 'hh-neuron' for now. Other cases will be added later.
    :param kwargs:
    :return:
    """

    model = model_dict['model']
    n_e = model_dict['n_e']
    n_i = model_dict['n_i']
    duration = model_dict['duration']
    ## BUILD NEURONS AND SYNAPSES
    # INITIALIZE NEURONS------------------------------------------------------------------------------------------------------------------------------------------------
    if exc_params_mods is not None:
        exc_params = varargin(lpc5_parameters(model=model + '_exc', T_ramp=duration*second, **kwargs), exc_params_mods)
    elif exc_params_mods is None:
        #print('KLAPPT')
        exc_params = lpc5_parameters(model=model + '_exc', T_ramp=duration*second, **kwargs)
        # print(exc_params)
    if inh_params_mods is not None:
        inh_params = varargin(lpc5_parameters(model=model + '_inh', T_ramp=duration*second, **kwargs), inh_params_mods)
    elif inh_params_mods is None:
        inh_params = lpc5_parameters(model=model + '_inh', T_ramp=duration*second, **kwargs)
    
    if model == 'hh-neuron' or model == 'hh-ecs':
        if n_e > 0 and n_i > 0:
            exc_neurons = lpc5_neuron(n_e, exc_params, model=model, name='NGexc', dt=0.1*us)
            inh_neurons = lpc5_neuron(n_i, inh_params, model=model, name='NGinh', dt=0.1*us)
        elif not n_i > 0:
            exc_neurons = lpc5_neuron(n_e, exc_params, model=model, name='NGexc', dt=0.1*us)
        elif not n_e > 0:
            inh_neurons = lpc5_neuron(n_i, inh_params, model=model, name='NGinh', dt=0.1*us)
    elif model == 'hh-neuron-synapse' or model == 'hh-ecs-synapse':
        exc_neurons = lpc5_neuron(n_e , exc_params, model=model, name='NGexc', dt=0.1*us)
        inh_neurons = lpc5_neuron(n_i, inh_params, model=model, name='NGinh', dt=0.1*us)

    # PARAMETER PREPARATION
    # I_max_mods = [params['I_max']/(namp/cm**2)] * n
    # I_max_mods = np.multiply(modifier, I_max_mods) * (namp/cm**2)
    print("Parameter initialization successful---------------------------------------------------------------------------")

    ## INITIALIZE SYNAPSES------------------------------------------------------------------------------------------------------------------------------------------------
    if model == 'hh-neuron-synapse' or model == 'hh-ecs-synapse':
        exc_to_exc_synapses = exc_synapse(exc_neurons, exc_neurons)
        exc_to_inh_synapses = exc_synapse(exc_neurons, inh_neurons)
        inh_to_exc_synapses = inh_synapse(inh_neurons, exc_neurons)
        inh_to_inh_synapses = inh_synapse(inh_neurons, inh_neurons)
        # exc_to_exc_synapses.delay = 'sqrt((x_pre-x_post)**2 + (y_pre-y_post)**2 + (z_pre-z_post)**2) * 40 * ms' 
        # exc_to_inh_synapses.delay = 'sqrt((x_pre-x_post)**2 + (y_pre-y_post)**2 + (z_pre-z_post)**2) * 40 * ms' 
        # inh_to_exc_synapses.delay = 'sqrt((x_pre-x_post)**2 + (y_pre-y_post)**2 + (z_pre-z_post)**2) * 40 * ms' 
        # inh_to_inh_synapses.delay = 'sqrt((x_pre-x_post)**2 + (y_pre-y_post)**2 + (z_pre-z_post)**2) * 40 * ms' 
        print('Synapses and delay successfully initialized-------------------------------------------------------------------')

    ## SET MONITORS-------------------------------------------------------------------------------------------------------------------------------------------------------
    # Define monitored variables
    if model=='hh-neuron':
        neuron_variables = ['v', 'I_inj', 'I_Na', 'I_Kd', 'I_M', 'I_L', 'm', 'n', 'h', 'p', 'a_m', 'b_m', 'p_inf', 'tau_m', 'tau_h', 'tau_n', 'tau_p', 'I_max'] # I_noise
    elif model=='hh-ecs':
        neuron_variables = ['v', 'E_K', 'E_Cl', 'E_Na', 'I_Na', 'I_Na_L', 'I_Kd', 'I_M', 'I_Cl_L', 'I_KCC', 'I_NKP', 'C_Cl_N','C_Na_N', 'C_K_N', 'C_Cl_E','C_Na_E', 'C_K_E', 'I_inj', 'Check_1', 'Check_2', 'Check_3', 'n_Na_N', 'n_K_N', 'n_Cl_N', 'm', 'n', 'h', 'p', 'a_m', 'b_m', 'tau_m', 'tau_h', 'tau_n', 'tau_p', 'I_max']
    elif model == 'hh-neuron-synapse':
        neuron_variables = ['v', 'I_inj', 'I_Na', 'I_Kd', 'I_M', 'I_L', 'm', 'n', 'h', 'p', 'a_m', 'b_m', 'p_inf', 'tau_m', 'tau_h', 'tau_n', 'tau_p', 'I_max'] # I_noise
        synapse_variables = ['r', 'x_syn']
    elif model=='hh-ecs-synapse':
        neuron_variables = ['v', 'E_K', 'E_Cl', 'E_Na', 'I_Na', 'I_Na_L', 'I_Kd', 'I_M', 'I_Cl_L', 'I_KCC', 'I_NKP', 'C_Cl_N','C_Na_N', 'C_K_N', 'C_Cl_E','C_Na_E', 'C_K_E', 'I_inj', 'Check_1', 'Check_2', 'Check_3', 'n_Na_N', 'n_K_N', 'n_Cl_N', 'm', 'n', 'h', 'p', 'a_m', 'b_m', 'tau_m', 'tau_h', 'tau_n', 'tau_p', 'I_max']
        synapse_variables = ['r', 'x_syn']
        
    # Define monitors
    if model == 'hh-neuron' or model =='hh-ecs':
        if n_e > 0 and n_i > 0:
            sv_mon_exc = StateMonitor(exc_neurons, variables=neuron_variables, record=np.arange(n_e,dtype=int), dt=0.01*ms, name='svmon_exc_add')
            sv_mon_inh = StateMonitor(inh_neurons, variables=neuron_variables, record=np.arange(n_i,dtype=int), dt=0.01*ms, name='svmon_inh_add')
            spike_mon_exc = SpikeMonitor(exc_neurons)
            spike_mon_inh = SpikeMonitor(inh_neurons)
        elif not n_i > 0:
            sv_mon_exc = StateMonitor(exc_neurons, variables=neuron_variables, record=np.arange(n_e,dtype=int), dt=0.01*ms, name='svmon_exc_add')
            spike_mon_exc = SpikeMonitor(exc_neurons)
        elif not n_e > 0:
            sv_mon_inh = StateMonitor(inh_neurons, variables=neuron_variables, record=np.arange(n_i,dtype=int), dt=0.01*ms, name='svmon_inh_add')
            spike_mon_inh = SpikeMonitor(inh_neurons)
    elif model == 'hh-neuron-synapse' or model == 'hh-ecs-synapse':
        sv_mon_exc = StateMonitor(exc_neurons, variables=neuron_variables, record=np.arange(n_e,dtype=int), dt=0.01*ms, name='svmon_exc_add')
        sv_mon_inh = StateMonitor(inh_neurons, variables=neuron_variables, record=np.arange(n_i,dtype=int), dt=0.01*ms, name='svmon_inh_add')
        rate_mon_exc = PopulationRateMonitor(exc_neurons)
        rate_mon_inh = PopulationRateMonitor(inh_neurons)
        spike_mon_exc = SpikeMonitor(exc_neurons)
        spike_mon_inh = SpikeMonitor(inh_neurons)
        # exc_synapse_monitor = StateMonitor(exc_synapses, synapse_variables, record= True)
        # inh_synapse_monitor = StateMonitor(inh_synapses, synapse_variables, record= True)
    print('Monitor setup successful---------------------------------------------------------------------------')
        
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
        print('check')
        network = Network([exc_neurons, inh_neurons, sv_mon_exc, sv_mon_inh, spike_mon_exc, spike_mon_inh, rate_mon_exc, rate_mon_inh, exc_to_exc_synapses, exc_to_inh_synapses, inh_to_inh_synapses, inh_to_exc_synapses]) #exc_synapse_monitor, inh_synapse_monitor
    
    # Run the simulator
    network.run(duration=duration*second,report='text')
    device.build(directory=code_dir, clean=True)
    # print(network.get_states())
    print('Network run complete-------------------------------------------------------------------------------')


    # Set show_monitor -> True if you want to print the monitored variables
    if show_monitor:
        show_monitor_vars(sv_mon_exc, neuron_variables)


    # bins, rates = rate_mon_exc.
    # bin_centers = bins + 10*ms/2
    fig, ax = plt.subplots(1, 2, figsize = (10, 5))
    ax[0].plot(sv_mon_exc.t/second, sv_mon_exc.v[:].T/mvolt)
    ax[1].plot(sv_mon_inh.t/second, sv_mon_inh.v[:].T/mvolt)
    #visualise_connectivity(exc_to_exc_synapses)
    # ax[0].plot(bin_centers, rates)
    # ax[1].plot(rate_mon_exc.t/second, rate_mon_exc.smooth_rate(window='flat', width=0.5*ms)/Hz)
    #ax[1].plot(sv_mon_inh.t/second, sv_mon_inh.v[:].T/mvolt)    plt.show 


    if model == 'hh-neuron-synapse':
        total_connectivity, inhibitory_connectivity = calculate_network_connectivities(exc_to_exc_synapses, exc_to_inh_synapses, inh_to_exc_synapses, inh_to_inh_synapses, n_e, n_i)

    # if model == 'hh-neuron-synapse':
    #     store_data(model, data_name, [sv_mon_exc, sv_mon_inh, spike_mon_exc, spike_mon_inh], n_e, n_i, total_connectivity, inhibitory_connectivity)
    #     print('data stored successfully---------------------------------------------------------------------------')
    # elif model == 'hh-ecs':
    #     store_data(model, data_name, [sv_mon_exc, sv_mon_inh, spike_mon_exc, spike_mon_inh], n_e, n_i, total_connectivity = None, inhibitory_connectivity = None)
    # elif model == 'hh-ecs-synapse':
    #     store_data('hh-ecs', data_name, [sv_mon_exc, sv_mon_inh, spike_mon_exc, spike_mon_inh], n_e, n_i)
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
    
    ## REQUIRED INFORMATION TO DEFINE MODEL: Model type (-> String), Number of exc neurons (-> Integer), Number of inh neurons (-> Integer), exc_params_adj (-> Dicitionary), duration 

    model_dict = {
        'model' : 'hh-ecs',
        'n_e' : 1,
        'n_i' : 1,
        'duration' : 1
    }

    lpc5_simulation(model_dict, 
                    data_name = 'data/hh-neuron-synapse_test'
                    )
    
    # device.reinit()
    # device.activate()
    # set_device('cpp_standalone',directory=code_dir,build_on_run=False)
    # prefs.devices.cpp_standalone.openmp_threads = 1 ## The number of threads used in the parallelization (machine-dependent)


    # model_dict['duration'] = 0.5
    # exc_params_mods = {
    #     'I_max' : 50000 * namp/cm**2 
    # }
    # lpc5_simulation(model_dict, 
    #                 data_name = 'data/malgucken', 
    #                 exc_params_mods= exc_params_mods
    #                 )
    plt.show()