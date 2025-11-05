#-----------------------------------------------------------------------------------------------------------------------
# Default Imports
#-----------------------------------------------------------------------------------------------------------------------
import numpy as np
import scipy.constants as spc

#-----------------------------------------------------------------------------------------------------------------------
# Brian2 import: we use Brian CPP-standalone code generation for fast parallelized simulations
#-----------------------------------------------------------------------------------------------------------------------
from brian2 import *
code_dir = './codegen'
prefs.GSL.directory = '/opt/anaconda3/envs/Brian2_NGILab/include/'   ## The directory where the GSL library headings are found
set_device('cpp_standalone',directory=code_dir,build_on_run=False)
prefs.devices.cpp_standalone.openmp_threads = 2 ## The number of threads used in the parallelization (machine-dependent)
prefs.logging.file_log = False
prefs.logging.delete_log_on_exit = True

import matplotlib.pyplot as plt

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
    return spc.R*(T+273.15)/spc.physical_constants['Faraday constant'][0]
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
                                                                  'include_dirs': ['/usr/include/gsl']})

def lpc5_parameters(model='hh-neuron',**kwargs):
    """
    Method to generate neuron parameters' dictionary.

    Input parameters:
    - model  : {'hh-neuron'} (or 'hh') | 'hh-ecs'   Model equations
    **kwargs : list of custom parameters specified by <parameter name>=<value in units>

    Return:
    - pars   : dict
    """

    # Define dictionary of default parameters values
    pars = {## Concentrations to setup reverse potentials
            'N0_i': 10*mmolar,
            'N0_e': 145*mmolar,
            'K0_i': 130*mmolar,
            'K0_e': 3*mmolar,
            'C0_i': 5*mmolar,
            'C0_e': 130*mmolar,
            ## Neuron Parameters and conductances
            'c_m' : 1*ufarad/cm**2,
            'g_Na': 2.04e6*usiemens/cm**2,
            'g_K' : 0.638e6*usiemens/cm**2,
            'g_Cl': 0.338e3*usiemens/cm**2,
            'v_thr' : 0*mvolt,
            ## Gating variables
            'U_m' : -38*mvolt,
            'U_h' : -66*mvolt,
            'U_n' : 18.7*mvolt,
            'W_m' : 6*mvolt,
            'W_h' : 6*mvolt,
            'W_n' : 9.7*mvolt,
            ## Temperature
            'T_exp' : 37, ## Body temperature of the animal
            ## External Stimulation
            'I_dc'  : 0*namp/cm**2,
            'I_ramp': 0*namp/cm**2,
            'T_ramp': 100*second,
            }

    pars = varargin(pars,**kwargs)

    # We define this quantities after any user-defined parameters are specified
    pars['T_adj'] = 2.3**(0.1*(pars['T_exp']-21))
    pars['E_Na'] = NernstPotential(pars['N0_e'],pars['N0_i'],1,pars['T_exp'])*volt
    pars['E_K']  = NernstPotential(pars['K0_e'],pars['K0_i'],1,pars['T_exp'])*volt
    pars['E_Cl'] = NernstPotential(pars['C0_e'],pars['C0_i'],-1,pars['T_exp'])*volt

    ## Physics constants
    pars['F'] = spc.physical_constants['Faraday constant'][0]*coulomb/mole

    return pars

def lpc5_neuron(N,params,model='hh-neuron',name='hh*',dt=None):
    """
    Method that generate the "neuron" model with given parameters.

    Input aarguments:
    - N      : int      Number of neurons to simulate
    - params : dict     Dictionary of neuron parameters
    - model  : {'hh-neuron'} | 'hh-ecs'  Model equations
    - name   : String   The name assigned to the NeuronGroup
    - dt     : float*second minimum time step for numerical integration

    Return:
    - neuron : NeuronGroup
    """

    eqs = Equations('''
    a_m=-0.182/second/mvolt*(v-U_m)/expm1((U_m-v)/W_m)   : hertz
    b_m=-0.124/second/mvolt*(U_m-v)/expm1((v-U_m)/W_m)   : hertz
    a_h=-0.015/second/mvolt*(U_h-v)/expm1((v-U_h)/W_h)   : hertz
    b_h=-0.015/second/mvolt*(v-U_h)/expm1((U_h-v)/W_h)   : hertz
    m_inf=a_m/(a_m+b_m)                                : 1
    h_inf=a_h/(a_h+b_h)                                : 1
    n_inf=1/(1+exp(-(v-U_n)/W_n))                      : 1
    tau_m=1e-3/(a_m+b_m)/T_adj                         : second
    tau_h=1e-3/(a_h+b_h)/T_adj                         : second
    tau_n=4*ms/(1+exp(-(v+56.56*mvolt)/(44.14*mvolt)))/T_adj      : second
    I_inj=I_dc+I_ramp*t/T_ramp                      : amp/meter**2
    dm/dt=(m_inf-m)/tau_m                           : 1
    dh/dt=(h_inf-h)/tau_h                           : 1
    dn/dt=(n_inf-n)/tau_n                           : 1
    ''')

    if model=='hh-neuron':
        eqs += Equations('''
            I_Na=g_Na*m**3*h*(v-E_Na)                       : amp/meter**2
            I_K=g_K*n*(v-E_K)                               : amp/meter**2
            I_Cl=g_Cl*(v-E_Cl)                              : amp/meter**2
            dv/dt=(I_inj-I_Na-I_K-I_Cl)/c_m                 : volt   
        ''')
    else:
        ## w/ ECS
        pass

    neuron = NeuronGroup(N,eqs,
                         threshold='v>v_thr',
                         namespace=params,
                         name=name,
                         method='rk4',
                         dt=dt)

    # Initialize gating variables
    neuron.m = 0.01
    neuron.h = 0.99
    neuron.n = 0.01

    return neuron

def lpc5_simulation(duration=1,model='hh-neuron',**kwargs):
    """
    This is the actual simulation of the Neuron. It is provided as a standalone method for modularity.

    Input parameters:
    - duration   : float*second   Simulated Time
    - model      : {'hh-neuron'} | 'hh-ecs'   Use only 'hh-neuron' for now. Other cases will be added later.
    :param kwargs:
    :return:
    """

    device.delete(force=True)
    ## Build the neuron model and monitors
    params = lpc5_parameters(model=model,T_ramp=duration*second,**kwargs)

    # Initialize model
    cell = lpc5_neuron(1,params,model=model,name='HH',dt=0.1*us)
    cell.v = NernstPotential(params['C0_e'],params['C0_i'],-1,params['T_exp'])*volt

    # Set monitors
    if model=='hh-neuron':
        variables = ['v']
    elif model=='hh-ecs':
        variables = ['v','C_i','N_i','K_i','I_K','I_Na','I_Cl','n_Cl','n_K','E_Cl']
    sv_mon = StateMonitor(cell,variables=variables,record=True,dt=0.1*ms,name='svmon')
    # Gene
    network = Network([cell,sv_mon])
    ## Run the simulator
    network.run(duration=duration*second,report='text')
    device.build(directory=code_dir, clean=True)

    ## Visualizing data
    if model=='hh-neuron':
        fig, ax = plt.subplots(1, 1)
        ax.plot(sv_mon.t, sv_mon.v[:].T, 'k-')
    elif model=='hh-ecs':
        fig, ax = plt.subplots(6, 1,sharex=True)
        ax[0].plot(sv_mon.t, sv_mon.v[:].T, 'k-')
        ax[0].plot(sv_mon.t, sv_mon.E_Cl[:].T, 'm-')
        # ax[1].plot(sv_mon.t, sv_mon.C_i[:].T, 'm-')
        # ax[1].plot(sv_mon.t, sv_mon.C_i[:].T, 'r-')
        # ax[1].plot(sv_mon.t, sv_mon.N_i[:].T, 'r-')
        # ax[1].plot(sv_mon.t,sv_mon.I_Na[:].T,'r-')
        # ax[2].plot(sv_mon.t,sv_mon.I_K[:].T,'b-')
        # ax[3].plot(sv_mon.t,sv_mon.I_Cl[:].T,'g-')
        # ax[4].plot(sv_mon.t,sv_mon.I_NKP[:].T,'y-')
        # ax[5].plot(sv_mon.t,sv_mon.I_KCC[:].T,'m-')
        ax[1].plot(sv_mon.t, sv_mon.n_Cl[:].T, 'm-')
        ax[1].hlines(0,*sv_mon.t[[0,-1]])
        ax[2].plot(sv_mon.t, sv_mon.I_Cl[:].T, 'r-')
        ax[3].plot(sv_mon.t, sv_mon.n_K[:].T, 'g-')

    device.delete(force=True)

if __name__=="__main__":
    """
    Uncomment individual sections to test different methods 
    """

    # #-------------------------------------------------------------------------------------------------------------------
    # # Verify Reverse Potentials
    # #-------------------------------------------------------------------------------------------------------------------
    # T= 37
    # print('E_Na\t',NernstPotential(126,9,1,T))
    # print('E_K\t',NernstPotential(2.5,130,1,T))
    # print('E_Cl\t',NernstPotential(130,10,-1,T))

    # #-------------------------------------------------------------------------------------------------------------------
    # # Generate Neuron Simulation
    # #-------------------------------------------------------------------------------------------------------------------
    print(lpc5_parameters('hh-neuron',I_dc=10**4.5*nA/cm**2))
    # Single Neuron
    lpc5_simulation(duration=1,model='hh-neuron',I_dc=10**4.2*nA/cm**2)

    # #-------------------------------------------------------------------------------------------------------------------
    ## Show figures
    # #-------------------------------------------------------------------------------------------------------------------
    plt.show()