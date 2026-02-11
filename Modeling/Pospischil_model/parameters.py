from brian2 import *
import math
import scipy.constants as spc


def Hill(x,K,n):
    return x**n/(x**n+K**n)
def ThermalVoltage(T):
    return 8.314472000000000e+00*(T+273.15)/9.648534297750000e+04
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


def return_initial_parameters(model):
    """
    Method to generate neuron parameters' dictionary.

    Input parameters:
    - model  : {'hh-neuron'} (or 'hh') | 'hh-ecs'   Model equations
    **kwargs : list of custom parameters specified by <parameter name>=<value in units>

    Return:
    - pars   : dict
    """
    ta = TimedArray([1] * 5, dt= 1 *second)
    # ta = TimedArray([0] * 2 + [1] * 2 + [0] * 16, dt=0.005*second)

    # Define dictionary of default parameters values
    pars = {## Concentrations to setup reverse potentials
            'C0_Na_N': 10*mmolar,
            'C0_Na_E': 145*mmolar,
            'C0_K_N': 130*mmolar,
            'C0_K_E': 3*mmolar,
            'C0_Cl_N': 5*mmolar,
            'C0_Cl_E': 130*mmolar,
            ## Neuron Parameters and conductances
            'c_m' : 1*ufarad/cm**2,
            'g_noise_L' : 0.2e6 *usiemens/cm**2,
            'v_thr' : - 50 * mvolt,     
            'resting_potential' : -70,
            ## Gating variables
            'U_m' : -38*mvolt,
            'U_h' : -66*mvolt,
            'U_n' : 18.7*mvolt,
            'W_m' : 6*mvolt,
            'W_h' : 6*mvolt,
            'W_n' : 9.7*mvolt,
            ## Temperature
            'T_exp' : 36, ## Body temperature of the animal
            ## External Stimulation
            'I_ramp': 1000*namp/cm**2,
            'I_max' : 500*namp/cm**2,
            'T_ramp': 0.4*second,
            'TimedArray' : ta,
            ## Cell surface and volume
            'S_N' : 150*um**2,
            'V_N' : 500*um**3,
            'V_E' : 500*um**3,
            # zetas
            'zeta_Na' : 13*mmolar,
            'zeta_K' : 0.2*mmolar, 
            # sigma_noise 
            'sigma_noise' : 500000 * namp/cm**2,
            # V_T
            'V_T' : NaN * volt,
            'E_Na': NaN *volt,
            'E_K' : NaN * volt,
            'E_Cl' : NaN * volt,
            'E_L' : NaN * volt,
            'g_Na' : NaN *siemens/cm**2, 
            'g_Kd' : NaN * siemens/cm**2,
            'g_M' : NaN * siemens/cm**2,        ##parameters need to be checked to be fitting to one specific model, could be 4 could be 97 usiems/cm**2
            'g_L' : NaN * siemens/cm**2,
            'tau_p_max' : NaN * second,
            'g_Na_L' : NaN *usiemens/cm**2,
            'g_K_L' : NaN * usiemens/cm**2,
            'g_Cl_L' : NaN * usiemens/cm**2,
            'g_KCC' : NaN * siemens/cm**2
            }

    # pars = varargin(pars,**kwargs)

    if model == 'hh-neuron_exc' or model == 'hh-neuron-synapse_exc':
        pars['E_Na'] = 50 * mvolt
        pars['E_K']  = -90 * mvolt
        pars['E_L'] = -70.3 * mvolt
        pars['g_Na'] = 5.6e4*usiemens/cm**2
        pars['g_Kd'] = 6e3*usiemens/cm**2
        pars['g_M'] = 75 * usiemens/cm**2        ##parameters need to be checked to be fitting to one specific model, could be 4 could be 97 usiems/cm**2
        pars['g_L'] = 20.5 * usiemens/cm**2
        pars['V_T'] = -56.2 * mvolt
        pars['tau_p_max'] = 0.608 * second
    elif model == 'hh-neuron_inh' or model == 'hh-neuron-synapse_inh':
        pars['E_Na'] = 50 * mvolt
        pars['E_K']  = -90 * mvolt
        pars['E_L'] = -56.2 * mvolt
        pars['g_Na'] = 1e4*usiemens/cm**2
        pars['g_Kd'] = 2.1e3*usiemens/cm**2
        pars['g_M'] = 98 * usiemens/cm**2        ##parameters need to be checked to be fitting to one specific model, could be 4 could be 97 usiems/cm**2
        pars['g_L'] = 13.3 * usiemens/cm**2
        pars['V_T'] = -67.9 * mvolt
        pars['tau_p_max'] = 0.934 * second
    elif model == 'hh-ecs_exc' or model == 'hh-ecs-synapse_exc':
        pars['E_Na'] = NernstPotential(pars['C0_Na_E'],pars['C0_Na_N'],1,pars['T_exp'])*volt
        pars['E_K']  = NernstPotential(pars['C0_K_E'],pars['C0_K_N'],1,pars['T_exp'])*volt
        pars['E_Cl'] = NernstPotential(pars['C0_Cl_E'],pars['C0_Cl_N'],-1,pars['T_exp'])*volt
        pars['E_L'] = -70.3 * mvolt
        pars['g_Na_L'] = NaN *usiemens/cm**2
        pars['g_K_L'] = 0*usiemens/cm**2
        pars['g_Cl_L'] = 0.5e2*usiemens/cm**2
        pars['g_KCC'] = pars['g_Cl_L'] * (-70*mvolt - pars['E_Cl'])/(pars['E_Cl'] - pars['E_K'])
        pars['V_T'] = -56.2 * mvolt 

    # We define this quantities after any user-defined parameters are specified
    pars['T_adj'] = 2.3**(0.1*(pars['T_exp']-21))
    
    ## Physics constants
    pars['F'] = spc.physical_constants['Faraday constant'][0]*coulomb/mole
    pars['k'] = spc.physical_constants['Boltzmann constant'][0]*joule/kelvin
    pars['q'] = spc.physical_constants['elementary charge'][0]*coulomb
    pars['e'] = math.e
    
    # pars['V_T'] = ThermalVoltage(pars['T_exp']) * 1000 * mvolt
    
    
    
    return pars




def return_synapse_pars():
    pars = {
        'E_syn' : 50 * mvolt,
        'tau_r_x' : 5*ms,
        'tau_x' : 50*ms,
        'J' : 1,
        'g_s_inh' : 400*usiemens/cm**2,
        'g_s_exc' : 100*usiemens/cm**2,
        'W' : 4,
        'A' : 0.6,
        #'c_m': 1*ufarad/cm**2,
        #'weight': 60*mvolt,
    }
    return pars