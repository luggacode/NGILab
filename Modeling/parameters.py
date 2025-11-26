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


def return_initial_parameters(model='hh-neuron',**kwargs):
    """
    Method to generate neuron parameters' dictionary.

    Input parameters:
    - model  : {'hh-neuron'} (or 'hh') | 'hh-ecs'   Model equations
    **kwargs : list of custom parameters specified by <parameter name>=<value in units>

    Return:
    - pars   : dict
    """

    ta = TimedArray([0] * 2 + [1] * 2 + [0] * 16, dt=0.005*second)

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
            'g_Na': 2.04e6*usiemens/cm**2,
            'g_Na_L': NaN *usiemens/cm**2,
            'g_K_L': 0*usiemens/cm**2,
            'g_K' : 0.693e6*usiemens/cm**2,
            'g_Cl_L': 0.5e2*usiemens/cm**2,
            'g_KCC' : NaN*usiemens/cm**2,
            'v_thr' : 59.3295*mvolt,
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
            'I_max' : 4000*namp/cm**2,
            'T_ramp': 100*second,
            'TimedArray' : ta,
            ## times
            't_step': 0.05 * second,
            't_end' : 0.1 * second,
            ## Cell surface and volume
            'S_N' : 150*um**2,
            'V_N' : 500*um**3,
            'V_E' : 500*um**3,
            # zetas
            'zeta_Na' : 13*mmolar,
            'zeta_K' : 0.2*mmolar, 
            }

    # pars = varargin(pars,**kwargs)

    # We define this quantities after any user-defined parameters are specified
    pars['T_adj'] = 2.3**(0.1*(pars['T_exp']-21))
    pars['E_Na'] = NernstPotential(pars['C0_Na_E'],pars['C0_Na_N'],1,pars['T_exp'])*volt
    pars['E_K']  = NernstPotential(pars['C0_K_E'],pars['C0_K_N'],1,pars['T_exp'])*volt
    pars['E_Cl'] = NernstPotential(pars['C0_Cl_E'],pars['C0_Cl_N'],-1,pars['T_exp'])*volt
    pars['g_KCC'] = pars['g_Cl_L'] * (-70*mvolt - pars['E_Cl'])/(pars['E_Cl'] - pars['E_K'])

    ## Physics constants
    pars['F'] = spc.physical_constants['Faraday constant'][0]*coulomb/mole
    pars['k'] = spc.physical_constants['Boltzmann constant'][0]*joule/kelvin
    pars['q'] = spc.physical_constants['elementary charge'][0]*coulomb
    pars['e'] = math.e
    
    # Thermal Voltage
    pars['V_T'] = ThermalVoltage(37) * volt


    return pars

return_initial_parameters()