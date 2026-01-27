## NEA params.py
from brian2 import *
from brian2.units.constants import faraday_constant as F
from brian2.units.constants import avogadro_constant as N_A



#-----------------------------------------------------------------------------------------------------------------------
## General-Purpose Utilities
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
## Generate Model Parameters
#-----------------------------------------------------------------------------------------------------------------------
def nea_parameters(**kwargs):
    pars_neu = {## Concentrations to setup reverse potentials
            'N0_i'   : 10*mmolar,
            'K0_i'   : 130*mmolar,
            'C0_i'   : 5*mmolar,
            'HBC0_i' : 15*mmolar,  # The HBC internally must be such that v_rest > E_GABA (so GABA is always inhibitory for the neuron)
            ## Neuron Parameters and conductances
            'c_m'    : 1*ufarad/cm**2,
            'g_Na'   : 2.04e6*usiemens/cm**2,
            'g_K'    : 0.693e6*usiemens/cm**2,
            'g_L_Na' : 32.7*usiemens/cm**2,
            'g_L_K'  : 0.0*usiemens/cm**2,
            'g_L_Cl' : 50*usiemens/cm**2,
            'v_thr'  : 0*mvolt,
            ## Gating variables
            'U_m'    : -38*mvolt,
            'U_h'    : -66*mvolt,
            'U_n'    : 18.7*mvolt,
            'W_m'    : 6*mvolt,
            'W_h'    : 6*mvolt,
            'W_n'    : 9.7*mvolt,
            ## Temperature
            'T_exp'  : 37, ## Body temperature of the animal
            ## External Stimulation
            'I_dc'   : 0*namp/cm**2,
            'I_ramp' : 0*namp/cm**2,
            'T_ramp' : 100*second,
            ## Geometry
            'S_n'      : 700*um**2,
            'Lambda_n' : 1750*um**3,
            ## Transport
            'I_NKA'  : 0*namp/cm**2,
            'zeta_Na': 13*mmolar,  ## Can be between 25--30 Clausen et al., Front Physiol. 2017
            'zeta_K' : 0.2*mmolar,  ## Clausen et al., Front Physiol. 2017
            ## I_KCC
            'g_KCC'  : 0*usiemens/cm**2,
            ## Permeability Ratios
            'P_K_Na' : 1.2,
            'P_HBC_Cl' : 0.4,
            }

    pars_ecs = {# Extracellular concentrations
            'N0_e'   : 145*mmolar,
            'K0_e'   : 3*mmolar,
            'C0_e'   : 130*mmolar,
            'HBC0_e'  : 35*mmolar,
            'H0_e'    : 50*nmolar,
            'G0_e'    : 25*nmolar,
            'GABA0_e' : 50*nmolar,
            # Diffusion Rates
            'D_Na_e' : 2/second,
            'D_K_e'  : 2/second,
            'D_Cl_e' : 2/second,
            'D_Glu_e'  : 5/second,
            'D_GABA_e' : 5/second,
            'D_Glu'  : 0.33*um**2/msecond,
            'D_GABA' : 0.33*um**2/msecond,
            # Geometry
            'Lambda_e' : 500*um**3,
            # Synaptic geometry
            'l_diff' : 0.13*um,
            't_cleft': 0.2*um
    }

    pars_astro = pars = {## Concentrations to setup reverse potentials
            'N0_a'    : 15*mmolar,
            'K0_a'    : 100*mmolar,
            'C0_a'    : 40*mmolar,
            'H0_a'    : 60*nmolar,
            'G0_a'    : 25*nmolar, ## Glutamate
            'GABA0_a' : 20*mmolar, ## GABA ()
            'HBC0_a'  : 10*mmolar,
            ## Astrocyte Parameters
            'c_m_a'   : 1*ufarad/cm**2,
            ## Temperature
            'T_exp'   : 37, ## Body temperature of the animal
            ## External Stimulation
            'I_dc'    : 0*namp/cm**2,
            'I_ramp'  : 0*namp/cm**2,
            'T_ramp'  : 100*second,
            ## Geometry
            'S_a'       : 850*um**2,
            'Lambda_a'  : 1000*um**3,
            ## Kir
            'g_Kir'     : 175*usiemens/cm**2,
            'zeta_Kir'  : 13*mmolar,
            ## EAAT
            'sigma_EAAT': 100/um**2,
            'g_EAAT'    : 100*usiemens/cm**2,
            'g_T_Cl_a'  : 100*usiemens/cm**2,
            'g_L_Cl_a'  : 60*usiemens/cm**2,
            'Omega_Glu' : 25/second,
            ## GAT
            'g_GAT'   : 6*usiemens/cm**2,
            ## NKCC
            'g_NKCC'  : 0*usiemens/cm**2,
            ## NKP
            'I_NKA_a' : 0*namp/cm**2,
            ## GABA
            'g_GABA'  : 0*usiemens/cm**2,
            'tau_GABA': 50*ms,
            'J_GABA'  : 1/umolar/second,
            ## Buffering
            'D_Na_a'  : 0.1/second,
            'D_K_a'   : 0.1/second,
            'D_Cl_a'  : 0.1/second,
            'D_Na_m'  : 0.01/second,
            'D_K_m'   : 0.1/second,
            'D_Cl_m'  : 0.01/second,
    }

    # Generate default dictionary
    pars = {**pars_neu,**pars_ecs,**pars_astro}

    # Custom-parameters
    pars = varargin(pars,**kwargs)

    # Retrieve partial permeabilities
    pars['pi_AMPA_Na'] = 1.0/(1+pars['P_K_Na'])
    pars['pi_AMPA_K'] = pars['P_K_Na']/(1+pars['P_K_Na'])
    pars['pi_GABA_Cl'] = 1.0/(1+pars['P_HBC_Cl'])

    # Extrapolate useful variables
    pars['T_adj'] = 2.3**(0.1*(pars['T_exp']-21))

    # Extrapolation of the diffusion conductances [volume/sec]
    pars['k_Glu'] = np.pi*pars['D_Glu']*(pars['l_diff']+2*pars['t_cleft'])
    pars['k_GABA'] = np.pi*pars['D_GABA']*(pars['l_diff']+2*pars['t_cleft'])

    return pars



def synapse_parameters(ttype='glu',**kwargs):
    assert any(ttype==t for t in ['glu','gaba']),"Allowed transmitter types (ttype): 'glu' or 'gaba' only"
    if ttype=='glu':
        pars = {'Nt_rel' : 0.1*mmolar,
                'J'      : 1/umolar/second,
                'tau_r'  : 10*msecond,
                'g'      : 10*nsiemens/cm**2,
                'D_Glu'  : 10/second,
        }
    elif ttype=='gaba':
        pars = {'Nt_rel' : 1.0*mmolar,
                'J'      : 1/umolar/second,
                'tau_r'  : 50*msecond,
                'g'      : 10*nsiemens/cm**2,
                'D_GABA' : 10/second,
                }

    pars['Lambda_s'] = 8.75e-3*um**3
    # The current S_s and sigma_R are temporarily taken from the reference below --> need to be better estimated
    # https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0130924
    pars['S_s'] = 200*nmeter
    pars['sigma_R'] = 1000*um**-2
    # Estimate receptor density
    pars['R_T'] = pars['sigma_R']/N_A/pars['Lambda_s']*pars['S_s']**2

    pars = varargin(pars,**kwargs)

    return pars