# Equations.py
import math
import scipy.constants as spc
import sympy as sp
from brian2 import *


def return_HH_equations(model = 'hh-neuron'):
    
    eqs = Equations('''
    x       : 1 (constant)
    y       : 1 (constant)
    z       : 1 (constant)
    a_m=-0.182/second/mvolt*(v-U_m)/(expm1((U_m-v)/W_m))   : hertz
    b_m=-0.124/second/mvolt*(U_m-v)/(expm1((v-U_m)/W_m))   : hertz
    a_h=-0.015/second/mvolt*(U_h-v)/(expm1((v-U_h)/W_h))   : hertz
    b_h=-0.015/second/mvolt*(v-U_h)/(expm1((U_h-v)/W_h))   : hertz
    m_inf=a_m/(a_m+b_m)                                : 1
    h_inf=a_h/(a_h+b_h)                                : 1
    n_inf=1/(1+exp(-(v-U_n)/W_n))                      : 1
    tau_m=1e-3/(a_m+b_m)/T_adj                         : second
    tau_h=1e-3/(a_h+b_h)/T_adj                         : second
    tau_n=4*ms/(1+exp(-(v+56.56*mvolt)/(44.14*mvolt)))/T_adj      : second
    mask                                            : 1
    mask_noise                                      : 1
    I_inj=mask * I_max * TimedArray(t)                   : amp/meter**2
    dI_noise/dt = mask_noise * sigma_noise * sqrt(c_m/g_noise_L) * xi/second           : amp/meter**2
    dm/dt=(m_inf-m)/tau_m                           : 1
    dh/dt=(h_inf-h)/tau_h                           : 1
    dn/dt=(n_inf-n)/tau_n                           : 1
    ''')

    if model == 'hh-neuron':
        eqs += Equations('''
            I_Na=g_Na*m**3*h*(v-E_Na)                       : amp/meter**2
            I_K=g_K*n*(v-E_K)                               : amp/meter**2
            I_Cl_L=g_Cl_L*(v-E_Cl)                              : amp/meter**2
            dv/dt=(I_inj-I_Na-I_K-I_Cl_L + I_noise)/c_m                 : volt
        ''')
    
    elif model=='hh-ecs':
        eqs += Equations('''
            dn_Na_N/dt = -S_N/F*(I_Na + I_Na_L + 3*I_NKP) : mol
            dn_Cl_N/dt = S_N/F*(I_Cl_L + I_KCC) : mol
            dn_K_N/dt = -S_N/F*(I_K + I_K_L - 2*I_NKP - I_KCC) : mol
            n_Na_E                              : mol
            n_K_E                               : mol
            n_Cl_E                              : mol                  
            C_Na_N = clip(C0_Na_N + n_Na_N/V_N, 0 * mM, 155 * mM) : mM
            C_Na_E = clip(C0_Na_E - n_Na_N/V_E, 0 * mM, 155 * mM) : mM
            C_Cl_N = clip(C0_Cl_N + n_Cl_N/V_N, 0 * mM, 133 * mM) : mM
            C_Cl_E = clip(C0_Cl_E - n_Cl_N/V_E, 0 * mM, 133 * mM) : mM
            C_K_N = clip(C0_K_N + n_K_N/V_N, 0 * mM, 135 * mM) : mM
            C_K_E = clip(C0_K_E - n_K_N/V_E, 0 * mM, 135 * mM)  : mM
            E_Na = NernstPotential(C_Na_E, C_Na_N, 1, T_exp)*volt : volt
            E_K = NernstPotential(C_K_E, C_K_N, 1, T_exp)*volt : volt
            E_Cl = NernstPotential(C_Cl_E, C_Cl_N, -1, T_exp)*volt : volt    
            g_Na_mod                                        : 1
            g_K_mod                                        : 1
            I_Na = g_Na_mod * g_Na*m**3*h*(v-E_Na)                       : amp/meter**2      
            I_K = g_K_mod * g_K*n*(v-E_K)                               : amp/meter**2
            I_Cl_L = g_Cl_L*(v-E_Cl)                              : amp/meter**2 
            I_Na_L = g_Na_L * (v-E_Na)                          : amp/meter**2 
            I_K_L = g_K_L * (v-E_K)                         : amp/meter**2 
            I_L = I_Na_L + I_Cl_L + I_K_L                     : amp/meter**2 
            sigma = 1/7 * expm1(C_Na_E/(67.3 * mM))              : 1
            f_NaK = 1/(1 + 0.1245 * exp(-0.1 * v/V_T) + 0.0365 * sigma * exp(-v/V_T)) : 1
            I_NKP = I_NKP_max * f_NaK * Hill(C_Na_N, zeta_Na, 1.5) * Hill(C_K_E, zeta_K, 1) : amp/meter**2
            I_KCC = g_KCC*(E_K - E_Cl)                          : amp/meter**2           
            dv/dt=(I_inj - I_Na - I_Na_L - I_K - I_K_L - I_Cl_L - I_NKP)/c_m       : volt
            Injected_current = I_max                           :amp/meter**2
            Check_1 = I_Na + I_Na_L + 3 * I_NKP                          :amp/meter**2
            Check_2 = I_K + I_K_L - 2 * I_NKP - I_KCC                   :amp/meter**2
            Check_3 = I_Cl_L + I_KCC                              :amp/meter**2
        ''')
    
    
    # a_u=-0.0033/second * exp(0.1 * (v + 35*mvolt)/mvolt)   : hertz
    # b_u=0.0033/second * expm1(0.1* (v + 35*mvolt)/mvolt)   : hertz
    # u_inf=a_m/(a_m+b_m)                                : 1
    # tau_u=1/(T_adj*(a_u + b_u))                       : second
    # du/dt = (u_inf - u)/tau_u                         : 1
    # g_Kv7_mod                                      : 1
    # I_Kv7 = g_Kv7_mod * g_Kv7 * u * (v-E_K)                      : amp/meter**2

    elif model=='hh-ecs-synapse':
        eqs += Equations('''
            dn_Na_N/dt = -S_N/F*(I_Na + I_Na_L + 3*I_NKP) : mol
            dn_Cl_N/dt = S_N/F*(I_Cl_L + I_KCC) : mol
            dn_K_N/dt = -S_N/F*(I_K + I_K_L - 2*I_NKP - I_KCC) : mol
            C_Na_N = clip(C0_Na_N + n_Na_N/V_N, 0 * mM, 155 * mM) : mM
            C_Na_E = clip(C0_Na_E - n_Na_N/V_E, 0 * mM, 155 * mM) : mM
            C_Cl_N = clip(C0_Cl_N + n_Cl_N/V_N, 0 * mM, 133 * mM) : mM
            C_Cl_E = clip(C0_Cl_E - n_Cl_N/V_E, 0 * mM, 133 * mM) : mM
            C_K_N = clip(C0_K_N + n_K_N/V_N, 0 * mM, 135 * mM) : mM
            C_K_E = clip(C0_K_E - n_K_N/V_E, 0 * mM, 135 * mM)  : mM
            E_Na = NernstPotential(C_Na_E, C_Na_N, 1, T_exp)*volt : volt
            E_K = NernstPotential(C_K_E, C_K_N, 1, T_exp)*volt : volt
            E_Cl = NernstPotential(C_Cl_E, C_Cl_N, -1, T_exp)*volt : volt          
            I_Na = g_Na*m**3*h*(v-E_Na)                       : amp/meter**2      
            I_K = g_K*n*(v-E_K)                               : amp/meter**2
            I_Cl_L = g_Cl_L*(v-E_Cl)                              : amp/meter**2 
            I_Na_L = g_Na_L * (v-E_Na)                          : amp/meter**2 
            I_K_L = g_K_L * (v-E_K)                         : amp/meter**2 
            I_L = I_Na_L + I_Cl_L + I_K_L                     : amp/meter**2 
            sigma = 1/7 * expm1(C_Na_E/(67.3 * mM))              : 1
            f_NaK = 1/(1 + 0.1245 * exp(-0.1 * v/V_T) + 0.0365 * sigma * exp(-v/V_T)) : 1
            I_NKP = I_NKP_max * f_NaK * Hill(C_Na_N, zeta_Na, 1.5) * Hill(C_K_E, zeta_K, 1) : amp/meter**2
            I_KCC = g_KCC*(E_K - E_Cl)                          : amp/meter**2           
            dv/dt=(I_inj - I_Na - I_Na_L - I_K - I_K_L - I_Cl_L - I_NKP - I_syn + I_noise)/c_m       : volt
            g_syn_e                               : siemens/meter**2 
            g_syn_i                               : siemens/meter**2 
            E_AMPA = V_T * log((C_Na_E + 1.2*C_K_E)/(C_Na_N + 1.2 * C_K_N)) : volt
            E_GABA = -150 *mvolt                            : volt
            I_syn = g_syn_e * (v - E_AMPA) + g_syn_i * (v - E_GABA)  : amp/meter**2 
            Check_1 = I_Na + I_Na_L + 3 * I_NKP                          :amp/meter**2
            Check_2 = I_K + I_K_L - 2 * I_NKP - I_KCC                   :amp/meter**2
            Check_3 = I_Cl_L + I_KCC                              :amp/meter**2
        ''') 
    ## I_syn = g_syn*(v - E_syn)                       : amp/meter**2 
    ##  E_syn = V_T * log((C_Na_E + 1.2*C_K_E)/(C_Na_N + 1.2 * C_K_N)) : volt
    ## E_syn has to be E_AMPA (E_syn = V_T * log((N0_e + 1.2*K0_e)/(N0_i + 1.2 * K0_i)) : volt) for excitatory synapse and E_GABA () = E_Cl for inhibitory synapse
  
    elif model=='hh-ecs-synapse-astrocyte':
        eqs += Equations('''
            dn_Na_N/dt = -S_N/F*(I_Na + I_Na_L + 3*I_NKP) : mol
            dn_Cl_N/dt = S_N/F*(I_Cl_L + I_KCC) : mol
            dn_K_N/dt = -S_N/F*(I_K + I_K_L - 2*I_NKP - I_KCC) : mol
            C_Na_N = clip(C0_Na_N + n_Na_N/V_N, 0 * mM, 155 * mM) : mM
            C_Na_E = clip(C0_Na_E - n_Na_N/V_E, 0 * mM, 155 * mM) : mM
            C_Cl_N = clip(C0_Cl_N + n_Cl_N/V_N, 0 * mM, 133 * mM) : mM
            C_Cl_E = clip(C0_Cl_E - n_Cl_N/V_E, 0 * mM, 133 * mM) : mM
            C_K_N = clip(C0_K_N + n_K_N/V_N, 0 * mM, 135 * mM) : mM
            C_K_E = clip(C0_K_E - n_K_N/V_E, 0 * mM, 135 * mM)  : mM
            E_Na = NernstPotential(C_Na_E, C_Na_N, 1, T_exp)*volt : volt
            E_K = NernstPotential(C_K_E, C_K_N, 1, T_exp)*volt : volt
            E_Cl = NernstPotential(C_Cl_E, C_Cl_N, -1, T_exp)*volt : volt          
            I_Na = g_Na*m**3*h*(v-E_Na)                       : amp/meter**2      
            I_K = g_K*n*(v-E_K)                               : amp/meter**2
            I_Cl_L = g_Cl_L*(v-E_Cl)                              : amp/meter**2 
            I_Na_L = g_Na_L * (v-E_Na)                          : amp/meter**2 
            I_K_L = g_K_L * (v-E_K)                         : amp/meter**2 
            I_L = I_Na_L + I_Cl_L + I_K_L                     : amp/meter**2 
            sigma = 1/7 * expm1(C_Na_E/(67.3 * mM))              : 1
            f_NaK = 1/(1 + 0.1245 * exp(-0.1 * v/V_T) + 0.0365 * sigma * exp(-v/V_T)) : 1
            I_NKP = I_NKP_max * f_NaK * Hill(C_Na_N, zeta_Na, 1.5) * Hill(C_K_E, zeta_K, 1) : amp/meter**2
            I_KCC = g_KCC*(E_K - E_Cl)                          : amp/meter**2           
            dv/dt=(I_inj - I_Na - I_Na_L - I_K - I_K_L - I_Cl_L - I_NKP - I_syn + I_noise)/c_m       : volt
            g_syn_e                               : siemens/meter**2 
            g_syn_i                               : siemens/meter**2 
            E_AMPA = V_T * log((C_Na_E + 1.2*C_K_E)/(C_Na_N + 1.2 * C_K_N)) : volt
            E_GABA = -150 *mvolt                            : volt
            I_syn = g_syn_e * (v - E_AMPA) + g_syn_i * (v - E_GABA)  : amp/meter**2 
            I_NKP_A                                                 : amp/meter**2
            I_Kir_A                                              : amp/meter**2
            I_Na_L_A                                           : amp/meter**2
            Check_1 = I_Na + I_Na_L + 3 * I_NKP                          :amp/meter**2
            Check_2 = I_K + I_K_L - 2 * I_NKP - I_KCC                   :amp/meter**2
            Check_3 = I_Cl_L + I_KCC                              :amp/meter**2
        ''') 

    else:
        pass

    return eqs

def return_astrocyte_equations(model):
    if model == 'hh-ecs-synapse':
        eqs = Equations('''
        n_Na_A                                              : mol
        n_K_A                                               : mol
        n_Cl_A                                              : mol
        C_Na_A                                              : mM
        C_K_A                                               : mM
        C_Cl_A                                              : mM
        E_Na_A                                              : volt
        E_K_A                                               : volt
        E_Cl_A                                              : volt
        v_A                                                 : volt
        ''')
    
    return

# elif model == 'hh-neuron-synapse':
#         eqs += Equations('''
#             I_Na=g_Na*m**3*h*(v-E_Na)                       : amp/meter**2
#             I_K=g_K*n*(v-E_K)                               : amp/meter**2
#             I_Cl_L=g_Cl*(v-E_Cl)                            : amp/meter**2
#             dv/dt=(I_inj-I_Na-I_K-I_Cl_L-I_syn)/c_m         : volt   
#             g_syn                               : siemens/meter**2
#             E_syn = E_Cl : volt
#             I_syn = g_syn*(v - E_syn)                       : amp/meter**2
#         ''')
#         ## E_syn has to be E_AMPA (E_syn = V_T * log((N0_e + 1.2*K0_e)/(N0_i + 1.2 * K0_i)) : volt) for excitatory synapse and E_GABA () for inhibitory synapse

def return_plotting_list(model):
    ## To guarantee an error free code execution make sure that the last list item contains the highest number out of all plot_number keys
    
    if model == 'hh-ecs':
        plotting_list = [
            {
                'variable': 'I_inj', 
                'axis': 'I_inj in nA/cm^2',
                'plot_number' : 0,
                'unit' : namp/cm**2
            },
            {
                'variable': 'I_Na',
                'axis': 'I_Na in A/m^2',
                'plot_number' : 1,
                'unit': namp/cm**2
            },
            {
                'variable': 'I_Na_L',
                'axis': 'I_Na_L in A/m^2',
                'plot_number' : 2,
                'unit': namp/cm**2
            },
            {
                'variable': 'I_K', 
                'axis': 'I_K (nA/m^2)',
                'plot_number' : 3,
                'unit': namp/cm**2
            },
            {
                'variable': 'I_Cl_L', 
                'axis': 'I_Cl (nA/m^2)',
                'plot_number' : 4,
                'unit': namp/cm**2 
            },
            {
                'variable': 'I_NKP',
                'axis': 'I_NKP in nA/cm^2',
                'plot_number' : 5,
                'unit': namp/cm**2
            },
            {
                'variable': 'I_KCC',
                'axis': 'I_KCC in a/m^2',
                'plot_number' : 6,
                'unit': namp/cm**2
            },
            {
                'variable': 'v', 
                'axis': 'v (mV)',
                'plot_number' : 7,
                'unit' : mV
            },
            {
                'variable': 'E_K',
                'axis': 'E_K (mV)',
                'plot_number' : 8,
                'unit': mV
            },
            {
                'variable': 'E_Cl',
                'axis': 'E_Cl (mV)',
                'plot_number' : 9,
                'unit': mV
            },
            {
                'variable': 'E_Na',
                'axis': 'E_Na (mV)',
                'plot_number' : 10,
                'unit': mV
            },
            {
                'variable': 'C_Na_N',
                'axis': 'C_Na_N in mol/m^3',
                'plot_number' : 11,
                'unit': mmolar
            },
            {
                'variable': 'C_Na_E',
                'axis': 'C_Na_E in mol/m^3',
                'plot_number' : 12,
                'unit': mmolar
            },
            {
                'variable': 'C_K_N', 
                'axis': 'C_K_N in mol/m^3',
                'plot_number' : 13,
                'unit': mmolar
            },
            {
                'variable': 'C_K_E', 
                'axis': 'C_K_E in mol/m^3',
                'plot_number' : 14,
                'unit': mmolar
            },
            {
                'variable': 'C_Cl_N',
                'axis': 'C_Cl_N in mol/m^3',
                'plot_number' : 15,
                'unit': mmolar
            },
            {
                'variable': 'C_Cl_E',
                'axis': 'C_Cl_E in mol/m^3',
                'plot_number' : 16,
                'unit': mmolar
            },
            {
                'variable': 'Check_1',
                'axis': 'Check_1 in nA/cm^2',
                'plot_number' : 17,
                'unit': namp/cm**2
            },
            {
                'variable': 'Check_2',
                'axis': 'Check_2 in nA/cm^2',
                'plot_number' : 18,
                'unit': namp/cm**2
            },
            {
                'variable': 'Check_3',
                'axis': 'Check_3 in nA/cm^2',
                'plot_number' : 19,
                'unit': namp/cm**2
            },
            {
                'variable': 'I_Kv7',
                'axis': 'I_Kv7 (nA/cm^2)',
                'plot_number' : 20,
                'unit': namp/cm**2
            }]
            # },
            # {
            #     'variable': 'I_syn',
            #     'axis': 'I_syn (nA/cm^2)',
            #     'plot_number' : 20,
            #     'unit': namp/cm**2
            # }]
        
    # {
    #             'variable': 'E_syn',
    #             'axis': 'E_Syn (mV)',
    #             'plot_number' : 22,
    #             'unit': mV
    #         }

    # {
    #             'variable': 'g_syn',
    #             'axis': 'g_syn (usiemens/cm^2)',
    #             'plot_number' : 21,
    #             'unit': usiemens/cm**2
    #         }

    # elif model == 'hh-neuron':
    #     plotting_list = [
    #         {
    #             'variable': 'v', 
    #             'axis': 'v (mV)',
    #             'plot_number' : 0,
    #             'unit' : mvolt
    #         },
    #         {
    #             'variable': 'm', 
    #             'axis': 'gating variable m',
    #             'plot_number' : 1,
    #             'unit' : 1
    #         },
    #         {
    #             'variable': 'n', 
    #             'axis': 'gating variable n',
    #             'plot_number' : 2,
    #             'unit' : 1
    #         },
    #         {
    #             'variable': 'h', 
    #             'axis': 'gating variable h',
    #             'plot_number' : 3,
    #             'unit' : 1
    #         }
    elif model == 'hh-neuron':
        plotting_list = [
            {
                'variable': 'v', 
                'axis': 'v (mV)',
                'plot_number' : 0,
                'unit' : mvolt
            },
            {
                'variable': 'I_noise', 
                'axis': 'namp/cm**2',
                'plot_number' : 1,
                'unit' : namp/cm**2
            }
        ]
    return plotting_list

