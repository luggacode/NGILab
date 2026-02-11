# Equations.py
import  math
import scipy.constants as spc
import sympy as sp
from brian2 import *


def return_HH_equations(model = 'hh-neuron'):
    eqs = Equations('''
    x       : 1 (constant)
    y       : 1 (constant)
    z       : 1 (constant)
    a_m = -0.32 * (v - V_T - 13 * mvolt)/(exp(-(v - V_T - 13*mvolt)/(4*mvolt)) - 1)/(mvolt*ms)   : hertz
    b_m =  0.28*(v-V_T-40*mvolt)/(exp((v-V_T-40*mvolt)/(5*mvolt))-1)/(mvolt*ms) : hertz
    a_h = 0.128 * exp(-(v - V_T - 17*mvolt)/(18*mvolt))/ms                                    : hertz
    b_h = 4/(1 + exp(-(v- V_T - 40*mvolt)/(5*mvolt)))/ms                                      : hertz
    a_n = -0.032 * (v - V_T - 15*mvolt)/(exp(-(v - V_T - 15*mvolt)/(5*mvolt))-1) /(mvolt*ms)  : hertz
    b_n = 0.5 * exp(-(v - V_T - 10 * mvolt)/(40*mvolt))/ms                                    : hertz
    p_inf = 1/(1 + exp(-(v + 35*mvolt)/(10*mvolt)))                                           : 1
    tau_m = 1/(a_m + b_m)                                                                     : second
    tau_h = 1/(a_h + b_h)                                                                     : second
    tau_n = 1/(a_n + b_n)                                                                     : second
    tau_p = tau_p_max/(3.3 * exp((v + 35 *mvolt)/(20 * mvolt)) + exp(- (v + 35 * mvolt)/(20 * mvolt)))      : second
    mask_noise                                      : 1
    I_max                                          : amp/meter**2
    I_inj= I_max * TimedArray(t)                   : amp/meter**2
    dm/dt= a_m * (1 - m) - b_m * m                  : 1
    dh/dt= a_h * (1 - h) - b_h * h                  : 1
    dn/dt=a_n * (1 - n) - b_n * n                   : 1
    dp/dt = (p_inf - p)/tau_p                       : 1 
    ''')
#  mask                                            : 1
# I_inj_mod                                       : 1
# I_inj= I_max * TimedArray(t)                   : amp/meter**2
# I_inj = I_ramp/T_ramp * (t - 0.2*second) * TimedArray(t)                      : amp/meter**2
# dI_noise/dt = mask_noise * sigma_noise * sqrt(c_m/g_noise_L) * xi/second           : amp/meter**2

    if model == 'hh-neuron':
        eqs += Equations('''
            g_Na_mod                                        : 1
            g_Kd_mod                                        : 1
            g_M_mod                                         : 1
            g_L_mod                                         : 1
            I_Na = g_Na_mod * g_Na*m**3*h*(v-E_Na)                       : amp/meter**2
            I_Kd = g_Kd_mod * g_Kd*n**4*(v-E_K)                               : amp/meter**2
            I_M = g_M_mod * g_M*p*(v-E_K)                               : amp/meter**2
            I_L = g_L_mod * g_L * (v - E_L)                       : amp/meter**2
            dv/dt=(I_inj - I_Na - I_Kd - I_M - I_L)/c_m                 : volt                
        ''')

    elif model == 'hh-neuron-synapse':
        eqs += Equations('''
            g_Na_mod                                        : 1
            g_Kd_mod                                        : 1
            g_M_mod                                         : 1
            g_L_mod                                         : 1
            I_Na = g_Na_mod * g_Na*m**3*h*(v-E_Na)                       : amp/meter**2
            I_Kd = g_Kd_mod * g_Kd*n**4*(v-E_K)                               : amp/meter**2
            I_M = g_M_mod * g_M*p*(v-E_K)                               : amp/meter**2
            I_L = g_L_mod * g_L * (v - E_L)                       : amp/meter**2
            dv/dt=(I_inj - I_Na - I_Kd - I_M - I_L - I_syn)/c_m                 : volt
            g_syn_e                               : siemens/meter**2 
            g_syn_i                               : siemens/meter**2 
            E_AMPA = V_T * log((C0_Na_E + 1.2*C0_K_E)/(C0_Na_N + 1.2 * C0_K_N)) : volt
            E_GABA = -150 *mvolt                            : volt
            I_syn = g_syn_e * (v - E_AMPA) + g_syn_i * (v - E_GABA)  : amp/meter**2                 
        ''')

# - I_L
# + I_noise

    elif model== 'hh-ecs':
        eqs += Equations('''
            dn_Na_N/dt = -S_N/F*(I_Na + I_Na_L + 3*I_NKP) : mol
            dn_Cl_N/dt = S_N/F*(I_Cl_L + I_KCC) : mol
            dn_K_N/dt = -S_N/F*(I_Kd + I_M + I_K_L - 2*I_NKP - I_KCC) : mol
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
            g_Kd_mod                                        : 1
            g_M_mod                                        : 1
            I_Na = g_Na_mod * g_Na*m**3*h*(v-E_Na)                       : amp/meter**2
            I_Kd = g_Kd_mod * g_Kd*n**4*(v-E_K)                               : amp/meter**2
            I_M = g_M_mod * g_M*p*(v-E_K)                               : amp/meter**2
            I_Cl_L = g_Cl_L*(v-E_Cl)                              : amp/meter**2 
            I_Na_L = g_Na_L * (v-E_Na)                          : amp/meter**2 
            I_K_L = g_K_L * (v-E_K)                         : amp/meter**2 
            I_L = I_Na_L + I_Cl_L + I_K_L                     : amp/meter**2 
            sigma = 1/7 * expm1(C_Na_E/(67.3 * mM))              : 1
            f_NaK = 1/(1 + 0.1245 * exp(-0.1 * v/V_T) + 0.0365 * sigma * exp(-v/V_T)) : 1
            I_NKP = I_NKP_max * f_NaK * Hill(C_Na_N, zeta_Na, 1.5) * Hill(C_K_E, zeta_K, 1) : amp/meter**2
            I_KCC = g_KCC*(E_K - E_Cl)                          : amp/meter**2           
            dv/dt=(I_inj - I_Na - I_Na_L - I_Kd - I_M - I_K_L - I_Cl_L - I_NKP )/c_m       : volt
            Injected_current = I_max                           :amp/meter**2
            Check_1 = I_Na + I_Na_L + 3 * I_NKP                          :amp/meter**2
            Check_2 = I_Kd + I_K_L - 2 * I_NKP - I_KCC                   :amp/meter**2
            Check_3 = I_Cl_L + I_KCC                              :amp/meter**2     
        ''')
 
    ## CURRENTLY WITHOUT I_NOISE


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
            I_Kd = g_K*n*(v-E_K)                               : amp/meter**2
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

    else:
        pass

    return eqs
