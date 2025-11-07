# Equations.py
import math
import scipy.constants as spc
import sympy as sp
from brian2 import *


def calculate_gating_variable_m(equilibrium_potential):
    x = (182*(equilibrium_potential - (-38))/(1-math.exp((-38 - equilibrium_potential)/6)))/((182*(equilibrium_potential - (-38))/(1-math.exp((-38 - equilibrium_potential)/6)) + (124 *(-38 - equilibrium_potential))/1-math.exp((equilibrium_potential - (-38))/6)))
    return x

def calculate_gating_variable_h(equilibrium_potential):
    x = (15*(-66 - equilibrium_potential)/(1 - math.exp((equilibrium_potential - (-66))/6)))/(15*(-66 - equilibrium_potential)/(1 - math.exp((equilibrium_potential - (-66))/6))+(15*(equilibrium_potential - (-66))/(1 - math.exp((-66 - equilibrium_potential)/6))))
    return x

def calculate_gating_variable_n(equilibrium_potential):
    x = 1/(1+math.exp((18.7 - equilibrium_potential)/9.7))
    return x

def calculate_resting_state_potential(P, T, C):
    resting_potential = spc.R * T/spc.physical_constants['Faraday constant'][0] * math.log((P['P_K'] * C['C_K_E'] + P['P_Na'] * C['C_Na_E'] + P['P_Cl'] * C['C_Cl_N'])/ P['P_K'] * C['C_K_N'] + P['P_Na'] * C['C_Na_N'] + P['P_Cl'] * C['C_K_E'])
    return resting_potential

def calculate_cell_potential(gate, g):
    potential = sp.symbols('potential')
    gleichung = - g['Na'] * gate['m'] * gate['h'] * (potential - 71) - g['K'] * gate['n'] * (potential - 100) - g['Cl'] * (potential - (-87))
    solution = sp.solve(gleichung, potential)
    return solution

def return_HH_equations(model = 'hh-neuron'):
    
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
    I_inj=I_max * TimedArray(t)                   : amp/meter**2
    dm/dt=(m_inf-m)/tau_m                           : 1
    dh/dt=(h_inf-h)/tau_h                           : 1
    dn/dt=(n_inf-n)/tau_n                           : 1
    ''')

    if model=='hh-ecs':
        eqs += Equations('''
            dn_Na_N/dt = -S_N/F*(I_Na + 3*I_NKP) : mol
            dn_Na_E/dt = S_N/F*(I_Na + 3*I_NKP) : mol
            dn_Cl_N/dt = S_N/F*(I_Cl + I_KCC) : mol
            dn_Cl_E/dt = -(S_N/F*(I_Cl + I_KCC)) : mol
            dn_K_N/dt = -S_N/F*(I_K - 2*I_NKP - I_KCC) : mol
            dn_K_E/dt = S_N/F*(I_K - 2*I_NKP - I_KCC) : mol
            C_Na_N = clip(C0_Na_N + n_Na_N/V_N, 0 * mM, 155 * mM) : mM
            C_Na_E = clip(C0_Na_E + n_Na_E/V_E, 0 * mM, 155 * mM) : mM
            C_Cl_N = clip(C0_Cl_N + n_Cl_N/V_N, 0 * mM, 133 * mM) : mM
            C_Cl_E = clip(C0_Cl_E + n_Cl_E/V_E, 0 * mM, 133 * mM) : mM
            C_K_N = clip(C0_K_N + n_K_N/V_N, 0 * mM, 135 * mM) : mM
            C_K_E = clip(C0_K_E + n_K_E/V_E, 0 * mM, 135 * mM)  : mM
            E_Na = NernstPotential(C_Na_E, C_Na_N, 1, T_exp)*volt : volt
            E_K = NernstPotential(C_K_E, C_K_N, 1, T_exp)*volt : volt
            E_Cl = NernstPotential(C_Cl_E, C_Cl_N, -1, T_exp)*volt : volt          
            I_Na = (g_Na*m**3*h + g_Na_L)*(v-E_Na)                       : amp/meter**2
            I_K = g_K*n*(v-E_K)                               : amp/meter**2
            I_Cl = g_Cl*(v-E_Cl)                              : amp/meter**2 
            sigma = 1/7 * (expm1(C_Na_E/(67.3 * mM))-1)              : 1
            f_NaK = 1/(1 + 0.1245 * exp(-0.1 * v/V_T) + 0.0365 * sigma * expm1(-v/V_T)) : 1
            I_NKP = I_NKP_max * f_NaK * 1/(1 + (zeta_Na/C_Na_N)**1.5)*Hill(C_K_E, zeta_K, 1) : amp/meter**2
            I_KCC = g_KCC*(E_K - E_Cl)                          : amp/meter**2           
            dv/dt=(I_inj - I_Na - I_K - I_Cl + I_NKP)/c_m       : volt
            Check_1 = I_Na + 3 * I_NKP                          :amp/meter**2
            Check_2 = I_K + 2 * I_NKP - I_KCC                   :amp/meter**2
            Check_3 = I_Cl + I_KCC                              :amp/meter**2
        ''')
    elif model == 'hh-neuron':
        eqs += Equations('''
            I_Na=g_Na*m**3*h*(v-E_Na)                       : amp/meter**2
            I_K=g_K*n*(v-E_K)                               : amp/meter**2
            I_Cl=g_Cl*(v-E_Cl)                              : amp/meter**2
            dv/dt=(I_inj-I_Na-I_K-I_Cl)/c_m                 : volt   
        ''')
    else:
        ## w/ ECS
        pass

    return eqs

#  I_NKP = I_NKP_max*(1 + 0.1245 * exp(-v/(10*V_T)) - 0.0052 * exp(-v/V_T) * (1 - exp(C_Na_E/(67.3*mM)))) * Hill(C_Na_N, zeta_Na, 1.5) * Hill(C_K_E, zeta_K, 1) : amp/meter**2
# f_NaK = 1/(1 + 0.1245 * exp(-0.1 * v/V_T) + 0.0365 * sigma * exp(-v/V_T))

gating_variables = {
    'm': 0.01,
    'n': 0.01,
    'h': 0.99
}

conductances = {
    'K': 6.93,
    'Na': 2.04,
    'Cl': 0.0000338
}
