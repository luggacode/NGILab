# Equations.py
import math
import scipy.constants as spc
from scipy.optimize import brentq
from scipy.optimize import root
import sympy as sp
from brian2 import *


## Initial variables

F = 9.648534297750000e+04


conductances = {
    'K': 0.693e6,
    'Na': 2.04e6,
    'Cl': 50,
}

concentrations = {## Concentrations to setup reverse potentials
            'C0_Na_N': 10,
            'C0_Na_E': 145,
            'C0_K_N': 130,
            'C0_K_E': 3,
            'C0_Cl_N': 5,
            'C0_Cl_E': 130
}

ion_equilibria = {
}

zetas = {
    'zeta_Na' : 13,
    'zeta_K' : 0.2
}

## Basic functions
def ThermalVoltage(T):
    return 8.314472000000000 *(T+273.15)/9.648534297750000e+04
def NernstPotential(x_e,x_i,z,T):
    V_T = ThermalVoltage(T)
    return V_T/z*np.log(x_e/x_i)
def Hill(x,K,n):
    return x**n/(x**n+K**n)

V_T = ThermalVoltage(37) * 1000

## Concentrations and ionic variations
def C_Cl_N(C0_Cl_N, n_Cl_N, Volume):
    return C0_Cl_N + n_Cl_N/Volume
def C_Na_N(C0_Na_N, n_Na_N, Volume):
    return C0_Na_N + n_Na_N/Volume
def C_Na_N(C0_K_N, n_K_N, Volume):
    return C0_K_N + n_K_N/Volume
def dn_Na_N(I_Na, I_NKP, surface):
    return -surface/F*(I_Na + 3*I_NKP)
def dn_K_N(I_K, I_NKP, I_KCC, surface):
    return -surface/F*(I_K - 2*I_NKP - I_KCC)    
def dn_Cl_N(I_Cl, I_KCC, surface):
    return surface/F*(I_Cl + I_KCC)

## Calculation reversal potential for each ion
ion_equilibria['E_K'] = NernstPotential(concentrations['C0_K_E'],concentrations['C0_K_N'],1, 37) * 1000
ion_equilibria['E_Na'] = NernstPotential(concentrations['C0_Na_E'],concentrations['C0_Na_N'],1, 37) * 1000
ion_equilibria['E_Cl'] = NernstPotential(concentrations['C0_Cl_E'],concentrations['C0_Cl_N'],-1, 37) * 1000


## Calculation gating variables 
def alpha_m(V):
    return 182*(V - (-38))/(1 - math.exp((-38 - V)/6))
def beta_m(V):
    return 124*(-38 - V)/(1 - math.exp((V - (-38))/6))
def alpha_h(V):
    return 15*(-66 - V)/(1 - math.exp((V - (-66))/6))
def beta_h(V):
    return 15*(V - (-66))/(1 - math.exp((-66 - V)/6))
def gating_variable_m(V):
    return alpha_m(V)/(alpha_m(V) + beta_m(V))
def gating_variable_h(V):
    return alpha_h(V)/(alpha_h(V) + beta_h(V))
def gating_variable_n(V):
    return 1/(1+math.exp((18.7 - V)/9.7))   

## Calculations ionic currents 1 - VARIABLE VOLTAGE
def I_Cl(V, Reversal_potential): 
    return conductances['Cl'] * (V - Reversal_potential)
def I_Na(V, Reversal_potential):   
    return conductances['Na'] * (gating_variable_m(V)**3) * gating_variable_h(V) * (V - Reversal_potential)
def I_K(V, Reversal_potential):    
    return conductances['K'] * (gating_variable_n(V)) * (V - Reversal_potential)
def sigma():
    return 1/7 * (exp(concentrations['C0_Na_E']/67.3)-1)
# TODO: make sigma generic by handing the concentration to the sigma function
def f_NaK(V, V_T):  
    return 1/(1 + 0.1245 * exp(-0.1 * V/V_T) + 0.0365 * sigma() * exp(-V/V_T))

def I_NKP(I_NKP_max, f_NaK, concentrations, zetas):
    V_T = ThermalVoltage(37) * 1000
    return I_NKP_max * f_NaK  * Hill(concentrations['C0_Na_N'], zetas['zeta_Na'], 1.5) * Hill(concentrations['C0_K_E'], zetas['zeta_K'], 1)
def I_KCC(conductances, ion_equilibria):
    return conductances['g_KCC']*(ion_equilibria['E_K'] - ion_equilibria['E_Cl'])



## Calculations ionic currents_inf 2 - AT RESTING POTENTIAL (ik that they aren't necessary with the functions above already being defined but hey)
def I_Na_inf(leakage_conductance):
    return (conductances['Na'] * gating_variable_m(-70)**3 * gating_variable_h(-70) + leakage_conductance) * (-70 - ion_equilibria['E_Na'])
def I_Na_inf_ohne_lc():
    return (conductances['Na'] * gating_variable_m(-70)**3 * gating_variable_h(-70)) * (-70 - ion_equilibria['E_Na'])
def I_K_inf():
    return conductances['K'] * gating_variable_n(-70) * (-70 - (ion_equilibria['E_K']))
# TODO: include leakage conductance potassium channel
def I_Cl_inf():
    return conductances['Cl'] * (-70 - ion_equilibria['E_Cl'])

## Calculations I_NKP_max - BASED ON THE RELATION BETWEEN I_Na_inf and I_NKP <- I_Na_inf + 3*I_NKP = 0
def calculate_I_NKP_max(I_Na_inf, f_NaK, Hill_Na, Hill_K):
    return -I_Na_inf/(3 * f_NaK * Hill_Na * Hill_K)

## Consistency equation for the calculation of the necessary leakage conductance of the sodium-channel
def equilibrium_current(leakage_conductance):
    return 2/3 * I_Na_inf(leakage_conductance) + I_K_inf() + I_Cl_inf()


## Calculation of the total currents in HH and HH-ECS
def total_current_hh(V):
    return   - I_K_inf() - I_Na_inf(leakage_conductance) - I_Cl(V, ion_equilibria['E_Cl'])
# TODO: include dynamic potential for the I_Na_inf
def total_current_hhecs(leakage_conductance, x):
    return   - I_K_inf() - I_Na_inf(leakage_conductance) - I_Cl(-70, ion_equilibria['E_Cl']) - I_NKP(-70, f_NaK(-70, V_T, sigma()), concentrations, zetas)

## Calculation of the resting state Potential in HH + Calculation conductance of KCC + Calculation leakage conductance of the sodium-channel
def calculate_resting_state_potential_hh():
    resting_potential = brentq(total_current_hh, -100, 50)
    return resting_potential
def conductance_KCC(V, Nernst_Cl, Nernst_K, conductance_Cl):
    return conductance_Cl * (V - Nernst_Cl)/(Nernst_Cl - Nernst_K)
def calc_leakage_conductance():
    leakage_conductance = brentq(equilibrium_current, 0, 50000)
    return leakage_conductance


conductances['g_KCC'] = conductance_KCC(-70, ion_equilibria['E_Cl'], ion_equilibria['E_K'], conductances['Cl'])

## ACTUAL VALUE CALCULATIONS

g_KCC = conductances['g_KCC']
print(f'CONDUCTANCE KCC: {g_KCC}')

# 2 variables missing to achieve equilibrium for the HH-ECS: 
# 1. additional leakage conductance of the sodium-channel
# 2. I_NKP_max for the sodium-potassium pump (derived from the I_Na_inf in this script but also possible from I_K_inf and I_KCC_inf - see 2.27)

## Calculations of VALUE 1 (g_Na)
leakage_conductance = calc_leakage_conductance()
print(f'SODIUM LEAKAGE CONDUCTANCE: {leakage_conductance}')

## Calculation of VALUE 2 (I_NKP_max) 
I_Na_inf_calc = I_Na_inf(leakage_conductance)
f_NaK_calc = f_NaK(-70, V_T)
Hill_Na = Hill(concentrations['C0_Na_N'], zetas['zeta_Na'], 1.5)
Hill_K = Hill(concentrations['C0_K_E'], zetas['zeta_K'], 1)

I_NKP_max = calculate_I_NKP_max(I_Na_inf_calc, f_NaK_calc, Hill_Na, Hill_K)

print(f'I_NKP_max: {I_NKP_max}')



#### CHECKS for error fixing

## CURRENT-CHECK_1 -> expected to be 0 as I_Na_inf = - 3 * I_NKP
CURRENT_CHECK_1 = I_Na_inf(leakage_conductance) + 3 * I_NKP(I_NKP_max, f_NaK(-70, V_T), concentrations, zetas)
print(f'CURRENT CHECK 1: {CURRENT_CHECK_1}')

## CURRENT-CHECK_2 -> expected to be 0 as I_K_inf = 2 * I_NKP + I_KCC
CURRENT_CHECK_2 = I_K_inf() - 2 * I_NKP(I_NKP_max, f_NaK(-70, V_T), concentrations, zetas) - I_KCC(conductances, ion_equilibria)
print(f'CURRENT CHECK 2: {CURRENT_CHECK_2}')

## CURRENT-CHECK_3 -> expected to be 0 as I_K_inf = 2 * I_NKP + I_KCC
CURRENT_CHECK_3 = I_Cl_inf() + I_KCC(conductances, ion_equilibria)
print(f'CURRENT CHECK 3: {CURRENT_CHECK_3}')

## POTENIAL_CHECK -> expected to show that the calculated conductance matches the potential
POTENTIAL_CHECK = -I_KCC(conductances, ion_equilibria)/conductances['Cl'] + ion_equilibria['E_Cl']
print(f'POTENTIAL CHECK: {POTENTIAL_CHECK}')


## TOTAL CURRENT CHECK
## Check whether the total ionic current is 0 and EQUILIBRIUM is achieved
total_current = I_Na_inf(leakage_conductance) + I_K_inf() + I_Cl_inf() + I_NKP(I_NKP_max, f_NaK(-70, V_T), concentrations, zetas)
print(f'TOTAL CURRENT: {total_current}')

print('NERNST Cl')
print(NernstPotential(concentrations['C0_Cl_E'], concentrations['C0_Cl_N'], -1, 37))
print('NERNST K')
print(NernstPotential(concentrations['C0_K_E'], concentrations['C0_K_N'], 1, 37))
print('Na_Nernst')
print(NernstPotential(concentrations['C0_Na_E'], concentrations['C0_Na_N'], 1, 37))

print(f'I_Na_inf {I_Na_inf_calc}')

I_K_CURRENT = I_K_inf()
print(f'I_K: {I_K_CURRENT}')
print('n')
print(gating_variable_n(-70))
print(ion_equilibria['E_K'])

I_Cl_CURRENT = I_Cl_inf()
print(f'I_Cl_CURRENT: {I_Cl_CURRENT}')

I_KCC_CURRENT = I_KCC(conductances, ion_equilibria)
print(f'I_KCC_CURRENT: {I_KCC_CURRENT}')

NKP_current = I_NKP(I_NKP_max, f_NaK(-70, V_T), concentrations, zetas)
print(f'NKP_CURRENT: {NKP_current}')


# resting_potential = calculate_resting_state_potential_hh()
# print(f'RESTING POTENTIAL: {resting_potential}')






# def calculate_resting_state_potential_hhecs():
#     amplitude_current = brentq(total_current_hhecs, -100, 50)
#     return amplitude_current
# amplitude = calculate_resting_state_potential_hhecs()
# print(f"Equilibrium (resting) potential: {amplitude:.8f} ")
# Variance_K = dn_K_N(I_K(-89.89025359, ion_equilibria['E_K']), I_NKP(-89.89025359, 37, concentrations, zetas), I_KCC(conductances, ion_equilibria), 1)
# Variance_Na = dn_Na_N(I_Na(-89.89025359, ion_equilibria['E_Na']), I_NKP(-89.89025359, 37, concentrations, zetas), 1)
# Variance_Cl = dn_Cl_N(I_Cl(-89.89025359, ion_equilibria['E_Cl']), I_KCC(conductances, ion_equilibria), 1)
# current = total_current_hhecs()

