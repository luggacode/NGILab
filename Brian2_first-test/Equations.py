# Equations.py
import math
import scipy.constants as spc
from sympy import symbols, solve


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

relative_permeabilties = {
    'P_K' : 10 ,
    'P_Na' : 9,
    'P_Cl' : 8
}

Concentrations_ions = {
    'C_K_N' : 10 ,
    'C_K_E' : 9,
    'C_Na_N' : 8,
    'C_Na_E' : 10 ,
    'C_Cl_N' : 9,
    'C_Cl_E' : 8
}

Temperature = 37 + 273.15

def calculate_cell_potential(gate, g):
    potential = symbols('potential')
    equation = - g['Na'] * gate['m'] * gate['h'] * (potential - 71) - g['K'] * gate['n'] * (potential - 100) - g['Cl'] * (potential - (-87))
    solution = solve(equation, potential)
    return solution



# print(calculate_resting_state_potential(relative_permeabilties, Temperature, Concentrations_ions))
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

solution = calculate_cell_potential(gating_variables, conductances)
print(solution)
print(calculate_gating_variable_m(-87.51590043))
print(calculate_gating_variable_h(-87.51590043))
print(calculate_gating_variable_n(-87.51590043))