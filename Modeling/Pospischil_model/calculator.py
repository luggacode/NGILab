from brian2 import *
import math
import parameters
from scipy.optimize import brentq
from scipy.optimize import root


class Calculator:
    def __init__(self, model):
        self.model = model  
        self.params = parameters.return_initial_parameters(model)
        self.F = 9.648534297750000e+04
        self.V_T = parameters.return_initial_parameters(model)['V_T']


    #-----------------------------------------------------------------------------------------------------------------------
    ## Basic functions 
    #-----------------------------------------------------------------------------------------------------------------------
    def ThermalVoltage(self, T, F):
        return 8.314472000000000 *(T+273.15)/F
    def NernstPotential(self, x_e,x_i,z, params):
        return params['V_T']/z*np.log(x_e/x_i)
    def Hill(self,x,K,n):
        return x**n/(x**n+K**n)

    #-----------------------------------------------------------------------------------------------------------------------
    ## Gating variables
    #-----------------------------------------------------------------------------------------------------------------------
    def alpha_m(self, V, V_T):
        return -0.32 * (V - V_T - 13 *mvolt)/(math.exp(-(V - V_T - 13 * mvolt)/4/mvolt) - 1)/mvolt
    def beta_m(self, V, V_T):
        return 0.28 * (V - V_T - 40 * mvolt)/(math.exp((V - V_T - 40 * mvolt)/5/mvolt) - 1)/mvolt
    def alpha_h(self, V, V_T):
        return 0.128 * math.exp((-V - V_T - 17 * mvolt)/18/mvolt)
    def beta_h(self, V, V_T):
        return 4/(1 + math.exp(-(V - V_T - 40 * mvolt)/5/mvolt))
    def alpha_n(self, V, V_T):
        return -0.032 * (V - V_T - 15 * mvolt)/(math.exp((-(V - V_T - 15 * mvolt)/5/mvolt))- 1)/mvolt
    def beta_n(self, V, V_T):
        return 0.5 * math.exp(-(V -V_T - 10 * mvolt)/40/mvolt)
    def gating_variable_m(self, V, V_T):
        return self.alpha_m(V, V_T)/(self.alpha_m(V, V_T) + self.beta_m(V, V_T))
    def gating_variable_h(self, V, V_T):
        return self.alpha_h(V, V_T)/(self.alpha_h(V, V_T) + self.beta_h(V, V_T))
    def gating_variable_n(self, V, V_T):
        return self.alpha_n(V, V_T)/(self.alpha_n(V, V_T) + self.beta_n(V, V_T))
    def gating_variable_p(self, V):
        return 1/(1 + math.exp(- (V + 35 * mvolt)/10/mvolt))

    #-----------------------------------------------------------------------------------------------------------------------
    ## Ionic Current calculations
    #-----------------------------------------------------------------------------------------------------------------------
    def I_Na_inf(self, leakage_conductance):
        V_T = self.V_T
        return (self.params['g_Na'] * self.gating_variable_m(self.params['resting_potential'] * mvolt, V_T)**3 * self.gating_variable_h(self.params['resting_potential'] * mvolt, V_T) + leakage_conductance) * (self.params['resting_potential'] * mvolt - self.params['E_Na'])

    def I_Kd(self):
        V_T = self.V_T
        return self.params['g_Kd']*self.gating_variable_n(self.params['resting_potential'] * mvolt, V_T)**4*(self.params['resting_potential'] * mvolt-self.params['E_K'])

    def I_M(self):
        return self.params['g_M']*self.gating_variable_p(self.params['resting_potential'] * mvolt)*(self.params['resting_potential'] * mvolt-self.params['E_K'])

    def I_K_inf(self, leakage_coductance):
        return self.I_Kd() + self.I_M()

    def I_Cl_inf(self):
        return self.params['g_Cl_L'] * (self.params['resting_potential'] * mvolt - self.params['E_Cl'])

    def I_NKP(self, I_NKP_max, f_NaK):
        return I_NKP_max * f_NaK  * self.Hill(self.params['C0_Na_N'], self.params['zeta_Na'], 1.5) * self.Hill(self.params['C0_K_E'], self.params['zeta_K'], 1)

    def I_KCC(self):
        return self.params['g_KCC']*(self.params['E_K'] - self.params['E_Cl'])

    #-----------------------------------------------------------------------------------------------------------------------
    ## Setup equilibrium conditions - parameter calculations 
    #-----------------------------------------------------------------------------------------------------------------------

    ## Consistency equation for the calculation of the necessary leakage conductance of the sodium-channel
    def equilibrium_current(self, leakage_conductance):
        return 2/3 * self.I_Na_inf(leakage_conductance) + self.I_Kd() + self.I_M() + self.I_Cl_inf()

    def calc_leakage_conductance(self):
        f_single_variable = lambda x: self.equilibrium_current(x * usiemens/cm**2)
        solution = brentq(f_single_variable, 0, 100000) * usiemens/cm**2
        return solution

    def sigma(self):
        return 1/7 * (exp(self.params['C0_Na_E']/67.3/mmolar)-1)
    def f_NaK(self):  
        return 1/(1 + 0.1245 * exp(-0.1 * self.params['resting_potential'] * mvolt/self.V_T) + 0.0365 * self.sigma() * exp(-self.params['resting_potential'] * mvolt/self.V_T))

    ## Calculations I_NKP_max - BASED ON THE RELATION BETWEEN I_Na_inf and I_NKP <- I_Na_inf + 3*I_NKP = 0
    def calculate_I_NKP_max(self, I_Na_inf, f_NaK, Hill_Na, Hill_K):
        return -I_Na_inf/(3 * f_NaK * Hill_Na * Hill_K)


#-----------------------------------------------------------------------------------------------------------------------
## Value calculations
#-----------------------------------------------------------------------------------------------------------------------

# calc = Calculator('hh-ecs_exc')
# I_Kd = calc.I_Kd()
# I_M = calc.I_M()
# I_Cl = calc.I_Cl_inf()
# n = calc.gating_variable_n(calc.params['resting_potential'] * mvolt, calc.V_T) 
# p = calc.gating_variable_p(calc.params['resting_potential'] * mvolt)
# print(f'potassium delayed rectifier {I_Kd}')
# print(f'potassium muscarinic current {I_M}')
# print(f'chloride leakage: {I_Cl}')
# print(f'gating variable n: {n}')
# print(f'gating variable p: {p}')
# print(calc.params['E_K'])
# print(calc.params['g_Kd'])
# print(calc.params['g_M'])

# print(f"Thermal Voltage: {V_T}")

# ## Calculation leakage conductance 
# leakage_conductance = calc_leakage_conductance()
# ## calculation mafgnitude I_Na_inf
# I_Na_inf_1 = I_Na_inf(leakage_conductance, params['resting_potential'] * mvolt)

# ## Calculation Hill-functions, fNaK and I_NKP_max
# Hill_Na = Hill(params['C0_Na_N'], params['zeta_Na'], 1.5)
# Hill_K = Hill(params['C0_K_E'], params['zeta_K'], 1)
# f_NaK_1 = f_NaK(params['resting_potential'] * mvolt, V_T)
# I_NKP_max = calculate_I_NKP_max(I_Na_inf_1, f_NaK_1, Hill_Na, Hill_K)
# sigma_1 = sigma()

# ## Additional calc
# n = gating_variable_n(-70 * mvolt, V_T)  
# beta_n_eq = beta_n(params['resting_potential'] * mvolt, V_T)
# alpha_n_eq = alpha_n(params['resting_potential'] * mvolt, V_T)
# I_Na_eq = I_Na_inf(leakage_conductance, params['resting_potential'] * mvolt)
# I_Kd_eq = I_Kd(params['resting_potential'] * mvolt)

# print(f"leakage conductance: {leakage_conductance}")
# print(f"Sodium equilibrium current: {I_Na_inf_1}")
# print(f"sigma: {sigma_1}")
# print(f"f_NaK: {f_NaK_1}")
# print(f"I_NKP_max: {I_NKP_max}")
# print(f"beta_n: {beta_n_eq}")
# print(f"alpha_n: {alpha_n_eq}")
# print(f"gating variable n: {n}")
# print(f"I_Na_eq: {I_Na_eq}")
# print(f"I_Kd_eq: {I_Kd_eq}")
# print(f"Thermal Voltage: {V_T}")

# CURRENT_CHECK_2 = I_Kd(params['resting_potential'] * mvolt) + I_M(params['resting_potential'] * mvolt) - 2 * I_NKP(I_NKP_max, f_NaK(params['resting_potential'] * mvolt, params['V_T']), params) - I_KCC(params)
# print(f'CURRENT CHECK 2: {CURRENT_CHECK_2}')

# p = gating_variable_p(params['resting_potential'] * mvolt)
# print(f"gating_variable p: {p}")