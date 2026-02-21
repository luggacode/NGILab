from brian2 import *

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
            # {
            #     'variable': 'I_Na_L',
            #     'axis': 'I_Na_L in A/m^2',
            #     'plot_number' : 2,
            #     'unit': namp/cm**2
            # },
            {
                'variable': 'I_Kd', 
                'axis': 'I_Kd (nA/m^2)',
                'plot_number' : 2,
                'unit': namp/cm**2
            },
            {
                'variable': 'I_M', 
                'axis': 'I_M (nA/m^2)',
                'plot_number' : 3,
                'unit': namp/cm**2
            },
            # {
            #     'variable': 'I_Cl_L', 
            #     'axis': 'I_Cl (nA/m^2)',
            #     'plot_number' : 5,
            #     'unit': namp/cm**2 
            # },
            {
                'variable': 'I_NKP',
                'axis': 'I_NKP in nA/cm^2',
                'plot_number' : 4,
                'unit': namp/cm**2
            },
            {
                'variable': 'I_KCC',
                'axis': 'I_KCC in a/m^2',
                'plot_number' : 5,
                'unit': namp/cm**2
            },
            {
                'variable': 'v', 
                'axis': 'v (mV)',
                'plot_number' : 6,
                'unit' : mV
            },
            # {
            #     'variable': 'E_K',
            #     'axis': 'E_K (mV)',
            #     'plot_number' : 9,
            #     'unit': mV
            # },
            # {
            #     'variable': 'E_Cl',
            #     'axis': 'E_Cl (mV)',
            #     'plot_number' : 10,
            #     'unit': mV
            # },
            # {
            #     'variable': 'E_Na',
            #     'axis': 'E_Na (mV)',
            #     'plot_number' : 11,
            #     'unit': mV
            # },
            {
                'variable': 'C_Na_N',
                'axis': 'C_Na_N in mol/m^3',
                'plot_number' : 7,
                'unit': mmolar
            },
            # {
            #     'variable': 'C_Na_E',
            #     'axis': 'C_Na_E in mol/m^3',
            #     'plot_number' : 13,
            #     'unit': mmolar
            # },
            {
                'variable': 'C_K_N', 
                'axis': 'C_K_N in mol/m^3',
                'plot_number' : 8,
                'unit': mmolar
            },
            # {
            #     'variable': 'C_K_E', 
            #     'axis': 'C_K_E in mol/m^3',
            #     'plot_number' : 15,
            #     'unit': mmolar
            # },
            {
                'variable': 'C_Cl_N',
                'axis': 'C_Cl_N in mol/m^3',
                'plot_number' : 9,
                'unit': mmolar
            },
            # {
            #     'variable': 'C_Cl_E',
            #     'axis': 'C_Cl_E in mol/m^3',
            #     'plot_number' : 17,
            #     'unit': mmolar
            # },
            {
                'variable': 'Check_1',
                'axis': 'Check_1 in nA/cm^2',
                'plot_number' : 10,
                'unit': namp/cm**2
            },
            {
                'variable': 'Check_2',
                'axis': 'Check_2 in nA/cm^2',
                'plot_number' : 11,
                'unit': namp/cm**2
            },
            {
                'variable': 'Check_3',
                'axis': 'Check_3 in nA/cm^2',
                'plot_number' : 12,
                'unit': namp/cm**2
            },
            # {
            #     'variable': 'm',
            #     'axis': 'gating variable m',
            #     'plot_number' : 21,
            #     'unit': 1
            # },
            # {
            #     'variable': 'n',
            #     'axis': 'gating variable n',
            #     'plot_number' : 22,
            #     'unit': 1
            # },
            # {
            #     'variable': 'h',
            #     'axis': 'gating variable h',
            #     'plot_number' : 23,
            #     'unit': 1
            # },
            # {
            #     'variable': 'p',
            #     'axis': 'gating variable p',
            #     'plot_number' : 24,
            #     'unit': 1
            # },
            # {
            #     'variable': 'a_m',
            #     'axis': 'activation rate m a_m',
            #     'plot_number' : 25,
            #     'unit': 1
            # },
            # {
            #     'variable': 'b_m',
            #     'axis': 'deactivation rate m b_m',
            #     'plot_number' : 26,
            #     'unit': 1
            # },
            ]
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
                'variable': 'I_inj', 
                'axis': 'I_inj (nA/cm^2)',
                'plot_number' : 0,
                'unit' : namp/cm**2,
                'color': '#C81D25'
            },
            {
                'variable': 'I_Na',
                'axis': 'I_Na (nA/cm^2)',
                'plot_number' : 1,
                'unit': namp/cm**2,
                'color': '#225FCF'
            },
            {
                'variable': 'I_Kd', 
                'axis': 'I_Kd (nA/cm^2)',
                'plot_number' : 2,
                'unit': namp/cm**2,
                'color': '#225FCF'
            },
            {
                'variable': 'I_M', 
                'axis': 'I_M (nA/cm^2)',
                'plot_number' : 3,
                'unit': namp/cm**2,
                'color': '#225FCF'
            },
            {
                'variable': 'I_L', 
                'axis': 'I_L (nA/cm^2)',
                'plot_number' : 4,
                'unit': namp/cm**2,
                'color': '#225FCF'
            },
            {
                'variable': 'v', 
                'axis': 'v (mV)',
                'plot_number' : 5,
                'unit' : mvolt,
                'color': '#225FCF'
            },
            # {
            #     'variable': 'n',
            #     'axis': 'gating variable n',
            #     'plot_number' : 6,
            #     'unit': 1
            # },
            # {
            #     'variable': 'I_noise', 
            #     'axis': 'namp/cm**2',
            #     'plot_number' : 2,
            #     'unit' : namp/cm**2
            # },
            # {
            #     'variable': 'm',
            #     'axis': 'gating variable m',
            #     'plot_number' : 2,
            #     'unit': 1
            # },
            # {
            #     'variable': 'n',
            #     'axis': 'gating variable n',
            #     'plot_number' : 3,
            #     'unit': 1
            # },
            # {
            #     'variable': 'h',
            #     'axis': 'gating variable h',
            #     'plot_number' : 4,
            #     'unit': 1
            # },
            # {
            #     'variable': 'p',
            #     'axis': 'gating variable p',
            #     'plot_number' : 5,
            #     'unit': 1
            # },
            # {
            #     'variable': 'a_m',
            #     'axis': 'activation rate m a_m',
            #     'plot_number' : 6,
            #     'unit': 1
            # },
            # {
            #     'variable': 'b_m',
            #     'axis': 'deactivation rate m b_m',
            #     'plot_number' : 7,
            #     'unit': 1
            # },
            # {
            #     'variable': 'p_inf',
            #     'axis': 'p_inf',
            #     'plot_number' : 8,
            #     'unit': 1
            # }
        ]
    elif model == 'hh-neuron-addons':
        plotting_list = [
            {
                'variable': 'm',
                'axis': 'gating variable m',
                'plot_number' : 0,
                'unit': 1
            },
            {
                'variable': 'h',
                'axis': 'gating variable h',
                'plot_number' : 0,
                'unit': 1
            },
            {
                'variable': 'n',
                'axis': 'gating variable n',
                'plot_number' : 0,
                'unit': 1
            },
            {
                'variable': 'p',
                'axis': 'gating variable p',
                'plot_number' : 0,
                'unit': 1
            },
            {
                'variable': 'tau_m',
                'axis': 'time constant tau_m',
                'plot_number' : 1,
                'unit': ms
            },
            {
                'variable': 'tau_h',
                'axis': 'time constant tau_h',
                'plot_number' : 1,
                'unit': ms
            },
            {
                'variable': 'tau_n',
                'axis': 'time constant tau_n',
                'plot_number' : 1,
                'unit': ms
            },
            {
                'variable': 'tau_p',
                'axis': 'time constant tau_p',
                'plot_number' : 1,
                'unit': ms
            },

            # {
            #     'variable': 'a_m',
            #     'axis': 'activation rate m a_m',
            #     'plot_number' : 6,
            #     'unit': 1
            # },
            # {
            #     'variable': 'b_m',
            #     'axis': 'deactivation rate m b_m',
            #     'plot_number' : 7,
            #     'unit': 1
            # },
            # {
            #     'variable': 'p_inf',
            #     'axis': 'p_inf',
            #     'plot_number' : 8,
            #     'unit': 1
            # }
        ]
    return plotting_list

