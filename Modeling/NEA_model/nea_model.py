# Import necessary modules
import matplotlib.pyplot as plt
import scipy.constants as spc
from brian2 import *
from brian2.units.constants import faraday_constant as F
from brian2.units.constants import avogadro_constant as N_A
from params import nea_parameters, synapse_parameters
from model_eqs import get_model_eqs

# Import Warnings to avoid unnecessary warnings
import warnings as wrn
wrn.filterwarnings("ignore")
BrianLogger.suppress_name('resolution_conflict')

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
## Build User-defined convenience functions to be also called by equations in neuron models by Brian
#-----------------------------------------------------------------------------------------------------------------------
def Hill(x,K,n):
    return x**n/(x**n+K**n)
Hill = Function(Hill,arg_units=[mmolar,mmolar,1], return_unit=1,auto_vectorise=False)
Hill_cpp = '''
    #include <math.h>
    double Hill(double x,double K,double n)
    {
        return pow(x,n)/(pow(x,n)+pow(K,n));
    };
    '''
Hill.implementations.add_implementation('cpp',Hill_cpp,compiler_kwds={'headers': ['"math.h"']})

def ThermalPotential(T):
    return spc.R*(T+273.15)/spc.physical_constants['Faraday constant'][0]
ThermalPotential = Function(ThermalPotential,arg_units=[1], return_unit=1,auto_vectorise=False)
ThermalVoltage_cpp = '''
    #include <gsl/gsl_const_mksa.h>
    double ThermalPotential(const double T)
    {
        const double R = GSL_CONST_MKSA_MOLAR_GAS;
        const double F = GSL_CONST_MKSA_FARADAY;
        return R*(T+273.15)/F;
    }
    '''
ThermalPotential.implementations.add_implementation('cpp',ThermalVoltage_cpp,
                                                  compiler_kwds={'headers': ['"gsl_const_mksa.h"'],
                                                                 'include_dirs': ['/opt/anaconda3/envs/Brian2_NGILab/include/gsl/']})

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
    # V_T = ThermalVoltage(T)
    V_T = ThermalPotential(T)
    return V_T/z*np.log(x_e/x_i)
NernstPotential = Function(NernstPotential,arg_units=[mmolar,mmolar,1,1], return_unit=1,auto_vectorise=False)
NernstPotential_cpp = '''
    #include <gsl/gsl_const_mksa.h>
    const double R = GSL_CONST_MKSA_MOLAR_GAS;
    const double F = GSL_CONST_MKSA_FARADAY;
    double ThermalVoltage(const double T)
    {
        return R*(T+273.15)/F;
    };
    double NernstPotential(double x_e,double x_i,double z,const double T)
    {
        return ThermalVoltage(T)*log(x_e/x_i)/z;
    };
    '''
NernstPotential.implementations.add_implementation('cpp',NernstPotential_cpp,
                                                   dependencies={'log': DEFAULT_FUNCTIONS['log']},
                                                   compiler_kwds={'headers': ['"gsl_const_mksa.h"'],
                                                                  'include_dirs': ['/opt/anaconda3/envs/Brian2_NGILab/include/gsl/']})

def IversonBrackets(x,eps):
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
    return np.abs(x)/(np.abs(x)+eps)
IversonBrackets = Function(IversonBrackets,arg_units=[1,1], return_unit=1,auto_vectorise=False)
IversonBrackets_cpp = '''
    #include <gsl/gsl_const_mksa.h>
    double IversonBrackets(double x, const double eps)
    {
        double z, e;
        /* Handle degenerate epsilon */
        if (eps <= 0.0) {
            return (x >= 0.0) ? 1.0 : 0.0;
        }
    
        z = x / eps;
        /* Numerically stable sigmoid */
        if (z >= 0.0) {
            e = exp(-z);
            return 1.0 / (1.0 + e);
        } else {
            e = exp(z);
            return e / (1.0 + e);
        }
    };
    '''
IversonBrackets.implementations.add_implementation('cpp',IversonBrackets_cpp,
                                                   dependencies={'abs': DEFAULT_FUNCTIONS['abs']})

def DiffusionFlux(delta_c,c_theta,omega_c,D_coeff):
    """
    Simmetric diffusion flux with saturation (from Lallouette et al., Comput. Gliosci. 2019).

    Input parameters (w/out units):
    - delta_c   : float    Concentration gradient between compartment
    - c_theta   : float    Concentration threshold for diffusion
    - omega_c   : float    Scaling factor
    - D_coeff   : float    Max diffusion rate

    Return:
    - J_diff    : Diffusion Flux
    """
    return -D_coeff/2 * (1 + np.tanh((np.abs(delta_c) - c_theta)/omega_c))*np.sign(delta_c)
DiffusionFlux = Function(DiffusionFlux,arg_units=[mmolar,mmolar,mmolar,mmolar/second], return_unit=mole/second,auto_vectorise=False)
DiffusionFlux_cpp = '''
    #include <math.h>
    int sgn(double val) {
        return (0 < val) - (val < 0);
    };
    double DiffusionFlux(double delta_c, double c_theta, double omega_c, double D_coeff)
    {
        return -D_coeff/2 * (1 + tanh((abs(delta_c) - c_theta)/omega_c))*sgn(delta_c);
    };
    '''
DiffusionFlux.implementations.add_implementation('cpp',DiffusionFlux_cpp,
                                                   dependencies={'abs': DEFAULT_FUNCTIONS['abs'],
                                                                 'tanh': DEFAULT_FUNCTIONS['tanh']})

# -----------------------------------------------------------------------------------------------------------------------
# Dummy cell stimulation models
# -----------------------------------------------------------------------------------------------------------------------
def periodic_nodes(N,rates,name='period*',dt=None):
    if isinstance(rates,input.timedarray.TimedArray):
        eqs = Equations('dv/dt = stimulus(t) : 1')
        nmspace = {'stimulus': rates}
    else:
        eqs = Equations('''
                        rate : Hz
                        dv/dt = rate : 1
                        ''')
        nmspace = None
    cells = NeuronGroup(N,eqs,
        threshold='v>=1.0',
        reset='v=0.0',
        name=name,
        namespace=nmspace,
        method='euler',
        dt=dt)
    if not isinstance(rates,str): cells.rate = rates
    cells.v = 0.0
    return cells

@check_units(trp=second)
def poisson_nodes(N,rates,trp=0.0*second,name='poiss*',dt=None):
    # Create equations to also handle case of tme-variant stimulation
    if isinstance(rates,input.timedarray.TimedArray):
        eqs = Equations('rate = stimulus(t) : Hz')
        nmspace = {'stimulus': rates}
    else:
        eqs = Equations('rate : Hz')
        nmspace = None
    cells = NeuronGroup(N,eqs,
        threshold='rand()<rate*dt',
        refractory=trp,
        name=name,
        namespace=nmspace,
        dt=dt)
    if not isinstance(rates,input.timedarray.TimedArray): cells.rate = rates
    return cells

# -----------------------------------------------------------------------------------------------------------------------
# Neuron-ECS-Astrocyte (NEA) Model
# -----------------------------------------------------------------------------------------------------------------------
def nea_node(params,sinput='glu',name='nea*',dt=0.1*us):

    eqs_string = get_model_eqs(sinput)
    ## Generate Equations
    eqs = Equations(eqs_string)

    ## Generate the neuron group
    nea = NeuronGroup(1, eqs,
                      # events=events,
                      namespace=params,
                      name=name,
                      method='euler',
                      dt=dt,
                      order=10)

    # Constants
    nea.GABA0_e = params['GABA0_e']
    nea.G0_e = params['G0_e']

    # Initialize variables
    nea.n_Na_n = 0*mole
    nea.n_K_n = 0*mole
    nea.n_Cl_n = 0*mole

    nea.n_Na_a = 0*mole
    nea.n_K_a = 0*mole
    nea.n_Cl_a = 0*mole
    nea.N_a_d = params['N0_a']
    nea.K_a_d = params['K0_a']
    nea.C_a_d = params['C0_a']
    nea.r_GABA = 0.

    nea.n_Na_e = 0*mole
    nea.n_K_e = 0*mole
    nea.n_Cl_e = 0*mole
    nea.N_e_d = params['N0_e']
    nea.K_e_d = params['K0_e']
    nea.C_e_d = params['C0_e']

    if sinput=='glu':
        nea.n_Glu_e = 0*mole
        nea.n_Glu_s = 0*mole
        nea.n_GABA_e = 0*mole
    else:
        nea.n_Glu_e = 0*mole
        nea.n_GABA_s = 0*mole
        nea.n_GABA_e = 0*mole

    nea.v = -70*mV
    nea.v_a = -90*mV

    return nea

# -----------------------------------------------------------------------------------------------------------------------
# Synaptic connections
# -----------------------------------------------------------------------------------------------------------------------
def synaptic_connection(stim_source, snc_target, params, sinput='glu', name='syn*', dt=0.1*us, delay=None):
    eqs = Equations('''
        # Bound fraction of postsynaptic receptors (assuming no activation at Nt_s=0)
        Nt_s : mmolar
        J_rec = -r_clip/tau_r+J*(1-r_clip)*Nt_s : 1/second
        J_r_post = J_rec : 1/second (summed)
        dr/dt = J_rec   : 1 (clock-driven)
        r_clip = clip(r,0,1)        : 1        
        ''')
    if sinput=='glu':
        eqs += Equations('''
                         G_AMPA_post = g*r_clip : siemens/meter**2 (summed)
                         ''')
        on_pre = '''
                 Nt_s = G_s_post              
                 n_Glu_s_post += Nt_rel*Lambda_s
                 '''
    else:
        eqs += Equations('''
                         G_GABA_post = g*r_clip : siemens/meter**2 (summed)
                         ''')
        on_pre = '''
                 Nt_s = GABA_s_post              
                 n_GABA_s_post += Nt_rel*Lambda_s
                 '''

    synapse = Synapses(stim_source, snc_target, eqs,
                       on_pre=on_pre,
                       namespace=params,
                       method='euler',
                       name=name,
                       delay=delay,
                       dt=dt,
                       order=0)

    return synapse

def nea_simulator(N_synapses,duration,
                  protocol='periodic',sinput='glu',
                  code_dir='./codegen/'):
    # Clean memory from previous builds (allowing multiple runs)
    device.delete(force=True)
    dt_sim = 0.05*us
    # Make sure that the number of incoming synaptic connections is an integer
    N_synapses = int(N_synapses)

    # Generate parameters
    pars_nea = nea_parameters()
    D = pars_nea['D_Glu_e'] if sinput=='glu' else pars_nea['D_GABA_e']
    pars_syn = synapse_parameters(D=D)
    pars_syn['Nt0'] = pars_nea['G0_e'] if sinput=='glu' else pars_nea['GABA0_e']
    pars_nea['R_T'] = pars_syn['R_T']

    # Generate synaptic connections
    if protocol=='periodic':
        stim = periodic_nodes(1,10*Hz,name='period*',dt=1*ms)
        i_pre = np.atleast_1d([0]).astype(int)
        j_pst = np.atleast_1d([0]).astype(int)
    elif protocol=='poisson':
        stim = poisson_nodes(N_synapses,10*Hz,trp=0.0*second,name='poiss*',dt=None)
        i_pre = np.arange(N_synapses).astype(int)
        j_pst = np.zeros(N_synapses).astype(int)
    assert sinput in ['glu', 'gaba'], "Transmitter type (ttype) can only be of 'glu' or 'gaba'"

    # Generate NEA node
    nea = nea_node(pars_nea,sinput=sinput,name='nea*',dt=dt_sim)
    nea.Lambda_e = pars_nea['Lambda_e']
    nea.Lambda_s = pars_syn['Lambda_s']

    # Set up synaptic connection
    syn = synaptic_connection(stim, nea, pars_syn, sinput=sinput, name='syn*', dt=dt_sim)
    syn.connect(i=i_pre, j=j_pst)
    syn.r = 0
    syn.Nt_s = 0*mmolar

    # Set up useful monitors
    mon_spk = SpikeMonitor(stim, record=True, name='spk')
    vnea = ['v','v_a','I_AMPA'] if sinput=='glu' else ['v','v_a','I_GABA']
    vnea += ['K_e','K_a','K_e_d','K_a_d']
    mon_nea = StateMonitor(nea,variables=vnea,record=True,dt=0.1*ms) if sinput=='glu' else StateMonitor(nea,variables=vnea,record=True,dt=0.1*ms)

    # Build Network
    network = Network([stim, nea, syn, mon_spk, mon_nea])

    ## Run the simulator
    network.run(duration=duration*second, report='text')
    device.build(directory=code_dir, clean=True)

    _, axs = plt.subplots(4, 1, figsize=(6, 10))
    axs[0].vlines(mon_spk.t_,mon_spk.i[:].T,mon_spk.i[:].T+0.9)
    axs[0].set(xticklabels='')
    axs[0].set_ylabel('Stimulation')

    # Synaptic current
    try:
        axs[1].plot(mon_nea.t_,mon_nea.I_AMPA[:].T/(nA/cm**2),'k-')
        axs[1].set(xticklabels='')
        axs[1].set_ylabel(r'$I_{AMPA}$ (nA/$\mu$m$^2$)')
    except:
        axs[1].plot(mon_nea.t_,mon_nea.I_GABA[:].T/(nA/cm**2),'k-')
        axs[1].set(xticklabels='')
        axs[1].set_ylabel(r'$I_{GABA}$ (nA/$\mu$m$^2$)')

    # Membrane Potentials
    axs[2].plot(mon_nea.t_,mon_nea.v[:].T/mV,'k-',label='N')
    axs[2].plot(mon_nea.t_,mon_nea.v_a[:].T/mV,'g-',label='A')
    axs[2].set_ylabel('v (mV)')
    axs[2].legend(loc='upper right')

    # Potassium Concentrations
    axs[3].plot(mon_nea.t_,mon_nea.K_e[:].T/mM,'b-',label='E')
    axs[3].plot(mon_nea.t_,mon_nea.K_e_d[:].T/mM,'b--',label=r'E$_\infty$')
    axs[3].plot(mon_nea.t_,mon_nea.K_a[:].T/mM,'c-',label='A')
    axs[3].plot(mon_nea.t_,mon_nea.K_a[:].T/mM,'c--',label=r'A$_\infty$')
    axs[3].set_xlabel('Time (s)')
    axs[3].set_xlabel('K$^+$ (mM)')
    axs[3].legend(loc='upper right')

if __name__=="__main__":
    # -----------------------------------------------------------------------------------------------------------------------
    # Imports and Preamble
    # -----------------------------------------------------------------------------------------------------------------------
    import matplotlib.pyplot as plt

    code_dir = './codegen'
    prefs.GSL.directory = '/opt/anaconda3/envs/Brian2_NGILab/include/'  ## The directory where the GSL library headings are found
    set_device('cpp_standalone',directory=code_dir,build_on_run=False)
    prefs.devices.cpp_standalone.openmp_threads = 1  ## The number of threads used in the parallelization (machine-dependent)
    prefs.logging.file_log = False
    prefs.logging.delete_log_on_exit = True

    # -----------------------------------------------------------------------------------------------------------------------
    # Testing
    # -----------------------------------------------------------------------------------------------------------------------
    nea_simulator(5,1.0,protocol='poisson',sinput='glu',code_dir='./codegen/')

    # -----------------------------------------------------------------------------------------------------------------------
    # Visualize
    # -----------------------------------------------------------------------------------------------------------------------
    plt.show()