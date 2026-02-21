import numpy as np
from brian2 import *


def combine_mons(mon1, mon2, axis):
    return np.concatenate((np.asarray(mon1), np.asarray(mon2)), axis=axis)

def data_saver_basic(NeuronGroup , filename):
     np.savez(filename + '.npz',
                v = np.asarray(NeuronGroup.v / mvolt),     # shape: (n_neurons, n_timepoints)
                t = np.asarray(NeuronGroup.t / second),        # shape: (n_timepoints,)
                I_inj = np.asarray(NeuronGroup.I_inj / (namp/cm**2)),
                # I_max = np.asarray(I_max_mods/(namp/cm**2)),
                I_Na = np.asarray(NeuronGroup.I_Na / (namp/cm**2)),
                I_Kd = np.asarray(NeuronGroup.I_Kd / (namp/cm**2)),
                I_M = np.asarray(NeuronGroup.I_M / (namp/cm**2)),
                I_L = np.asarray(NeuronGroup.I_L / (namp/cm**2)),
                m = np.asarray(NeuronGroup.m),
                h = np.asarray(NeuronGroup.h),
                n = np.asarray(NeuronGroup.n),
                p = np.asarray(NeuronGroup.p),
                tau_m = np.asarray(NeuronGroup.tau_m / second),
                tau_h = np.asarray(NeuronGroup.tau_h / second),
                tau_n = np.asarray(NeuronGroup.tau_n / second),
                tau_p = np.asarray(NeuronGroup.tau_p / second)
                # spike_times = np.asarray(spike_mon.spike_trains()[0])
            )

def data_saver_extended(model, NeuronGroup1, NeuronGroup2, NG1spike_mon, NG2spike_mon, n_e, n_i, filename, total_connectivity = None, inhibitory_connectivity = None): 
    print(filename)
    np.savez(filename + ".npz",
            t= np.asarray(NeuronGroup1.t/second),
            #I_max = I_max,
            I_inj = combine_mons(NeuronGroup1.I_inj/(namp/cm**2), NeuronGroup2.I_inj/(namp/cm**2), 0),
            I_Na = combine_mons(NeuronGroup1.I_Na/(namp/cm**2), NeuronGroup2.I_Na/(namp/cm**2), 0),
            I_Kd = combine_mons(NeuronGroup1.I_Kd/(namp/cm**2), NeuronGroup2.I_Kd/(namp/cm**2), 0),
            I_M = combine_mons(NeuronGroup1.I_M/(namp/cm**2), NeuronGroup2.I_M/(namp/cm**2), 0),
            v = combine_mons(NeuronGroup1.v/(mvolt), NeuronGroup2.v/(mvolt), 0),
            m = combine_mons(NeuronGroup1.m, NeuronGroup2.I_inj, 0),
            h = combine_mons(NeuronGroup1.h, NeuronGroup2.I_inj, 0),
            n = combine_mons(NeuronGroup1.n, NeuronGroup2.I_inj, 0),
            p = combine_mons(NeuronGroup1.p, NeuronGroup2.I_inj, 0),
            tau_m = combine_mons(NeuronGroup1.tau_m/(second), NeuronGroup2.tau_m/(second), 0),
            tau_h = combine_mons(NeuronGroup1.tau_h/(second), NeuronGroup2.tau_h/(second), 0),
            tau_n = combine_mons(NeuronGroup1.tau_n/(second), NeuronGroup2.tau_n/(second), 0),
            tau_p = combine_mons(NeuronGroup1.tau_p/(second), NeuronGroup2.tau_p/(second), 0),
            spike_mon_exc_t = np.asarray(NG1spike_mon.t),
            spike_mon_inh_t = np.asarray(NG2spike_mon.t),
            spike_mon_exc_i = np.asarray(NG1spike_mon.i),
            spike_mon_inh_i = np.asarray(NG2spike_mon.i), 
            spike_mon_exc_count = np.asarray(NG1spike_mon.count),
            spike_mon_inh_count = np.asarray(NG2spike_mon.count),
            number_neurons = np.asarray(n_e + n_i),
            number_exc_neurons = np.asarray(n_e),
            number_inh_neurons = np.asarray(n_i),
    )
    if not total_connectivity == None:
        data = np.load(filename + '.npz')
        data_dict = dict(data)
        data_dict['total_connectivity'] = np.asarray(total_connectivity)
        np.savez(filename + '.npz', **data_dict)
    
    if not inhibitory_connectivity == None:
        data = np.load(filename + '.npz')
        data_dict = dict(data)
        data_dict['inhibitory_connectivity'] = np.asarray(inhibitory_connectivity)
        np.savez(filename + '.npz', **data_dict)
    
    if model == 'hh-neuron':
        data = np.load(filename + 'npz.')
        data_dict = dict(data)
        data_dict['I_L'] = combine_mons(NeuronGroup1.I_L, NeuronGroup2.I_L, 0)

    if model == 'hh-ecs':
        data = np.load(filename + '.npz')
        data_dict = dict(data)
        data_dict['I_NKP'] = combine_mons(NeuronGroup1.I_NKP, NeuronGroup2.I_NKP, 0)
        data_dict['I_KCC'] = combine_mons(NeuronGroup1.I_KCC, NeuronGroup2.I_KCC, 0)
        data_dict['Check_1'] = combine_mons(NeuronGroup1.Check_1, NeuronGroup2.Check_1, 0)
        data_dict['Check_2'] = combine_mons(NeuronGroup1.Check_2, NeuronGroup2.Check_2, 0)
        data_dict['Check_3'] = combine_mons(NeuronGroup1.Check_3, NeuronGroup2.Check_3, 0)
        data_dict['C_Cl_N'] = combine_mons(NeuronGroup1.C_Cl_N, NeuronGroup2.C_Cl_N, 0)
        data_dict['C_Na_N'] = combine_mons(NeuronGroup1.C_Na_N, NeuronGroup2.C_Na_N, 0)
        data_dict['C_K_N'] = combine_mons(NeuronGroup1.C_K_N, NeuronGroup2.C_K_N, 0)
        data_dict['C_Cl_E'] = combine_mons(NeuronGroup1.C_Cl_E, NeuronGroup2.C_Cl_E, 0)
        data_dict['C_Na_E'] = combine_mons(NeuronGroup1.C_Na_E, NeuronGroup2.C_Na_E, 0)
        data_dict['C_K_E'] = combine_mons(NeuronGroup1.C_K_E, NeuronGroup2.C_K_E, 0)
        np.savez(filename + '.npz', **data_dict)
    
    
    
    
def store_data(model, filename, monitors, n_e, n_i, total_connectivity = None, inhibitory_connectivity = None, **kwargs):
    if model == 'hh-neuron' or model == 'hh-ecs':
        if n_e > 0 and n_i > 0:
            exc_mon = monitors[0]
            inh_mon = monitors[1]
            spike_mon_exc = monitors[2]
            spike_mon_inh = monitors[3]
            data_saver_extended(model, exc_mon, inh_mon, spike_mon_exc, spike_mon_inh, n_e, n_i, filename)
        elif not n_i > 0:
            exc_mon = monitors[0]
            data_saver_basic(exc_mon, filename)
        elif not n_e > 0:
            inh_mon = monitors[0]
            data_saver_basic(inh_mon, filename)
    elif model == 'hh-neuron-synapse':
        exc_mon = monitors[0]
        inh_mon = monitors[1]
        exc_spike_mon = monitors[2]
        inh_spike_mon = monitors[3]
        print('es passiert')
        data_saver_extended(exc_mon, inh_mon, exc_spike_mon, inh_spike_mon, n_e, n_i, filename, total_connectivity, inhibitory_connectivity)


