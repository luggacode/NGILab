import numpy as np
from brian2 import *

def store_data(filename, monitors, n_e, n_i, **kwargs):
    
    data = {}
    if n_e > 0 and n_i > 0:
        monitor_names = ['sv_mon_exc','sv_mon_inh', 'spike_mon_exc', 'spike_mon_inh']
        exc_mon = monitors[0]
        inh_mon = monitors[1]
        np.savez(filename + ".npz",
         t= np.asarray(exc_mon.t/second),
         #I_max = I_max,
         I_inj = np.concatenate((np.asarray(exc_mon.I_inj/(namp/cm**2)), np.asarray(inh_mon.I_inj/(namp/cm**2))), axis=0),
         I_Na= np.concatenate((np.asarray(exc_mon.I_Na/(namp/cm**2)), np.asarray(inh_mon.I_Na/(namp/cm**2))), axis = 0), 
         I_Kd= np.concatenate((np.asarray(exc_mon.I_Kd/(namp/cm**2)), np.asarray(inh_mon.I_Kd/(namp/cm**2))), axis = 0),
         I_M = np.concatenate((np.asarray(exc_mon.I_M/(namp/cm**2)), np.asarray(inh_mon.I_M/(namp/cm**2))), axis = 0),
         I_L = np.concatenate((np.asarray(exc_mon.I_L/(namp/cm**2)), np.asarray(inh_mon.I_L/(namp/cm**2))), axis = 0),
         v = np.concatenate((np.asarray(exc_mon.v/mvolt), np.asarray(exc_mon.v/mvolt)), axis = 0),
         m = np.concatenate((np.asarray(exc_mon.m), np.asarray(inh_mon.m)), axis = 0),
         h = np.concatenate((np.asarray(exc_mon.h), np.asarray(inh_mon.h)), axis = 0),
         n = np.concatenate((np.asarray(exc_mon.n), np.asarray(inh_mon.n)), axis = 0),
         p = np.concatenate((np.asarray(exc_mon.p), np.asarray(inh_mon.p)), axis = 0),
         tau_m = np.concatenate((np.asarray(exc_mon.tau_m/second), np.asarray(inh_mon.tau_m/second)), axis = 0),
         tau_h = np.concatenate((np.asarray(exc_mon.tau_h/second), np.asarray(inh_mon.tau_h/second)), axis = 0),
         tau_n = np.concatenate((np.asarray(exc_mon.tau_n/second), np.asarray(inh_mon.tau_n/second)), axis = 0),
         tau_p = np.concatenate((np.asarray(exc_mon.tau_p/second),np.asarray(inh_mon.tau_p/second)), axis = 0),
         number_neurons = np.asarray(n_e + n_i),
         number_exc_neurons = np.asarray(n_e),
         number_inh_neurons = np.asarray(n_i)
        #  spike_times = spike_times,
        #  spike_frequency = spike_frequency,
         ## below just for synaptic connections
        #  spike_mon_exc_t = spike_mon_exc_t,
        #  spike_mon_exc_i = spike_mon_exc_i,
        #  spike_mon_inh_t = spike_mon_inh_t,
        #  spike_mon_inh_i = spike_mon_inh_i,
        #  number_neurons = number_neurons,
        #  number_exc_neurons = number_exc_neurons,
        #  number_inh_neurons = number_inh_neurons,
        #  total_connectivity = total_connectivity
         )
    elif not n_i > 0:
        monitor_names = ['sv_mon_exc', 'spike_mon_exc']
        exc_mon = monitors[0]
        for name in exc_mon.record_variables:
            arr = getattr(exc_mon, name)
            data[name] = arr
        data['t'] = getattr(exc_mon, 't')
    elif not n_e > 0:
        monitor_names = ['sv_mon_inh', 'spike_mon_inh']
        inh_mon = monitors[0]
        for name in inh_mon.record_variables:
            arr = getattr(inh_mon, name)
            data[name] = arr
        data['t'] = getattr(exc_mon, 't')


