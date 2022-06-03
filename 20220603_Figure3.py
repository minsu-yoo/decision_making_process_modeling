



import brian2 as b2
from brian2 import NeuronGroup, Synapses, PoissonInput, PoissonGroup, network_operation
from brian2.monitors import StateMonitor, SpikeMonitor, PopulationRateMonitor
from random import sample
import numpy.random as rnd
from neurodynex3.tools import plot_tools
import numpy
import matplotlib.pyplot as plt
from math import floor
import time

b2.defaultclock.dt = 0.10 * b2.ms

%matplotlib inline
from neurodynex3.competing_populations import decision_making

"""
A simple example to get started.
Returns:

"""
stim_start = 0. * b2.ms
stim_duration = 1000 * b2.ms
print("stimulus start: {}, stimulus end: {}".format(stim_start, stim_start+stim_duration))

results = decision_making.sim_decision_making_network(N_Excit=341, N_Inhib=85, weight_scaling_factor=6.0,
                                      t_stimulus_start=stim_start, t_stimulus_duration=stim_duration,
                                      coherence_level=+0.0, w_pos=2.0, mu0_mean_stimulus_Hz=250 * b2.Hz,
                                      max_sim_time=2000. * b2.ms)
plot_tools.plot_network_activity(results["rate_monitor_A"], results["spike_monitor_A"],
                                 results["voltage_monitor_A"], t_min=0. * b2.ms, avg_window_width=20. * b2.ms,
                                 sup_title="Left")
plot_tools.plot_network_activity(results["rate_monitor_B"], results["spike_monitor_B"],
                                 results["voltage_monitor_B"], t_min=0. * b2.ms, avg_window_width=20. * b2.ms,
                                 sup_title="Right")

plt.show()



def get_spike_train_ts_indices(spike_train):
    """
    Helper. Extracts the spikes within the time window from the spike train
    """
    
    t_min = 0
    t_max = 2000
    
    ts = spike_train/b2.ms
    # spike_within_time_window = (ts >= t_min) & (ts <= t_max)
    # idx_spikes = numpy.where(spike_within_time_window)
    idx_spikes = (ts >= t_min) & (ts <= t_max)
    ts_spikes = ts[idx_spikes]
    return idx_spikes, ts_spikes


    

fig, (ax_A, ax_B) = plt.subplots(2, 1)


neuron_counter = 0
spike_train_idx_list =numpy.arange(0,85)

for neuron_index in spike_train_idx_list:
    idx_spikes, ts_spikes = get_spike_train_ts_indices(all_spike_trains[neuron_index])
    ax_A.scatter(ts_spikes, neuron_counter * numpy.ones(ts_spikes.shape),
                      marker=".", c="b", s=15, lw=0)
    neuron_counter += 1
ax_A.set_ylim([0, neuron_counter])

#ax.spines['right'].set_visible(False)
ax_A.spines['top'].set_visible(False)        # 왼쪽 축을 가운데 위치로 이동
ax_A.spines['right'].set_visible(False)    
ax_A.spines['left'].set_visible(False)    
ax_A.spines['bottom'].set_visible(False)    


ax_A.get_xaxis().set_visible(False)
ax_A.get_yaxis().set_visible(False)


all_spike_trains_B = spikemon_B.spike_trains()
neuron_counter = 0
for neuron_index in spike_train_idx_list:
    idx_spikes, ts_spikes = get_spike_train_ts_indices(all_spike_trains_B[neuron_index])
    ax_B.scatter(ts_spikes, neuron_counter * numpy.ones(ts_spikes.shape),
                      marker=".", c="r", s=15, lw=0)
    neuron_counter += 1
ax_B.set_ylim([0, neuron_counter])

#ax.spines['right'].set_visible(False)
ax_B.spines['top'].set_visible(False)        # 왼쪽 축을 가운데 위치로 이동
ax_B.spines['right'].set_visible(False)    
ax_B.spines['left'].set_visible(False)    
ax_B.spines['bottom'].set_visible(False)    


ax_B.get_xaxis().set_visible(False)
ax_B.get_yaxis().set_visible(False)

