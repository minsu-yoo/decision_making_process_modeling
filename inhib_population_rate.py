rate_monitor_512 =  results_512['rate_monitor_inhib']
window_width = 100.1 * b2.ms
(ax_rate) = plt.plot(figsize=(10,4))
t_max  = 500
t_min = 0 
ts = rate_monitor.t / b2.ms
idx_rate = (ts >= t_min) & (ts <= t_max)
smoothed_rates = rate_monitor_512.smooth_rate(window="flat", width=window_width)/b2.Hz

plt.plot(ts[idx_rate], smoothed_rates[idx_rate])

