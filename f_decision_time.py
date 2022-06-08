import numpy
import brian2 as b2

def get_decision_time(rate_monitors, avg_window_width=120.1*b2.ms,  rate_threshold=45.6*b2.Hz):
    
    
    #find which pop wins at last

    rateA = rate_monitors["rate_monitor_A"]
    rateB = rate_monitors["rate_monitor_B"]

    #window_width = 100.1 * b2.ms
    #(ax_rate) = plt.plot(figsize=(5,5))
    t_max  = 1500
    t_min = 0 
    ts = rateA.t / b2.ms
    idx_rate = (ts >= t_min) & (ts <= t_max)
    smoothed_rateA = rateA.smooth_rate(window="flat", width=avg_window_width)/b2.Hz
    smoothed_rateB = rateB.smooth_rate(window="flat", width=avg_window_width)/b2.Hz

    if (smoothed_rateB.max() < smoothed_rateA.max()):

        A_win = True
    else:

        A_win = False


    #get the firing rate differences between A and B

    if A_win:
        rate_diff = smoothed_rateA - smoothed_rateB
        winner = 'A'

     
    
    else:
        rate_diff = smoothed_rateB - smoothed_rateA
        winner = 'B'
    
    
    # find the time when the divergen starts happening
    threshold = rate_threshold/b2.Hz

    above_thre = (rate_diff > threshold)
    idx_over_threshold = numpy.where(above_thre==True)
    decision_time = idx_over_threshold[0][0];
    (decision_time/10)
    
    decision_time_in_ms = (decision_time/10)*b2.ms

    if A_win:
        decision_time_A = decision_time_in_ms
        decision_time_B = 0*b2.ms

    else:
        decision_time_A = 0*b2.ms
        decision_time_B = decision_time_in_ms


    return decision_time_A, decision_time_B



    