import numpy
import brian2 as b2


def get_decision_time(rateA, rateB, avg_window_width=120.1*b2.ms,  rate_threshold=45.6*b2.Hz):
    
    
    # Remove which group wins 


    #def get_decision_time(rateA, rateB, avg_window_width=120.1*b2.ms,  rate_threshold=45.6*b2.Hz):
        
        
        #find which pop wins at last
    #rate_threshold= 40*b2.Hz
    #rateA = results["rate_monitor_A"]
    #rateB = results["rate_monitor_B"]

    #avg_window_width = 100.1 * b2.ms
    #(ax_rate) = plt.plot(figsize=(5,5))
    #t_max  = 1500
    #t_min = 0 
    ts = rateA.t / b2.ms
    #idx_rate = (ts >= t_min) & (ts <= t_max)
    smoothed_rateA = rateA.smooth_rate(window="flat", width=avg_window_width)/b2.Hz
    smoothed_rateB = rateB.smooth_rate(window="flat", width=avg_window_width)/b2.Hz


    # monitor each firing rate at every moment 
    decision_time_A = 0*b2.ms

    decision_time_B = 0*b2.ms

    for i in range(len(ts)):
        smoothed_rateA[i]



    # find the time when the divergen starts happening
    threshold = rate_threshold/b2.Hz

    above_thre_A = (smoothed_rateA > threshold)
    idx_over_threshold_A = numpy.where(above_thre_A==True)


    above_thre_B = (smoothed_rateB > threshold)
    idx_over_threshold_B = numpy.where(above_thre_B==True)



    if len(idx_over_threshold_A[0]) > 0:
        decision_time_A = idx_over_threshold_A[0][0];
        decision_time_A = (decision_time_A/10)*b2.ms
        
    else: 
        pass

    if len(idx_over_threshold_B[0]) > 0:
        decision_time_B = idx_over_threshold_B[0][0];
        decision_time_B = (decision_time_B/10)*b2.ms

    else: 
        pass


    print(decision_time_A, decision_time_B)


