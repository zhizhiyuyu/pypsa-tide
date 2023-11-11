from matplotlib import rc
from numpy import array
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from ..modules.parameterisations import hill_chart_parametrisation_h
import os

# line_split and line_area_points are for thetis simulations.

def line_split(start, end, segments):
    x_delta = (end[0] - start[0]) / float(segments)
    y_delta = (end[1] - start[1]) / float(segments)
    points = []
    for i in range(1, segments):
        points.append([start[0] + i * x_delta, start[1] + i * y_delta])
    return array([start] + points + [end])


def plot_lagoon_output(data, file = "outputs/output.png", options = {"t_s": 0, "plot_t" : 180}):
    index = {"t": 0, "h_o": 1, "h_i": 2, "DZ": 3, "P": 4, "E": 5, "m": 6, "Q_t": 7, "Q_s": 8,
             "m_dt": 9, "m_t": 10, "f_r": 11}

    f, axarr = plt.subplots(3, sharex = "all", sharey = "none", figsize=(6,4), dpi =200)

    axarr[0].plot(data[:,index["t"]]/3600 - options["t_s"], data[:,index["h_i"]],label = " $\eta_{i}$", linewidth = 0.45,color='red')
    axarr[0].plot(data[:,index["t"]]/3600 - options["t_s"], data[:,index["h_o"]],label = " $\eta_{o}$",  linewidth = 0.25,color='blue',ls='dashed')
    axarr[0].legend(loc = 2, frameon=False, fontsize='x-small', labelspacing = 0.0)
    axarr[0].set_ylabel(r"$\eta (m)$ ", labelpad = 10 , fontsize = 8, weight = "bold")
    axarr[0].set_ylim([-7, 7])

    axarr[1].plot(data[:,index["t"]]/3600 - options["t_s"], data[:,index["Q_t"]]/1000,label = "$Q_{turb}$ ", linewidth = 0.45,color='red')
    axarr[1].plot(data[:,index["t"]]/3600 - options["t_s"], data[:,index["Q_s"]]/1000,label = "$Q_{sl}$ ", linewidth = 0.45,color='black',ls='dashed')
    axarr[1].legend(loc = 2, frameon=False, fontsize='x-small', labelspacing = 0.0)
    axarr[1].set_ylabel(r"$Q (m^3/s) \times 10^3$ ", labelpad = 5, fontsize = 8, weight = "bold")
    axarr[1].set_ylim([-250, 300])

    axarr[2].plot(data[:,index["t"]]/3600 - options["t_s"], data[:,index["P"]]/1000,label = "Power (GW)", linewidth = 0.45,color='black')  
    #axarr[2].legend(loc = 1, frameon=False, fontsize='x-small', labelspacing = 0.0)
    axarr[2].set_ylabel(r"$Power (GW)$ ", labelpad = 2, fontsize = 8, weight = "bold")
    axarr[2].set_ylim([0, 13])
    
    ax2 = axarr[2].twinx()  
    ax2.plot(data[:,index["t"]]/3600 - options["t_s"], data[:,index["E"]]/1000,label = "Electricity (GWh)", linewidth = 0.45,color='red', ls='--')  
    ax2.set_ylabel(r"$Electricity (GWh)$ ", labelpad = 2, fontsize = 8, weight = "bold")
    ax2.set_ylim([0, 420])
    #ax2.legend(loc=1, fontsize='x-small')

    # add minor ticks
 
    axarr[0].minorticks_on()
    plt.setp(axarr[0].get_yticklabels(), fontsize=8)
    axarr[1].minorticks_on()
    plt.setp(axarr[1].get_yticklabels(), fontsize=8)
    axarr[2].minorticks_on()
    plt.setp(axarr[2].get_yticklabels(), fontsize=8)
    ax2.minorticks_on()
    plt.setp(ax2.get_yticklabels(), fontsize=8)  
    # xticks
    plt.setp(axarr[2].get_xticklabels(), fontsize=8)
    
    plt.xlabel("Time $t$ (hours)", weight = "bold", fontsize = 8)
    plt.xlim(0, options["plot_t"])
    
    f.subplots_adjust(hspace=0)
    plt.savefig(file)
    plt.show()


def plot_twin_lagoon_output(data, file = "outputs/twin_lagoon_output.png", options = {"t_s": 0, "plot_t" : 180}):
    index = {"t": 0, "h_sl": [[1, 2],[3, 4]], "h_tb": [5,6], "DZ": [7,8,9],  "m": [10,11,12], "Q_s": [13,14], "Q_t": 15,
             "P": 16, "E": 17}

    f, axarr = plt.subplots(3, sharex = "all", sharey = "none", figsize=(6,6), dpi =200)

    axarr[0].plot(data[:,index["t"]]/3600 - options["t_s"], data[:,index["h_sl"][0][1]],label = " $WL_{o_1}$", linewidth = 0.45,color='blue')
    axarr[0].plot(data[:,index["t"]]/3600 - options["t_s"], data[:,index["h_sl"][1][1]],label = " $WL_{o_2}$", linewidth = 0.25,color='blue')
    axarr[0].plot(data[:,index["t"]]/3600 - options["t_s"], data[:,index["h_tb"][0]],label = " $WL_{HW}$",  linewidth = 0.25,color='red',ls='dashed')
    axarr[0].plot(data[:,index["t"]]/3600 - options["t_s"], data[:,index["h_tb"][1]],label = " $WL_{LW}$",  linewidth = 0.25,color='red',ls='dashed')
    axarr[0].legend( loc=2, ncol=3, borderaxespad=0., mode="expand", frameon=False, fontsize='x-small', labelspacing = 0.0)

    axarr[0].set_ylabel(r"$\eta (m)$ ", labelpad = 10 , fontsize = 10)
    # axarr[0].set_ylim([-10, 10])

    axarr[1].plot(data[:,index["t"]]/3600 - options["t_s"], data[:,index["Q_t"]]/1000,label = "$Q_{turb}$ ", linewidth = 0.45,color='red')
    axarr[1].plot(data[:,index["t"]]/3600 - options["t_s"], data[:,index["Q_s"]]/1000,label = "$Q_{sl}$ ", linewidth = 0.45,color='black',ls='dashed')
    axarr[1].plot(data[:, index["t"]] / 3600 - options["t_s"], data[:, index["Q_s"]] / 1000, label="$Q_{sl}$ ", linewidth=0.45, color='grey', ls='dashed')

    axarr[1].legend(ncol=3, loc=2, mode="expand", borderaxespad=0., frameon=False, fontsize='x-small', labelspacing = 0.0)
    axarr[1].set_ylabel(r"$Q  (m^3/s) \times 10^3$ ", labelpad = 5, fontsize = 10)
    # axarr[1].set_ylim([-80, 80])

    axarr[2].plot(data[:,index["t"]]/3600 - options["t_s"], data[:,index["P"]],label = "Power 0D", linewidth = 0.45,color='black')
    axarr[2].legend(ncol=3, loc=3, borderaxespad=0., frameon=False, mode="expand", fontsize='x-small', labelspacing = 0.0)
    axarr[2].set_ylabel(r"$P (MW)$ ", labelpad = 2, fontsize = 10)
    # axarr[2].set_ylim([-2000, 2000])

    plt.xlabel("Time $t$ (h)")
    plt.xlim(0, options["plot_t"])
    plt.setp(axarr[0].get_yticklabels(), fontsize=10)
    plt.setp(axarr[1].get_yticklabels(), fontsize=10)
    plt.setp(axarr[2].get_yticklabels(), fontsize=10)

    f.subplots_adjust(hspace=0)
    # plt.savefig(file)
    plt.show()

def add_point(array):
    # 1. 计算时间差值并找到新的开始点
    start_time_diff = array[0, 2] - array[0, 1]
    start_time = array[0, 0] - start_time_diff
    start_height = array[1, 0]

    # 2. 将新点添加到首端
    new_start_point = np.array([[start_time], [start_height]])
    array = np.hstack((new_start_point, array))

    # 3. 计算时间差值并找到新的结束点
    end_time_diff = array[0, -2] - array[0, -3]
    end_time = array[0, -1] + end_time_diff
    end_height = array[1, -1]

    # 4. 将新点添加到末端
    new_end_point = np.array([[end_time], [end_height]])
    array = np.hstack((array, new_end_point))

    return array


def extract_tidal_peaks_and_troughs(time, series):
    #     Use scipy to extract peaks
    peaks, _ = find_peaks(series)
    troughs, _ = find_peaks(-series)

    #     Ensure same index for peaks and trough
    array_size = min(peaks.shape[0], troughs.shape[0])
    peaks, troughs = peaks[:array_size], troughs[:array_size]

    peak_array = np.array([time[peaks], series[peaks]])
    trough_array = np.array([time[troughs], series[troughs]])

    peak_array = add_point(peak_array)
    trough_array = add_point(trough_array)
    
    return peak_array, trough_array


def extract_hill_chart(h_array, turbine_specs):
    """
    Produces the data to plot hill charts
    """
    p,q = [],[]
    for idx,h_val in enumerate(h_array):
        tmp = hill_chart_parametrisation_h(h_val, turbine_specs=turbine_specs)
        p.append(tmp[0])
        q.append(tmp[1])
    return np.array([h_array,np.array(p),np.array(q)]).T



def determine_mean_head(time_series, elevation_series, turbine_params=None):
    # Extract peaks and troughts
    peaks, troughs = extract_tidal_peaks_and_troughs(time_series, elevation_series)
    # Tidal range variation array
    H = peaks[1, 1:] - troughs[1, 1:]

    # Average time stamps between peaks and troughs
    t_p = np.mean(np.array([peaks[0, 1:], troughs[0, 1:]]), axis=0)

    # Set h_cap as mean tidal range
    if turbine_params == None:
        turbine_params = {"f_g": 50, "g_p": 95, "g": 9.807, "t_d": 7.35,
                          "h_cap": np.mean(H), "dens": 1025, "h_min": 1.00,
                          "eta": [0.93, 0.83], "options": 0}

    h_array = np.arange(0.1,10,0.1)
    hill_chart = extract_hill_chart(h_array, turbine_params)

    # summarise rated head and capacity
    rated_head, capacity = np.mean(H), np.max(hill_chart[:, 1])
    return rated_head, capacity


def ensure_dir(file_path):
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)