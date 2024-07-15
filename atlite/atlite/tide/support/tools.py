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


def plot_lagoon_output(data, file = "outputs/output.png", options = {"t_s": 0, "plot_t" : 240}):
    index = {"t": 0, "h_o": 1, "h_i": 2, "DZ": 3, "P": 4, "E": 5, "m": 6, "Q_t": 7, "Q_s": 8,
             "m_dt": 9, "m_t": 10, "f_r": 11}

    f, axarr = plt.subplots(3, sharex = "all", sharey = "none", figsize=(6,4), dpi =200)

    axarr[0].plot(data[:,index["t"]]/3600 - options["t_s"], data[:,index["h_i"]],label = " $\eta_{i}$", linewidth = 0.45,color='red')
    axarr[0].plot(data[:,index["t"]]/3600 - options["t_s"], data[:,index["h_o"]],label = " $\eta_{o}$",  linewidth = 0.25,color='blue',ls='dashed')
    axarr[0].legend(loc = 1, frameon=False, fontsize='x-small', labelspacing = 0.0)
    axarr[0].set_ylabel(r"$\eta (m)$ ", labelpad = 10 , fontsize = 10)
    axarr[0].set_ylim([-10, 10])

    axarr[1].plot(data[:,index["t"]]/3600 - options["t_s"], data[:,index["Q_t"]]/1000,label = "$Q_{turb}$ ", linewidth = 0.45,color='red')
    axarr[1].plot(data[:,index["t"]]/3600 - options["t_s"], data[:,index["Q_s"]]/1000,label = "$Q_{sl}$ ", linewidth = 0.45,color='black',ls='dashed')
    axarr[1].legend(loc = 1, frameon=False, fontsize='x-small', labelspacing = 0.0)
    axarr[1].set_ylabel(r"$Q  (m^3/s) \times 10^3$ ", labelpad = 5, fontsize = 10)
    axarr[1].set_ylim([-10, 10])

    axarr[2].plot(data[:,index["t"]]/3600 - options["t_s"], data[:,index["P"]],label = "Power 0D", linewidth = 0.45,color='black')
    axarr[2].legend(loc = 1, frameon=False, fontsize='x-small', labelspacing = 0.0)
    axarr[2].set_ylabel(r"$P (MW)$ ", labelpad = 2, fontsize = 10)
    axarr[2].set_ylim([-400, 400])

    plt.xlabel("Time $t$ (h)")
    plt.xlim(0, options["plot_t"])
    plt.setp(axarr[0].get_yticklabels(), fontsize=10)
    plt.setp(axarr[1].get_yticklabels(), fontsize=10)
    plt.setp(axarr[2].get_yticklabels(), fontsize=10)

    f.subplots_adjust(hspace=0)
    plt.savefig(file)


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


def extract_tidal_peaks_and_troughs(time, series):
    #     Use scipy to extract peaks
    peaks, _ = find_peaks(series)
    troughs, _ = find_peaks(-series)

    #     Ensure same index for peaks and trough
    array_size = min(peaks.shape[0], troughs.shape[0])
    peaks, troughs = peaks[:array_size], troughs[:array_size]

    peak_array = np.array([time[peaks], series[peaks]])
    trough_array = np.array([time[troughs], series[troughs]])

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
    #t_p = np.mean(np.array([peaks[0, 1:], troughs[0, 1:]]), axis=0)

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