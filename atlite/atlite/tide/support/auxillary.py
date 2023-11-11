import math
from scipy.signal import find_peaks
from ..modules.parameterisations import hill_chart_parametrisation_h


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
    t_p = np.mean(np.array([peaks[0, 1:], troughs[0, 1:]]), axis=0)

    # Set h_cap as mean tidal range
    if turbine_params == None:
        turbine_params = {"f_g": 50, "g_p": 95, "g": 9.807, "t_d": 7.35,
                          "h_cap": np.mean(H), "dens": 1025, "h_min": 1.00,
                          "eta": [0.93, 0.83], "options": 0}

    h_array = np.arange(0.1, 10, 0.1)
    hill_chart = extract_hill_chart(h_array, turbine_params)

    # summarise rated head and capacity
    rated_head, capacity = np.mean(H), np.max(hill_chart[:, 1])
    return rated_head, capacity



def determine_capital_cost(capacity, area_required, lagoon_factor=0.75, efficiency=0.5, capacity_factor=0.20):
    if area_required == 0:
        return 0   
    turbine_cost, caisson_cost = 1000, 1000
    mean_depth = 10
    radius = math.sqrt(area_required / math.pi)
    impoundment_length = 2 * math.pi * radius
    embankment_cost_per_m3 = 100  # typically 94 pounds per m3

    impoundment_cost = lagoon_factor * impoundment_length * (
    mean_depth * mean_depth * 7.5 / 2) * embankment_cost_per_m3 / (
                       capacity * 1000)  # 7.5 factor for ratio between embankment depth and width. (gentle enough slope)
    lagoon_cost = turbine_cost + caisson_cost + impoundment_cost
    #print("Capital cost in pounds per KW:", lagoon_cost)
    #print("Overall tidal range structure cost:", lagoon_cost * capacity * 1000 / 1e9, "billions")
    pound_to_euro = 1.17
    return lagoon_cost * pound_to_euro

def determine_area(capacity, tidal_range, efficiency=0.5, capacity_factor=0.20):

    grav, dens = 9.807, 1025
    head_difference_factor = 1.25

    DH = tidal_range * head_difference_factor
    area_required = capacity_factor * capacity / ((0.5 * grav * dens * DH * DH / 1e6*0.000277778) * efficiency / 12.42)
    print("Area (km2): ", area_required/1e6)
    return area_required/1e6

from ..modules.input_0D import read_area_elevation_curve
from scipy.integrate import quad
import numpy as np

def lagoon_volume_calculation(file, depth_correction=4, elevation=0):
    depth_area = np.load(file)
    f = read_area_elevation_curve(file,depth_correction=depth_correction, conversion=1e6)
    volume = quad(f,depth_area[0,0],elevation, limit=50,)
    return volume


#a = lambda x: 1e6
#volume = quad(a,-10,0, maxp1=40, limlst=40)
#vol = lagoon_volume_calculation('../inputs/area_swansea.npy')
#print(volume, vol, volume[0]-vol[0])

