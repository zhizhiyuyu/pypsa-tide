# calculate capacity per square kilometer
# theoratical resource
# PE = 0.5 * rho * g * Hiˆ2
import numpy as np
from scipy.signal import find_peaks

def extract_tidal_ranges(series):

    peaks, _ = find_peaks(series)
    troughs, _ = find_peaks(-series)

    all_indices = np.sort(np.concatenate([peaks, troughs]))

    tidal_ranges = np.diff(series[all_indices])

    return tidal_ranges

def calculate_potential_energy(series):
    rho = 1025  # kg/m³
    g = 9.81  # m/s²
    tidal_ranges = extract_tidal_ranges(series)
    total_pe = 0.5 * rho * g * np.sum(tidal_ranges ** 2)

    # kWh/m² GWh/km²
    total_pe_kwh = total_pe * 2.7778e-7

    return total_pe_kwh

def calculate_capacity(series):
    if np.isnan(series).any():
        return 0
    tidal_ranges = extract_tidal_ranges(series)
    DH = np.mean(abs(tidal_ranges))
    g = 9.807
    rho = 1025
    efficiency = 0.5
    capacity_factor = 0.20
    # MW/km2
    capacity_per_sqkm=( 2 * (0.5 * g * rho * DH * DH / 1e6*0.000277778) * efficiency / 12.42)/capacity_factor*1e6 
    return capacity_per_sqkm