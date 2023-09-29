import math

import numpy as np
from scipy.interpolate import interp1d
import pandas as pd

def read_outer_elevations(file):
    """
    Reading elevation time series.
    :param file: input_file with array of format (t, h)
    :return: elevation function
    """

    time_elevation_series = np.load(file)
    return interp1d(time_elevation_series[:,0],time_elevation_series[:,1])


def read_cutout_elevations(eta, time_series):
    return interp1d(time_series, eta)

def time_series(ds):
    #time_offset = 12.42 * 3600 * 10  # time shift
    time_cutout = ds.time
    timestamps = pd.to_datetime(time_cutout.values)
    time_diffs_ns = timestamps - timestamps[0]
    #time_final = time_diffs_ns.total_seconds().values - time_offset
    time_final = time_diffs_ns.total_seconds().values
    return  time_final

def read_twin_outer_elevations(file):
    """
    Reading elevation time series.
    :param file: input_file with array of format (t, h)
    :return: elevation function
    """
    time_elevation_series = np.load(file)

    HW = interp1d(time_elevation_series[:, 0], time_elevation_series[:, 2])
    LW = interp1d(time_elevation_series[:,0],time_elevation_series[:,1] )

    return {"LW": LW, "HW": HW}


def sinusoidal_outer_elevation(amplitude = 5.0, period = 12.42 * 3600):
    """
    Sinusoidal outer elevation based on time in h.
    :param amplitude: amplitude of tidal wave
    :param period: period of tidal wave, assuming M2 period of approximately 12.42 h
    :return:
    """
    omega = (2*np.pi) / period
    elevation = lambda t : amplitude * np.sin(omega * t)
    return elevation

def sinusoidal_outer_elevation_twin(amplitude = 5.0, period = 12.42 * 3600):
    """
    Sinusoidal outer elevation based on time in h.
    :param amplitude: amplitude of tidal wave
    :param period: period of tidal wave, assuming M2 period of approximately 12.42 h
    :return:
    """
    omega = (2*math.pi) / period
    HW = lambda t : amplitude * math.sin(omega * t)
    LW = lambda t : amplitude * math.sin(omega * t)
    return {"LW": LW, "HW": HW}

def extract_tidal_signal_from_UTM_coordinate(grid_file, hf_file, coords_xy, year):
    #TODO check if this works and whether it can be done in a quicker way

    import support.utm as utm
    import uptide.tidal_netcdf
    import datetime

    print("Producing input signal")
    utm_zone = 30
    utm_band = 'V'
    lat, lon = utm.to_latlon(coords_xy[0], coords_xy[1], utm_zone, utm_band)
    range_forcing_coords = ((lon-0.1, lon+0.1), (lat-0.1, lat+0.1))

    start_time = -50 * 24 * 3600
    end_time = 500 * 24 * 3600
    increment = int(300)
    n_increment = int((end_time - start_time) / increment +1)
    t_0d = np.linspace(start_time, end_time, int(n_increment))
    index = np.arange(n_increment)
    elev= np.zeros(n_increment)

    constituents = ['Q1', 'O1', 'P1', 'K1', 'N2', 'M2', 'S2', 'K2']
    tide = uptide.Tides(constituents)
    tide.set_initial_time(datetime.datetime(year, 1, 1, 0, 0))
    tnci = uptide.tidal_netcdf.OTPSncTidalInterpolator(tide, grid_file, hf_file,ranges=range_forcing_coords)

    for i in np.nditer(index):
        tnci.set_time(t_0d[i])
        try:
            elev[i] = tnci.get_val((lon, lat))  # Adding initial a correction depth for LAT
        except uptide.netcdf_reader.CoordinateError:
            elev[i] = 0.
    print("Input signal produced!")
    return  interp1d(t_0d,elev)

def read_area_elevation_curve(file, depth_correction  = 4, area_penalty = 1.0 , conversion=1.0):
    """
    Reading depth-area curve.
    :param file: input_file with array of format (h, area)
    :param depth_correction : amplitude correction in m (depending on datum)
    :return: h-area curve function
    """

    depth_area = np.load(file)
    return interp1d(depth_area[:,0] - depth_correction,depth_area[:,1] * area_penalty * conversion)

def read_twin_area_elevation_curve(file, depth_correction  = 4):
    """
    Reading depth-area curve.
    :param file: input_file with array of format (h, area)
    :param depth_correction : amplitude correction in m (depending on datum)
    :return: h-area curve function
    """

    depth_area = np.load(file)
    HW = interp1d(depth_area[:, 0] - depth_correction, depth_area[:, 1])
    LW = interp1d(depth_area[:, 0] - depth_correction, depth_area[:, 2])

    return {"LW":LW, "HW": HW}


def lagoon_system_idealised_area_case(area = 50):
    """
    Constant area case for twin-basin lagoon
    :param area: total inner lagoon area - assuming they are split evenly.
    :return: h-area curve function
    """
    LW = lambda h : area / 2
    HW = lambda h : area / 2

    return {"LW":LW, "HW": HW}

def lagoon_idealised_area_case(area = 61):
    """
    Constant area case for standard lagoon
    :param area: total inner lagoon area
    :return: h-area curve function
    """
    return lambda h: area