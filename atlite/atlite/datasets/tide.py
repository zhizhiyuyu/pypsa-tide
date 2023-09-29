import uptide
from datetime import datetime
from uptide import tidal_netcdf
import numpy as np
import xarray as xr
import pandas as pd
import logging

logger = logging.getLogger(__name__)

crs = 4326
features = {"tide_elevation": ["tide_elevation"]}


def get_data_tide_height(cutout, fes_path):
    # Extracting coordinates and adjusting longitude
    xs = (cutout.coords['x'] + 360) % 360
    ys = cutout.coords['y']
    ts = cutout.coords['time']

    if not cutout.dt[0].isdigit():
        dt_value = '1' + cutout.dt
    else:
        dt_value = cutout.dt
        
    dt = int(pd.Timedelta(dt_value).total_seconds())    
    # convert time to seconds
    time = ts.data.astype('M8[s]').astype('O')

    # Initializing tide calculation
    tide = uptide.Tides(['Q1', 'O1', 'P1', 'K1', 'N2', 'M2', 'S2', 'K2'])
    tide.set_initial_time(time[0])
    tnci = uptide.FES2014TidalInterpolator(tide, fes_path)

    # Calculating tide elevations
    eta_list = []
    for t in range(0,int((time[-1]-time[0]).total_seconds())+1, dt): # Extracting the numpy array from the DataArray
        tnci.set_time(t)  # Cast numpy int64 to regular int
        eta_t = [[tnci.get_val((lat, lon)) for lon in xs.data] for lat in ys.data]
        eta_list.append(eta_t)

    eta = np.array(eta_list)
    eta = np.squeeze(eta)
    return xr.DataArray(eta, coords=[("time", ts.data), ("y", ys.data), ("x", cutout.coords['x'].data)], name="tide_elevation")


def get_data(cutout, feature, tmpdir, **creation_parameters):
    """
    Get the tide height data.

    Parameters
    ----------
    cutout : atlite.Cutout
    feature : str
        Takes no effect, only here for consistency with other dataset modules.
    tmpdir : str
        Takes no effect, only here for consistency with other dataset modules.
    **creation_parameters :
        Must include `fes_path`.

    Returns
    -------
    xr.Dataset
    """
    if "fes_path" not in creation_parameters:
        logger.error('Argument "fes_path" not defined')
    path = creation_parameters["fes_path"]

    # assign time dimension even if not used
    return (
        get_data_tide_height(cutout, path)
        .to_dataset()
        .assign_coords(cutout.coords)
    )