from .modules import input_0D, input_barrages, lagoon_operation
from .support import tools
import numpy as np
import xarray as xr

class tidalLagoonSimulation:
    def __init__(self, time_series):
        """
        time_series: the time steps of water level data
        """
        self.dt = 120.
        self.time_series = time_series
        self.output = None       
        self.output_Energy = None

    def run_simulate(self, WL, total_capacity, area):
        """
        WL:  the water level time series
        total_capacity: the total_capacity of a plant (MW)
        area: area is an constant in each site
        """
        if len(WL) == 0:
            print(total_capacity, area)
            raise ValueError("Tidal range variation array H is empty.")
        
        elevation_time_series = input_0D.read_cutout_elevations(self.time_series, WL)
        #
        elevations = elevation_time_series(self.time_series)
        if len(elevations) == 0:
            raise ValueError("elevation_time_series returned an empty array.")
        #
        area_elevation_curve = input_0D.lagoon_idealised_area_case(area = area)
        operation = 'two-way'

        """
        Simulation setup
        """
        
        cycles = (self.time_series[-1] + self.dt ) /(12.42*3600) 
        t_sim = cycles * 12.42 * 3600
        H_m, P_c = tools.determine_mean_head(np.arange(0, t_sim, self.dt), elevation_time_series(np.arange(0, t_sim, self.dt)))
        H_r = 0.8 * H_m

        # Suggested turbine parameters
        turbine_params = {"f_g": 50, "g_p": 96, "g": 9.807, "t_d": 7.35,
                        "t_cap": P_c, "h_cap": H_r, "dens": 1025, "h_min": 1.00,
                        "eta": [0.93, 0.83], "options": 1}

        # For known total capacity specify
        N_t = total_capacity/tools.determine_mean_head(np.arange(0, t_sim, self.dt), elevation_time_series(np.arange(0, t_sim, self.dt)), turbine_params=turbine_params)[1]
        N_s = N_t/2
        
        lagoon_control, lagoon_params = input_barrages.input_predefined_barrage_specs(N_t,N_s,operation=operation,turbine_parameters=turbine_params) # initialise control / design parameteres - reading LagoonSpecs.dat
                
        # initialise operational status
        lagoon_status = input_barrages.initialise_barrage(len(lagoon_control))
        
        simulation_main = {"t": cycles * 12.42 * 3600, "Dt" : self.dt ,"start_t" : 0.}
        main, self.output = lagoon_operation.tidal_lagoon_0d_model(simulation_main, elevation_time_series, area_elevation_curve,
                                                            lagoon_status[0], lagoon_control[0], lagoon_params[0],
                                                            export_output= True, adaptive_operation=None,
                                                            variable_tide_limits=None)


            
def cap_series(interval, plant_sim, total_capacity):
    """
    interval: the difference of stime steps of WL time seriese
    Dt: the difference of time steps in simulation

    """
    # index
    index = {"t": 0, "h_o": 1, "h_i": 2, "DZ": 3, "P": 4, "E": 5, "m": 6, "Q_t": 7, "Q_s": 8,
            "m_dt": 9, "m_t": 10, "f_r": 11}

    plant_time_series_seconds = plant_sim.output[:, index["t"]]
    
    ratio = interval / plant_sim.dt

    if ratio == int(ratio):
        ratio = int(ratio)
        matching_indices = np.arange(0, len(plant_time_series_seconds), step=ratio)
        cumulative_E_time_series = plant_sim.output[matching_indices, index["E"]] 
        E_time_series = np.diff(cumulative_E_time_series, prepend=0)
        cap = E_time_series /total_capacity
        return cap 
    else:
        raise ValueError("The time intervals do not match exactly.")


def run_simulation_for_grid_point(WL, area,total_capacity, time_series):

    if np.isnan(area) or np.isnan(total_capacity) or total_capacity==0 or area ==0:
        return np.zeros_like(time_series)
    else:
        
        cutout_interval_seconds = time_series[1] - time_series[0]
        lagoon_sim = tidalLagoonSimulation(time_series)
        lagoon_sim.run_simulate(WL, total_capacity, area)
        cap = cap_series(cutout_interval_seconds, lagoon_sim, total_capacity)
        return cap




def wrapper_run_simulation(WL, area, total_capacity, time_series ):
    return xr.apply_ufunc(
        run_simulation_for_grid_point,
        WL,
        area,
        total_capacity,
        input_core_dims=[["time"],[],[]],
        output_core_dims=[["time"]],
        kwargs={"time_series": time_series},
        vectorize=True,
        dask="parallelized",
        output_dtypes=[np.float64]
    )

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
    capacity_per_sqkm=( 2 * (0.5 * g * rho * DH * DH / 1e6 * 0.000277778) * efficiency / 12.42)/capacity_factor * 1e6 
    return capacity_per_sqkm