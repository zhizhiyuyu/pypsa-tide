import numpy as np
from scipy.interpolate import interp1d
from .modules import input_0D, input_barrages, lagoon_operation
from .support import tools, auxillary
import os
import copy
import pickle
import scipy.optimize as so
import xarray as xr



class TidalLagoonSimulation:
    def __init__(self, eta, time_series, barrage, area=None, capacity=None):
        self.elevation_time_series = None
        self.area_elevation_curve = None
        self.turbine_params = None
        self.peaks = None
        self.troughs = None
        self.lagoon_status = None
        self.lagoon_control = None
        self.lagoon_params = None
        self.simulation_output = None
        self.output_energy = None
        self.variable_tide_limits = None
        self.barrage = barrage
        self.optimized_controls = None
        self.optimized_energy = None
        self.eta = eta
        self.time_series = time_series
        self.area = area
        self.total_capacity = capacity
        self.cost = None

    def calculate_elevations_and_area(self, area=50):
        # 计算elivation_time_series
        self.elevation_time_series = input_0D.read_cutout_elevations(self.eta, self.time_series)
        # 计算area_elevation_curve
        self.area_elevation_curve = input_0D.lagoon_idealised_area_case(area = area)

    def determine_head_and_turbine_params(self, Dt, turbine_params=None):
        t_sim = self.elevation_time_series.x[-1]
        H_m, P_c = tools.determine_mean_head(np.arange(0, t_sim, Dt), self.elevation_time_series(np.arange(0, t_sim, Dt)))
        if turbine_params is None:
            H_r = 0.8 * H_m
            self.turbine_params = {"f_g": 50, "g_p": 96, "g": 9.807, "t_d": 7.35, "t_cap": 20, "dens": 1025, "h_min": 1.00, "eta": [0.93, 0.83], "options": 1}
            return P_c
        elif turbine_params == "20MW":
            self.turbine_params = {"f_g": 50, "g_p": 95, "g": 9.807, "t_d": 7.35, "t_cap": 20, "dens": 1025, "h_min": 1.00, "eta": [0.93, 0.83], "options": 0}
            return self.turbine_params['t_cap']
        elif turbine_params == "30MW":
            self.turbine_params = {"f_g": 50, "g_p": 113, "g": 9.807, "t_d": 8.9, "t_cap": 30, "dens": 1025, "h_min": 1.00, "eta": [0.93, 0.83], "options": 0} 
            return self.turbine_params['t_cap']
        elif turbine_params == "40MW":
            self.turbine_params = {"f_g": 50, "g_p": 142, "g": 9.807, "t_d": 9.0, "t_cap": 40, "dens": 1025, "h_min": 1.00, "eta": [0.93, 0.83], "options": 0}
            return self.turbine_params['t_cap']
        else:
            self.turbine_params = turbine_params[0]
            return self.turbine_params['t_cap']
        
    def extract_peaks_and_troughs(self, limit_coefficient=1.2):
        time_elevation_series = np.array([self.elevation_time_series.x, self.elevation_time_series.y]).T
        self.peaks, self.troughs = tools.extract_tidal_peaks_and_troughs(time_elevation_series[:,0],time_elevation_series[:,1])
        self.variable_tide_limits = [interp1d(self.peaks[0,:], limit_coefficient * self.peaks[1,:]),
                                    interp1d(self.troughs[0,:], limit_coefficient * self.troughs[1,:])]
    def adaptive_operation(self):
        file_path = "./outputs/optimisation/control.p"
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"The file {file_path} does not exist! Please run optimization first.")   
        return pickle.load(open(file_path, "rb"))

    def ramp_simulate(self, Dt=120, cycles=10):
        # Pre-warming phase (simulation_ramp)
        simulation_ramp = {"t": cycles * 12.42 * 3600, "Dt": Dt, "start_t": 0}
        
        self.lagoon_status = lagoon_operation.tidal_lagoon_0d_model(simulation_ramp, self.elevation_time_series, self.area_elevation_curve,
                                                       self.lagoon_status, self.lagoon_control, self.lagoon_params)
                                                       
    def simulate(self, operation, Dt, cycles, barrage=None, warm=True, start_t=0., export_output=None, adaptive_operation=False):
        t_sim = cycles * 12.42 * 3600
        # Setup the control and parameters
        if isinstance(barrage, dict):
            self.lagoon_control, self.lagoon_params = input_barrages.input_barrage(barrage)  
        else:
            # Determine optimum rated head for turbines
            # input barrage is None, "20MW", "30MW","40MW" or a list of turbine parameters ([{"key":value}])
            cap = self.determine_head_and_turbine_params(Dt , turbine_params=barrage)
            N_t = int(self.total_capacity / cap)
            N_s = int(N_t / 2)
            self.lagoon_control, self.lagoon_params = input_barrages.input_predefined_barrage_specs(N_t, N_s, operation=operation, turbine_parameters=self.turbine_params)
            self.lagoon_params[0]['sluice_specs']['a_s'] = 324

        self.lagoon_status = input_barrages.initialise_barrage(len(self.lagoon_control))[0]
        self.lagoon_control = self.lagoon_control[0]
        self.lagoon_params = self.lagoon_params[0]

        if warm is not None:
            self.ramp_simulate(Dt=Dt, cycles=10)
        # Main simulation phase
        
        if operation == 'two-way-pump':
            limit_coefficient = 1.2
            self.extract_peaks_and_troughs(limit_coefficient=limit_coefficient)
        
        adaptive_control = None    
        if adaptive_operation:
            adaptive_control = self.adaptive_operation()
        
        simulation_main = {"t": t_sim, "Dt": Dt, "start_t": start_t}
        if export_output is not None:
            self.lagoon_status, self.simulation_output = lagoon_operation.tidal_lagoon_0d_model(simulation_main, self.elevation_time_series, self.area_elevation_curve, self.lagoon_status, self.lagoon_control, self.lagoon_params, export_output=export_output, adaptive_operation=adaptive_control, variable_tide_limits=self.variable_tide_limits)
        else: 
            self.lagoon_status = lagoon_operation.tidal_lagoon_0d_model(simulation_main, self.elevation_time_series, self.area_elevation_curve, self.lagoon_status, self.lagoon_control, self.lagoon_params, export_output=export_output, adaptive_operation=adaptive_control, variable_tide_limits=self.variable_tide_limits)  


    def plot_results(self, file_name, options=None):
        if self.simulation_output is None:
            raise ValueError("Please run the simulation and set export_output to be true first!")
        if options is None:
            options = {"t_s": 0, "plot_t": 180}
        tools.plot_lagoon_output(self.simulation_output, file= file_name, options=options)
    
    def run_simulation(self, operation, Dt, cycles, warm=None, start_t=0, export_output=True, adaptive_operation=False):

        area = self.area if self.area is not None else 50
        self.calculate_elevations_and_area(area)

        self.simulate(operation, Dt, cycles, barrage=self.barrage, warm=warm, start_t=start_t, export_output=export_output, adaptive_operation=adaptive_operation)
        
        self.total_capacity = self.total_capacity if self.total_capacity is not None else self.lagoon_control['N_t'] * self.lagoon_params['turbine_specs']['t_cap']
        
        self.output_energy = self.lagoon_status["E"]

        #self.cost = auxillary.determine_capital_cost_per_KW(self.total_capacity, area, lagoon_factor=0.75, efficiency=0.5, capacity_factor=0.20)


    def _gauge_func(self, x, simulation_main, elevation_time_series, area_elevation_curve, status, ctrl, params):
        status_copy = status.copy()
        ctrl["h_t"][0], ctrl["h_t"][1], ctrl["t_p"][0], ctrl["t_p"][1], ctrl["g_t"][0], ctrl["g_t"][1] = x
        return lagoon_operation.tidal_lagoon_0d_model(simulation_main, elevation_time_series, area_elevation_curve,
                                                      status_copy, ctrl, params, variable_tide_limits=self.variable_tide_limits)

    def _obj_func(self, x, simulation_main, elevation_time_series, area_elevation_curve, status, ctrl, params):
        status_copy = status.copy()
        ctrl["h_t"][0], ctrl["h_t"][1], ctrl["t_p"][0], ctrl["t_p"][1], ctrl["g_t"][0], ctrl["g_t"][1] = x
        cost = lagoon_operation.tidal_lagoon_0d_model(simulation_main, elevation_time_series, area_elevation_curve,
                                                      status_copy, ctrl, params, variable_tide_limits=self.variable_tide_limits)
        retuarn -cost["E"]

    def _optimize(self, Dt, cycles, bounds, x0, status, control, params, maxiter=50):
        control_storage = []
        optimised_energy = 0

        for i in range(0, cycles, 1):
            simulation_main = {"t": 12.42 * 3600 * 2, "Dt": Dt, "start_t": i * 12.42 * 3600}
            args = (simulation_main, self.elevation_time_series, self.area_elevation_curve, status, control.copy(), params)
            x_fin = so.minimize(self._obj_func, x0, args=args, method='L-BFGS-B', bounds=bounds, options={'eps': 5e-01, 'ftol': 1e-07, 'maxiter': maxiter})

            b = copy.deepcopy(control)
            b["h_t"][0], b["h_t"][1], b["t_p"][0], b["t_p"][1], b["g_t"][0], b["g_t"][1] = x_fin["x"]
            control_storage.append(b)
            optimised_energy += -x_fin["fun"]
            status.update(self._gauge_func(x_fin["x"], *args))

        return control_storage, optimised_energy

    def run_optimization(self, Dt, cycles, warm=True, bounds=[(0.0, 3.5),(0.0, 3.5), (0., 1.0),(0., 1.0), (0.0, 3.0),(0.0, 3.0)], total_capacity=None,operation=None, maxiter=50):     

        t_sim = cycles * 12.42 * 3600

        if self.barrage is not None:
            self.lagoon_control, self.lagoon_params = input_barrages.input_barrage(self.barrage)  
        else:
            # Determine optimum rated head for turbines
            self.determine_head_and_turbine_params(Dt, cycles)
            N_t = total_capacity / tools.determine_mean_head(np.arange(0, t_sim, Dt), self.elevation_time_series(np.arange(0, t_sim, Dt)), turbine_params=self.turbine_params)[1]
            N_s = N_t / 2
        
            self.lagoon_control, self.lagoon_params = input_barrages.input_predefined_barrage_specs(N_t, N_s, operation=operation, turbine_parameters=self.turbine_params)

        self.lagoon_status = input_barrages.initialise_barrage(len(self.lagoon_control))[0]
        self.lagoon_control = self.lagoon_control[0]
        self.lagoon_params = self.lagoon_params[0]

        if warm == True:
            self.ramp_simulate(Dt=Dt, cycles=10)
        
        if operation == 'two-way-pump':
            limit_coefficient = 1.2
            self.extract_peaks_and_troughs(limit_coefficient=limit_coefficient)

        x0 = np.array([self.lagoon_control["h_t"][0],self.lagoon_control["h_t"][1], self.lagoon_control["t_p"][0], self.lagoon_control["t_p"][1], self.lagoon_control["g_t"][0], self.lagoon_control["g_t"][1],])

        self.optimized_controls, self.optimized_energy = self._optimize(Dt, cycles, bounds, x0, self.lagoon_status, self.lagoon_control, self.lagoon_params, maxiter)
        
        with open ('outputs/optimisation/control.p', 'wb') as fp:
            pickle.dump(self.optimized_controls, fp)

def cap_series(interval, Dt,  lagoon_sim):
    # index
    index = {"t": 0, "h_o": 1, "h_i": 2, "DZ": 3, "P": 4, "E": 5, "m": 6, "Q_t": 7, "Q_s": 8,
            "m_dt": 9, "m_t": 10, "f_r": 11}

    lagoon_time_series_seconds = lagoon_sim.simulation_output[:, index["t"]]
    
    ratio = interval / Dt

    if ratio == int(ratio):
        ratio = int(ratio)
        matching_indices = np.arange(0, len(lagoon_time_series_seconds), step=ratio)
        cumulative_E_time_series = lagoon_sim.simulation_output[matching_indices, index["E"]]
        E_time_series = np.diff(cumulative_E_time_series, prepend=0)
        cap = E_time_series /lagoon_sim.total_capacity
        return cap 
    else:
        raise ValueError("The time intervals do not match exactly.")


def run_simulation_for_grid_point(eta_grid_point, area, capacity, time_series, barrage, operation, Dt, cycles):
    if np.isnan(eta_grid_point).any():
        return np.zeros_like(eta_grid_point)
    cutout_interval_seconds = time_series[1] - time_series[0]
    lagoon_sim = TidalLagoonSimulation(eta_grid_point, time_series, barrage=barrage, area=area, capacity=capacity)
    lagoon_sim.run_simulation(operation, Dt, cycles)
    cap = cap_series(cutout_interval_seconds, Dt, lagoon_sim)
    cost = lagoon_sim.cost
    return cap


def wrapper_run_simulation(eta, area, capacity, time_series, barrage, operation, Dt, cycles):
    return xr.apply_ufunc(
        run_simulation_for_grid_point,
        eta,
        area,
        capacity,
        input_core_dims=[["time"], [], []],
        output_core_dims=[["time"]],
        kwargs={"time_series": time_series, "barrage": barrage, "operation": operation, "Dt": Dt, "cycles": cycles},
        vectorize=True,
        dask="parallelized",
        output_dtypes=[np.float64]
    )

    #return xr.Dataset({'cap': cap, 'cost': cost})

def calculate_cost(capacity, area):
    return xr.apply_ufunc(
            auxillary.determine_capital_cost, 
            capacity,
            area,
            vectorize=True, 
            input_core_dims=[[], []], 
            output_core_dims=[[]],
            output_dtypes=[np.float64]

        )
