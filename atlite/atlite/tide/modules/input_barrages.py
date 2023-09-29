import math
import numpy as np
import os
import yaml
from .parameterisations import *

def initialise_barrage(Barrage_Number=1):
    """
    :param Barrage_Number: Number of barrages considered.
    :return:
    """
    barrage_status = []
    for i in range(Barrage_Number):
        QTurbines, QSluices, PowerGT, SumPower, Mode, ModeTime, ModeDuration, DZ, rampf, hdn0 =\
            0., 0., 0., 0., 1., 0., 0., 0., 0., 0.
        barrage_status.append({"m": Mode, "m_t": ModeTime, "m_dt": ModeDuration, "DZ": DZ, "f_r": rampf, "Q_t": QTurbines,
                               "Q_s": QSluices, "P": PowerGT, "E": SumPower,})
    return barrage_status


def input_barrage(barrage=None):
    
    lagoon = barrage['Lagoon']
    turbine = barrage['turbine']

    input_2D = [{
        'h_t': [lagoon['holdtime']['ebb'], lagoon['holdtime']['flood']],
        'h_p': lagoon['pumping_head'],
        't_p': [lagoon['pumping_duration']['ebb'], lagoon['pumping_duration']['flood']],
        'g_t': [lagoon['generation_time']['ebb'], lagoon['generation_time']['flood']],
        'tr_l': [lagoon['tidal_range_limits']['min'], lagoon['tidal_range_limits']['max']],
        'N_t': lagoon['turbine_no'],
        'N_s': lagoon['sluice_no']
    }]

    params = [{
        'turbine_specs': {
            'f_g': turbine['grid_frequency'],
            'g_p': turbine['generation_poles'],
            'g': 9.807,
            't_d': turbine['diameter'],
            't_cap': turbine['capacity'],
            'dens': turbine['fluid_density'],
            'h_min': lagoon['min_head'],
            'eta': [turbine['efficiency']['ebb'], turbine['efficiency']['flood']],
            'options': turbine['hill_chart']
        },
        'sluice_specs': {
            'a_s': lagoon['sluice_area'],
            'c_d': lagoon['COED'],
            'c_t': 1.3655909902568713,  # This value should be calculated based on other parameters
            'g': 9.807
        }
    }]

    return input_2D, params


def input_predefined_barrage_specs(NumTB, NumSL, operation='two-way', turbine_parameters=None):
    """
    :param NumTB: Number of turbines
    :param NumSL: Number of sluice gates
    :param operation: operation options: Ebb-only generation        ==> "ebb"
                                         Ebb-pump generation        ==> "ebb-pump"
                                         Two-way generation         ==> "two-way"
                                         Two-way-pumping generation ==> "two-way-pump"
    :param chart: hill chart options (based on capacity [0], rated head [1], idealised [2]
    :return: control parameter array , turbine parameters
    """

    params, input_2D = [], []
    if turbine_parameters==None:
        turbine_params = {"f_g": 50, "g_p": 95, "g": 9.807, "t_d": 7.35,
                          "t_cap": 20, "h_cap": 5., "dens": 1025, "h_min": 1.00,
                          "eta": [0.93, 0.83], "options": 0}
    else:
        turbine_params = turbine_parameters

    Coed_t = turbine_parametrisation(turbine_params["h_min"],
                                     turbine_params)[1] / ((math.pi * (turbine_params["t_d"] / 2)**2) *
                                                           math.sqrt(2 * turbine_params["g"] * turbine_params["h_min"]))

    sluice_params = {"a_s": 100, "c_d": 1.0, "c_t": Coed_t, "g": turbine_params["g"]}

    params.append({"turbine_specs": turbine_params, "sluice_specs": sluice_params})

    if operation == "ebb":
        input_2D.append({"h_t": [3.5, 0.], "h_p": 2.5, "t_p": [0., 0.], "g_t": [6.0, 6.0], "tr_l": [7, -6],
                         "N_t": NumTB, "N_s": NumSL})
    elif operation == "ebb-pump":
        input_2D.append({"h_t": [3.5, 0.], "h_p": 2.5, "t_p": [1.0, 0.], "g_t": [6.0, 6.0], "tr_l": [7, -6],
                         "N_t": NumTB, "N_s": NumSL})
    elif operation == "two-way":
        input_2D.append({"h_t": [3.0, 3.0], "h_p": 2.5, "t_p": [0., 0.], "g_t": [3.0, 3.0], "tr_l": [7, -6],
                         "N_t": NumTB, "N_s": NumSL})
    elif operation == "two-way-pump":
        input_2D.append({"h_t": [3., 3.], "h_p": 2.5, "t_p": [0.5, 0.5], "g_t": [3.0, 3.0], "tr_l": [7, -6],
                         "N_t": NumTB, "N_s": NumSL})

    return input_2D, params

