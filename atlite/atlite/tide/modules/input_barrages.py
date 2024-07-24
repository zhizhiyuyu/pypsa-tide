import math
import numpy as np
import os

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


def initialise_lagoon_system(Barrage_Number):
    barrage_status = []
    for i in range(Barrage_Number):
        barrage_status.append({"m": [0, 0, 0], "m_t": [0., 0., 0.], "m_dt": [0., 0., 0.], "DZ" :[0., 0., 0.],
                               "f_r" :[0., 0., 0.], "Q_t": 0., "Q_s": [0., 0,], "P" : 0., "E": 0.})
    return barrage_status


def input_barrage(Barrage_file = None):
    if Barrage_file == None:
        return print ("No barrage or lagoon detected")
    else:
        f3 = open(Barrage_file, "r")
        # reading simulation points. specific for 0-D Model of Lagoons and barrage
        NumTB, NumSL, ASluice, Coed, gt, DZmin, holdt, hpump, tpump = [], [], [], [], [], [], [], [], [],
        DiameterTB, Capacity, Grid_f, Turb_Gp, Density = [], [], [], [], []
        ASluices, Coed_t , tr_lim, params , input_2D = [], [], [], [], []
        eta, turbine_parameterisation_option = [], []

        BarrageNo = int(f3.readline(5))
        f3.readline()
        for i in range(0,BarrageNo):
            for j in range(0, 3):
                f3.readline()          # General Operation info
            NumTB.append(float(f3.readline(5)))
            f3.readline()
            NumSL.append(float(f3.readline(5)))
            f3.readline()
            ASluice.append(float(f3.readline(5)))
            f3.readline()
            Coed.append(float(f3.readline(5)))
            f3.readline()
            gt.append([float(x) for x in f3.readline(10).split()])
            f3.readline()
            DZmin.append(float(f3.readline(5)))
            f3.readline()
            holdt.append([float(x) for x in f3.readline(10).split()])
            f3.readline()
            hpump.append(float(f3.readline(5)))
            f3.readline()
            tpump.append([float(x) for x in f3.readline(10).split()])
            f3.readline()
            tr_lim.append([float(x) for x in f3.readline(10).split()])
            f3.readline()
            for k in range(0, 3):
                f3.readline()   # Turbine Characteristics below
            DiameterTB.append(float(f3.readline(5)))
            f3.readline()
            Capacity.append(float(f3.readline(5)))
            f3.readline()
            Grid_f.append(float(f3.readline(5)))
            f3.readline()
            Turb_Gp.append(int(f3.readline(5)))
            f3.readline()
            Density.append(float(f3.readline(5)))
            f3.readline()
            eta.append([float(x) for x in f3.readline(10).split()])
            f3.readline()
            for k in range(0, 3):
                f3.readline()   # Turbine Characteristics below
            turbine_parameterisation_option.append(int(f3.readline(5)))
            f3.readline()

            # Preliminary calculations
            Grav = 9.807
            turbine_params = {"f_g": Grid_f[i], "g_p": Turb_Gp[i], "g": Grav, "t_d": DiameterTB[i],
                              "t_cap": Capacity[i], "dens": Density[i], "h_min": DZmin[i],
                              "eta": eta[i], "options": turbine_parameterisation_option[i]}

            Coed_t = turbine_parametrisation(DZmin[i], turbine_params)[1] / \
                     ((math.pi * (DiameterTB[i] / 2) ** 2) * math.sqrt(2 * Grav * DZmin[i]))

            sluice_params = {"a_s": ASluice[i], "c_d": Coed[i], "c_t": Coed_t, "g": Grav}

            params.append({"turbine_specs": turbine_params, "sluice_specs": sluice_params})

            input_2D.append({"h_t": holdt[i], "h_p": hpump[i], "t_p": tpump[i], "g_t": gt[i], "tr_l": tr_lim[i],
                             "N_t": NumTB[i], "N_s": NumSL[i]} )

        f3.readline()
        f3.close()
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
        #print("two-way")
        input_2D.append({"h_t": [3.0, 3.0], "h_p": 2.5, "t_p": [0., 0.], "g_t": [3.0, 3.0], "tr_l": [7, -6],
                         "N_t": NumTB, "N_s": NumSL})
    elif operation == "two-way-pump":
        #print("two-way-pump")
        input_2D.append({"h_t": [3., 3.], "h_p": 2.5, "t_p": [0.5, 0.5], "g_t": [3.0, 3.0], "tr_l": [7, -6],
                         "N_t": NumTB, "N_s": NumSL})

    return input_2D, params


def input_lagoon_system(Barrage_file=None):
    if Barrage_file == None:
        return print("No lagoon system detected")
    else:
        f3 = open(Barrage_file, "r")
        # reading simulation points. specific for 0-D Model of Lagoons and barrage
        NumTB, NumSL, ASluice, Coed, gt, DZmin, holdt, hpump, tpump = [], [], [], [], [], [], [], [], [],
        DiameterTB, Capacity, Grid_f, Turb_Gp, Density = [], [], [], [], []
        ASluices, tr_lim, params, turbine_params, sluice_params = [], [], [], [], [],
        input_2D, eta, turbine_parametrisation_option = [], [], []
        BarrageNo = int(f3.readline(5))
        f3.readline()
        for i in range(0, BarrageNo):
            for j in range(0, 3):
                f3.readline()  # General Operation info
            NumTB.append(float(f3.readline(5)))
            f3.readline()
            NumSL.append([int(x) for x in f3.readline(10).split()])
            f3.readline()
            ASluice.append(float(f3.readline(5)))
            f3.readline()
            Coed.append(float(f3.readline(5)))
            f3.readline()
            DZmin.append(float(f3.readline(5)))
            f3.readline()
            for k in range(0, 3):
                f3.readline()  # Turbine Characteristics below
            DiameterTB.append(float(f3.readline(5)))
            f3.readline()
            Capacity.append(float(f3.readline(5)))
            f3.readline()
            Grid_f.append(float(f3.readline(5)))
            f3.readline()
            Turb_Gp.append(int(f3.readline(5)))
            f3.readline()
            Density.append(float(f3.readline(5)))
            f3.readline()
            eta.append([float(x) for x in f3.readline(10).split()])
            f3.readline()
            for k in range(0, 3):
                f3.readline()   # Turbine Characteristics below
            turbine_parametrisation_option.append(int(f3.readline(5)))
            f3.readline()

            ASluices.append((np.array(NumSL[i]) * ASluice[i]).tolist())
            # Preliminary calculations
            Grav = 9.807

            if turbine_parametrisation_option[i] == 1: eta[i] = 1
            turbine_params = {"f_g": Grid_f[i], "g_p": Turb_Gp[i], "g": Grav, "t_d": DiameterTB[i],
                              "t_cap": Capacity[i],"dens": Density[i], "h_min": DZmin[i], "eta": eta[i] ,
                              "options": turbine_parametrisation_option[i]}

            Coed_t = turbine_parametrisation(DZmin[i], turbine_params)[1] / \
                     ((math.pi * (DiameterTB[i] / 2) ** 2) * math.sqrt(2 * Grav * DZmin[i]))

            sluice_params = {"a_s": ASluice[i], "c_d": Coed[i], "c_t": Coed_t, "g": Grav}
            input_2D.append({"N_t": NumTB[i], "N_s": NumSL[i]})

            params.append({"turbine_specs": turbine_params, "sluice_specs": sluice_params})

        f3.readline()
        f3.close()
        return input_2D, params



def input_predefined_system_specs(NumTB, NumSL):
    """
    :param NumTB: Number of turbines
    :param NumSL: Number of sluice gates
    :param operation: operation options: Ebb-only generation        ==> "ebb"
                                         Ebb-pump generation        ==> "ebb-pump"
                                         Two-way generation         ==> "two-way"
                                         Two-way-pumping generation ==> "two-way-pump"
    :return: control parameter array , turbine parameters
    """

    params, input_2D = [], []
    turbine_params = {"f_g": 50, "g_p": 95, "g": 9.807, "t_d": 7.35,
                      "t_cap": 20, "dens": 1025, "h_min": 1.00,
                      "eta": [0.93, 0.83], "options": 0}

    Coed_t = turbine_parametrisation(turbine_params["h_min"],
                                     turbine_params)[1] / ((math.pi * (turbine_params["t_d"] / 2)**2) *
                                                           math.sqrt(2 * turbine_params["g"] * turbine_params["h_min"]))

    sluice_params = {"a_s": 150, "c_d": 1.0, "c_t": Coed_t, "g": turbine_params["g"]}

    params.append({"turbine_specs": turbine_params, "sluice_specs": sluice_params})


    input_2D.append({"N_t": int(NumTB), "N_s": [int(NumSL/2), int(NumSL/2)]})

    return input_2D, params