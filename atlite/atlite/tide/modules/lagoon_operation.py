import math
import numpy as np
from .parameterisations import turbine_parametrisation, turbine_generation, gate_sluicing, turbine_sluicing, \
    ramping_turbine

def lagoon_operation(h_i, h_o, t, status, control, turbine_specs, sluice_specs, flux_limiter=0.2):
    # --------------------------------------------------------------------------------------------------
    # turbine_specs["h_min"] = minimum head, z_up = upstream water level, z_dn = downstream water level
    # z_dn0 = downstream WL at previous timestep, t = time, ht = maximum holdtime, mod_0 = mode of operation at [i-1]
    # mod_t = time the current mode changed, status["m_dt"] = duration of current mode, f_g = grid frequency (turbine)
    # g_p = generator poles, g = gravitational acceleration, t_d = turbine diameter, t_cap = turbine capacity
    # dens = turbine fluid density (kg/m3) , control["N_t"] = number of turbines, a_s = total sluice area,
    # status["f_r"] = ramp function,
    # cd_t = turbine operating in sluice mode discharge coefficient.

    # Variables for pumping:
    # control["h_p"] = dictates which discharge will be calculated on the hill chart
    # t_pump = pumping duration in the current mode, h_boost = head difference desired for pumping

    mod_0 = status["m"]
    status["DZ"] = h_i - h_o

    # Main Operation Algorithm for a two-way operation
    if mod_0 == 9 and status["DZ"] > 0:
        status["m"], status["m_t"] = 10, t
        if control["t_p"][0] <= 0.2:
            status["m"], status["m_t"] = 1, t
    elif mod_0 == 10 and status["m_dt"] >= control["t_p"][0]:
        status["m"], status["m_t"] = 1, t
    elif mod_0 == 1 and control["h_t"][0] <= 0.2 and h_i > h_o:
        status["m"], status["m_t"] = 4, t
    elif mod_0 == 1 and status["m_dt"] < control["h_t"][0]:
        status["m"] = 1
    elif mod_0 == 1 and status["m_dt"] >= control["h_t"][0] and status["DZ"] > turbine_specs["h_min"]:
        status["m"], status["m_t"] = 2, t
    elif mod_0 == 2 and status["DZ"] < turbine_specs["h_min"] and status["m_dt"] > 0.25:
        status["m"], status["m_t"] = 4, t
    elif mod_0 == 2 and status["m_dt"] > control["g_t"][0]:
        status["m"], status["m_t"] = 3, t
    elif mod_0 == 3 and status["DZ"] >= turbine_specs["h_min"]:
        status["m"] = 3
    elif mod_0 == 3 and status["DZ"] < turbine_specs["h_min"]:
        status["m"], status["m_t"] = 4, t
    elif mod_0 == 4 and status["DZ"] < turbine_specs["h_min"] and h_i >= h_o:
        status["m"] = 4
    elif mod_0 == 4 and h_i < h_o:
        status["m"], status["m_t"] = 5, t
        if control["t_p"][1] <= 0.2:
            status["m"], status["m_t"] = 6, t
    elif mod_0 == 5 and status["m_dt"] > control["t_p"][1]:
        status["m"], status["m_t"] = 6, t
    elif mod_0 == 6 and control["h_t"][1] <= 0.2 and h_i < h_o:
        status["m"], status["m_t"] = 9, t
    elif mod_0 == 6 and -status["DZ"] < turbine_specs["h_min"] and status["m_dt"] < control["h_t"][1]:
        status["m"] = 6
    elif mod_0 == 6 and status["m_dt"] > control["h_t"][1] and -status["DZ"] > turbine_specs["h_min"]:
        status["m"], status["m_t"] = 7, t
    elif mod_0 == 7 and -status["DZ"] < turbine_specs["h_min"] and status["m_dt"] > 0.25:
        status["m"], status["m_t"] = 9, t
    elif mod_0 == 7 and status["m_dt"] > control["g_t"][1]:
        status["m"], status["m_t"] = 8, t
    elif mod_0 == 8 and -status["DZ"] > turbine_specs["h_min"]:
        status["m"] = 8
    elif mod_0 == 8 and -status["DZ"] < turbine_specs["h_min"]:
        status["m"], status["m_t"] = 9, t
    elif mod_0 == 9 and -status["DZ"] < turbine_specs["h_min"]:
        status["m"] = 9
    else:
        status["m"] = mod_0

    # Special cases
    # 1. Check for holding modes
    if mod_0 == 1 and turbine_specs["h_min"] > -status["DZ"] > 0 and status["m_dt"] > 2.:
        status["m"] = 6
    elif mod_0 == 6 and turbine_specs["h_min"] > status["DZ"] > 0 and status["m_dt"] > 2.:
        status["m"] = 1
    # 2. Sluice if in holding mode by accident
    if mod_0 == 1 and -status["DZ"] > 0 and status["m_dt"] > 0.1:
        status["m"], status["m_t"] = 9, t
    elif mod_0 == 6 and status["DZ"] > 0 and status["m_dt"] > 0.1:
        status["m"], status["m_t"] = 4, t

    status["m_dt"] = t - status["m_t"]

    # Ramp function set-up (primarily for stability and opening/closing hydraulic structures
    if status["m"] != mod_0:
        status["f_r"] = 0

    # Special case for generating/sluicing and sluicing modes
    if status["m"] == 4 and mod_0 == 3:
        status["f_r"] = 1
    elif status["m"] == 9 and mod_0 == 8:
        status["f_r"] = 1

    # If hydraulic structures are opening still, increase status["f_r"] based on a sine function
    if mod_0 == status["m"] and status["m_dt"] < 0.2 and status["f_r"] < 1.0:
        status["f_r"] = math.sin(math.pi / 2 * ((t - status["m_t"]) / 0.2))
    elif 0.4 > status["m_dt"] >= 0.2:  # second condition added just in case.
        status["f_r"] = 1.0

    # Special case - trigger end of pumping  # empirical
    if status["m"] == 5 and h_i <= (control["tr_l"][1] + 0.50):
        status["f_r"] = math.sin(math.pi / 2 * abs(h_i - control["tr_l"][0]) / 0.5)
        if status["f_r"] <= 0.3:
            status["m"], status["m_t"], status["m_dt"] = 6, t, 0.0

    if status["m"] == 10 and h_i >= (control["tr_l"][0] - 0.50):
        status["f_r"] = math.sin(math.pi / 2 * abs(control["tr_l"][1] - h_i) / 0.5)
        if status["f_r"] <= 0.3:
            status["m"], status["m_t"], status["m_dt"] = 1, t, 0.0

    # Special case - trigger end of pumping  # empirical
    if status["m"] == 5 and h_i <= control["tr_l"][1]:
            status["m"], status["m_t"], status["m_dt"] = 6, t, 0.0

    if status["m"] == 10 and h_i >= (control["tr_l"][0]):
            status["m"], status["m_t"], status["m_dt"] = 1, t, 0.0


    # ramping down for pumping stability
    if status["m"] == 5 and control["t_p"][1] - status["m_dt"] <= 0.2:
        status["f_r"] = math.sin(math.pi / 2 * ((control["t_p"][1] - status["m_dt"]) / 0.2))

    if status["m"] == 10 and control["t_p"][0] - status["m_dt"] <= 0.2:
        status["f_r"] = math.sin(math.pi / 2 * ((control["t_p"][0] - status["m_dt"]) / 0.2))

    # Individual mode operations#
    if status["m"] == 1 or status["m"] == 6:
        status["Q_t"], status["Q_s"], status["P"] = 0.0, 0.0, 0.0  # Holding

    if status["m"] == 2:  # Generating            HW -> LW
        status["P"] = status["f_r"] * control["N_t"] * turbine_specs["eta"][0] * turbine_parametrisation(abs(status["DZ"]),
                                                                                                         turbine_specs)[0]
        status["Q_t"] = -status["f_r"] * control["N_t"] * turbine_parametrisation(abs(status["DZ"]), turbine_specs)[1]
        status["Q_s"] = 0.0

    if status["m"] == 3:  # Generating/sluicing    HW -> LW
        status["P"] = control["N_t"] * turbine_parametrisation(abs(status["DZ"]), turbine_specs)[0] * turbine_specs["eta"][0]
        status["Q_t"] = - control["N_t"] * turbine_parametrisation(abs(status["DZ"]), turbine_specs)[1]
        status["Q_s"] = gate_sluicing(status["DZ"], status["f_r"], control["N_s"],
                                      status["Q_s"], sluice_specs, flux_limiter=flux_limiter)

    if status["m"] == 4:  # sluicing          HW -> LW
        status["P"] = 0.0
        status["Q_t"] = turbine_sluicing(status["DZ"], 1.0, control["N_t"],
                                         status["Q_t"], sluice_specs, turbine_specs, flux_limiter=flux_limiter)
        status["Q_s"] = gate_sluicing(status["DZ"], status["f_r"], control["N_s"],
                                      status["Q_s"], sluice_specs, flux_limiter=flux_limiter)

    if status["m"] == 5:  # Ebb pumping
        status["Q_t"] = - max(
            status["f_r"] * control["N_t"] * turbine_parametrisation(control["h_p"], turbine_specs)[1],
            0)  # Re-formulate the discharge coefficient for turbines! Introduce cd_t
        status["P"] = -(abs(status["Q_t"]) * turbine_specs["dens"] * turbine_specs["g"] *
                        abs(status["DZ"]) / (10 ** 6)) / min(max(0.4, 0.28409853 * np.log(abs(status["DZ"])) + 0.60270881), 0.9)
        status["Q_s"] = 0.0

    if status["m"] == 7:  # Generating             LW -> HW
        status["P"] = status["f_r"] * control["N_t"] * turbine_specs["eta"][1] * turbine_parametrisation(abs(status["DZ"]),
                                                                                                         turbine_specs)[0]
        status["Q_t"] = status["f_r"] * control["N_t"] * turbine_parametrisation(abs(status["DZ"]), turbine_specs)[1]
        status["Q_s"] = 0.0

    if status["m"] == 8:  # Generating / Sluicing LW -> HW
        status["P"] = control["N_t"] * turbine_specs["eta"][1] * turbine_parametrisation(abs(status["DZ"]),
                                                                                         turbine_specs)[0]
        status["Q_t"] = control["N_t"] * turbine_parametrisation(abs(status["DZ"]), turbine_specs)[1]
        status["Q_s"] = gate_sluicing(status["DZ"], status["f_r"], control["N_s"],
                                      status["Q_s"], sluice_specs, flux_limiter=flux_limiter)

    if status["m"] == 9:  # sluicing              LW -> HW
        status["P"] = 0.0
        status["Q_t"] = turbine_sluicing(status["DZ"], 1.0, control["N_t"],
                                         status["Q_t"], sluice_specs, turbine_specs, flux_limiter=flux_limiter)
        status["Q_s"] = gate_sluicing(status["DZ"], status["f_r"], control["N_s"],
                                      status["Q_s"], sluice_specs, flux_limiter=flux_limiter)

    if status["m"] == 10:  # Flood pumping
        status["Q_t"] = max(status["f_r"] * control["N_t"] * turbine_parametrisation(control["h_p"], turbine_specs)[1],
                            0.0)  # Re-formulate the discharge coefficient for turbines! Introduce cd_t
        status["P"] = -(
            abs(status["Q_t"]) * turbine_specs["dens"] * turbine_specs["g"] * abs(status["DZ"]) / (10 ** 6)) / min(
            max(0.3, 0.28409853 * np.log(abs(status["DZ"])) + 0.60270881), 0.8)
        status["Q_s"] = 0.0



    return status


def two_lagoon_system_operation(h_sluice, h_turbine, t, status, control, turbine_specs, sluice_specs,
                                turb_buffer=0.5, flux_limiter=0.2):
    """
    :param h_sluice: water elevation $\eta$ array for each of the sluice gate sections.
                    (Convention [$\eta_{inner}$, $\eta_{outer}$])
    :param h_turbine: water elevation $\eta$ array for each of the turbine section.
                    (Convention [$\eta_{inner}$, $\eta_{outer}$])
    :param t: time in hours
    :param status: status of operation of current power plant (see initialisation function)
    :param control: control parameters (Turbine and sluice number)
    :param f_g: grid frequency (Hz)
    :param g_p: generator poles
    :param g: gravity ($m^2/s$)
    :param t_d: turbine diameter (m)
    :param t_cap: capacity of turbines (MW)
    :param dens: fluid density (kg/m$^3$)
    :param a_s: sluice gate area
    :param c_d: sluice discharge coefficient
    :param cd_t: turbine discharge coefficient
    :param turb_buffer: commencing generation at $h_{min}$+ turb_buffer (m)
    :param flux_limiter: setting a flux limiter due to the way head differences are sampled.
    :return:
    """
    status["DZ"][0] = h_sluice[0][0] - h_sluice[0][1]
    status["DZ"][1] = h_sluice[1][0] - h_sluice[1][1]
    status["DZ"][2] = h_turbine["h_HW"] - h_turbine["h_LW"]

    mod_0 = status["m"]

    i = 0   # LW sluices

    # Conditions
    if mod_0[i] == 0 and status["DZ"][i] > 0.0:
        status["m"][i], status["m_t"][i], status["f_r"][i] = 1, t, 0.
    elif mod_0[i] == 1 and status["DZ"][i] <= 0.0:
        status["m"][i], status["m_t"][i], status["f_r"][i] = 0, t, 0.

    status["m_dt"][i] = t - status["m_t"][i]

    # If hydraulic structures are opening still, increase f_r based on a sine function
    if mod_0[i] == status["m"][i] and status["m_dt"][i] < 0.3 and status["f_r"][i] < 1.0:
        status["f_r"][i] = math.sin(math.pi / 2 * ((t - status["m_t"][i]) / 0.3))
    elif status["m_dt"][i] >= 0.3 and status["m_dt"][i] < 0.4:  # second condition added just in case.
        status["f_r"][i] = 1.0

    # Calculate LW sluice fluxes
    if status["m"][i] == 0: status["Q_s"][i] = 0.
    if status["m"][i] == 1: status["Q_s"][i] = gate_sluicing(status["DZ"][i], status["f_r"][i], control["N_s"][i],
                                                             status["Q_s"][i], sluice_specs, flux_limiter=flux_limiter)

    i = 1   # HW sluices

    # Conditions
    if mod_0[i] == 0 and status["DZ"][i] < 0.0:
        status["m"][i], status["m_t"][i], status["f_r"][i] = 1, t, 0.
    elif mod_0[i] == 1 and status["DZ"][i] >= 0.0:
        status["m"][i], status["m_t"][i], status["f_r"][i] = 0, t, 0.

    status["m_dt"][i] = t - status["m_t"][i]

    # If hydraulic structures are opening still, increase f_r based on a sine function
    if mod_0[i] == status["m"][i] and status["m_dt"][i] < 0.3 and status["f_r"][i] < 1.0:
        status["f_r"][i] = math.sin(math.pi / 2 * ((t - status["m_t"][i]) / 0.3))
    elif status["m_dt"][i] >= 0.3 and status["m_dt"][i] < 0.4:  # second condition added just in case.
        status["f_r"][i] = 1.0

    # Calculate HW sluice fluxes
    if status["m"][i] == 0: status["Q_s"][i] = 0.
    if status["m"][i] == 1: status["Q_s"][i] = gate_sluicing(status["DZ"][i], status["f_r"][i], control["N_s"][i],
                                                             status["Q_s"][i], sluice_specs, flux_limiter=flux_limiter)

    i = 2  # Turbines

    # Conditions
    if mod_0[i] == 0 and status["DZ"][i] >= turbine_specs["h_min"] + turb_buffer:
        status["m"][i], status["m_t"][i], status["f_r"][i] = 1, t, 0.
    elif mod_0[i] == 1 and status["DZ"][i] < turbine_specs["h_min"]:
        status["m"][i], status["m_t"][i], status["f_r"][i] = 0, t, 0

    status["m_dt"][i] = t - status["m_t"][i]

    # If hydraulic structures are opening still, increase f_r based on a sine function
    if mod_0[i] == status["m"][i] and status["m_dt"][i] < 0.3 and status["f_r"][i] < 1.0:
        status["f_r"][i] = math.sin(math.pi / 2 * ((t - status["m_t"][i]) / 0.3))
    elif status["m_dt"][i] >= 0.3 and status["m_dt"][i] < 0.4:  # second condition added just in case.
        status["f_r"][i] = 1.0

    # Calculate Turbine fluxes
    if status["m"][i] == 0: status["Q_t"], status["P"] = 0., 0.
    if status["m"][i] == 1: status["Q_t"], status["P"] = turbine_generation(status["DZ"][i], status["f_r"][i],
                                                                        control["N_t"], turbine_specs)

    return status


def lagoon(t, Dt, h_i, h_o, status, control, params, boundaries):
    """
    Calculates based on the head differences and the operational conditions the fluxes of a normal tidal lagoon/barrage
    by calling lagoon__operation. Exports to file and imposes the boundary conditions.
    :param t: time (s)
    :param Dt: timestep (s)
    :param h_i: upstream water level (m)
    :param h_o: downstream water level (m)
    :param status: current status parameters of lagoon/barrage
    :param control: control parameters imposed for the operation
    :param boundaries: hydraulic structure location
    :param file: output file
    :return: lagoon_operation(h_i, h_o, t, status, control, turbine_specs, sluice_specs):
    """

    status = lagoon_operation(h_i, h_o,  t / 3600, status, control, params["turbine_specs"], params["sluice_specs"])
    status["E"] += status["P"] * Dt / 3600

    if boundaries != "None":
        boundaries["tb_o"].assign(status["Q_t"])
        boundaries["tb_i"].assign(-status["Q_t"])
        boundaries["sl_o"].assign(status["Q_s"])
        boundaries["sl_i"].assign(-status["Q_s"])

    return np.hstack(([t, h_o, h_i], [status["DZ"], status["P"], status["E"], status["m"], status["Q_t"], status["Q_s"],
                                      status["m_dt"], status["m_t"], status["f_r"]],))


def two_lagoon_system(t, dt, h_sl, h_tb, status, control, params, boundaries, transition_turbines = None):
    """
    Calculates based on the head differences and the operational conditions the fluxes of a tidal lagoon system
    by calling tidal_lagoon_system_operation. Exports to file and imposes the boundary conditions.
    :param t:  Time in sec
    :param dt: Timestep
    :param h_sl: average water level opposite sluice gates (Convention -- inner first, outer second)
    :param h_tb: average water level opposite turbines (Convention -- based on direction of turbines (assuming one-way))
    :param status: dictionary detailing the state of the lagoon operation
    :param control: control terms, in this case Number of turbines and sluice gates
    :param params: parameters associated with the operation
    :param boundaries: dictionary of boundaries specific to a two-basin system
    :param file: output file
    :param firedrake: flag on whether we are running a firedrake/Thetis simulation
    :return:
    """

    if transition_turbines is not None: control["N_t"] = ramping_turbine(transition_turbines,t/3600)
    status = two_lagoon_system_operation(h_sl, h_tb, t / 3600, status, control, params["turbine_specs"],
                                         params["sluice_specs"])
    status["E"] += status["P"] * dt / 3600

    # file = open("output.dat", "a")
    # file.write(" {:6.3f} ".format(t))
    # file.write(
    #     " {:9.3f} {:9.3f} {:6.3f} {:4.0f} {:9.3f}".format(*h_sl[0], status["DZ"][0], status["m"][0], status["Q_s"][0]))
    # file.write(
    #     " {:9.3f} {:9.3f} {:6.3f} {:4.0f} {:9.3f}".format(*h_sl[1], status["DZ"][1], status["m"][1], status["Q_s"][1]))
    # file.write(
    #     " {:9.3f} {:9.3f} {:6.3f} {:4.0f} {:9.3f}".format(h_tb["h_HW"], h_tb["h_LW"], status["DZ"][2], status["m"][2],
    #                                                       status["Q_t"]))
    # file.write(" {:9.3f} {:9.3f} ".format(status["P"], status["E"]))
    # file.write(" {:9.3f} ".format(control["N_t"]))
    # file.write(str(status["f_r"]) + "\n")
    # file.flush()


    if boundaries != "None":
        boundaries["tb_dn"].assign(status["Q_t"])
        boundaries["tb_up"].assign(-status["Q_t"])
        boundaries["sl_o_LW"].assign(status["Q_s"][0])
        boundaries["sl_i_LW"].assign(-status["Q_s"][0])
        boundaries["sl_o_HW"].assign(status["Q_s"][1])
        boundaries["sl_i_HW"].assign(-status["Q_s"][1])
    else:
        h_tb["h_LW"] = h_sl[0][0]
        h_tb["h_HW"] = h_sl[1][0]

    return np.hstack(
        ([t], h_sl[0], h_sl[1], list(h_tb.values()), status["DZ"], status["m"], status["Q_s"], status["Q_t"],
         status["P"], status["E"]))


def tidal_lagoon_system_0d_model(simulation, elev_time, area_elev, status, control, params, export_output=None,
                                 adaptive_operation = None):

    h_inner = {"h_LW": float(elev_time["LW"](simulation["start_t"]) + status["DZ"][0]),
               "h_HW": float(elev_time["HW"](simulation["start_t"]) + status["DZ"][1])}
    status["E"], cycle ,   = 0., -1,
    turbine_transition ={"initial": control["N_t"], "final": control["N_t"], "start_t":simulation["start_t"]}

    for i in range(len(status["m_t"])):
        status["m_t"][i] += simulation["start_t"] / 3600.

    if export_output is not None: data = np.empty((int(simulation["t"] / simulation["Dt"]), 18), dtype=object)

    for step in range(0, int(simulation["t"] / simulation["Dt"]), 1):
        h_sl = [[float(h_inner["h_LW"]), float(elev_time["LW"](step * simulation["Dt"] + simulation["start_t"]))],
                [float(h_inner["h_HW"]), float(elev_time["HW"](step * simulation["Dt"] + simulation["start_t"]))]]


        if adaptive_operation is not None:
            if cycle is not int((step * simulation["Dt"] + simulation["start_t"])/(12.42*3600)):
                turbine_transition["initial"], turbine_transition["start_t"] = control["N_t"], step * simulation["Dt"]/3600
                control =  adaptive_operation[int((step * simulation["Dt"])/(12.425*3600))]
                turbine_transition["final"] = control["N_t"]
                cycle += 1

        if export_output is not None:
            data[step] = two_lagoon_system(step * simulation["Dt"] + simulation["start_t"], simulation["Dt"], h_sl,
                                           h_inner,status, control, params, boundaries="None",
                                           transition_turbines = turbine_transition )


        else:
            two_lagoon_system(step * simulation["Dt"] + simulation["start_t"], simulation["Dt"], h_sl, h_inner,
                              status, control, params, boundaries="None", transition_turbines = turbine_transition)

        h_inner["h_LW"] += float((-status["Q_t"] + status["Q_s"][0]) * simulation["Dt"] / \
                                 (area_elev["LW"](h_inner["h_LW"]) * (10 ** 6)))
        h_inner["h_HW"] += float((status["Q_t"] + status["Q_s"][1]) * simulation["Dt"] / \
                                 (area_elev["HW"](h_inner["h_HW"]) * (10 ** 6)))

    if export_output is not None:
        return status, data
    else:
        return status


def tidal_lagoon_0d_model(simulation, elev_time, area_elev, status, control, params, export_output=None,
                          adaptive_operation=None, variable_tide_limits=None):
    h_inner = elev_time(simulation["start_t"]) + status["DZ"]
    h_inner += (status["Q_t"] + status["Q_s"]) * simulation["Dt"]/  (area_elev(h_inner) * (10 ** 6))
    status["E"], cycle = 0., -1

    if export_output is not None: data = np.empty((int(simulation["t"] / simulation["Dt"]), 12), dtype=object)
    f =  open("output.dat", "a")
    for step in range(0, int(simulation["t"] / simulation["Dt"]), 1):

        # Change control parameters if desirable
        if adaptive_operation is not None:
            if cycle is not int((step * simulation["Dt"] + simulation["start_t"])/(12.42*3600)):
                control =  adaptive_operation [int((step * simulation["Dt"] + simulation["start_t"])/(12.42*3600))]
                cycle += 1

        # Restrict pumping to tidal limits
        if variable_tide_limits is not None:
            control["tr_l"][0] = variable_tide_limits[0](step * simulation["Dt"])
            control["tr_l"][1] = variable_tide_limits[1](step * simulation["Dt"])


        if export_output is not None:
            data[step] = lagoon(step * simulation["Dt"] + simulation["start_t"], simulation["Dt"], h_inner,
                                elev_time(step * simulation["Dt"] + simulation["start_t"]), status, control, params,
                                boundaries="None")
            for item in data[step]:
                f.write(" {:8.3f} ".format(item))
            f.write("\n")
            f.flush()

        else:
            lagoon(step * simulation["Dt"] + simulation["start_t"], simulation["Dt"], h_inner,
                   elev_time(step * simulation["Dt"] + simulation["start_t"]), status, control, params, boundaries="None")

        h_inner += (status["Q_t"] + status["Q_s"]) * simulation["Dt"] / (area_elev(h_inner) * (10 ** 6))

    if export_output is not None:
        f.write("\n")
        return status, data
    else:
        return status
