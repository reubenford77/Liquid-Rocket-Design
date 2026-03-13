"""Regenerative cooling thermal analysis."""

import math
import numpy as np

from lrd.utils import specific_gas_constant
from lrd.nozzle import mach_from_area_ratio


def bartz_sigma(T_c, T_w, gamma, mach):
    """Bartz sigma correction factor.

    Accounts for the temperature difference between gas total temperature
    and wall temperature on transport properties.

    Parameters
    ----------
    T_c : float – combustion (total) temperature [K]
    T_w : float – gas-side wall temperature [K]
    gamma : float – ratio of specific heats
    mach : float – local Mach number

    Returns
    -------
    float – sigma (dimensionless)
    """
    g = gamma
    M2 = mach**2
    T_ratio = T_w / (2 * T_c)
    bracket = T_ratio * (1 + (g - 1) / 2 * M2) + 0.5
    return bracket**(-0.68) * (1 + (g - 1) / 2 * M2)**(-0.12)


def bartz_heat_flux(Pc, c_star, dt, rc, T_c, gamma, mol_weight, mach, T_w=None):
    """Full Bartz correlation for convective heat transfer at a nozzle station.

    Parameters
    ----------
    Pc : float – chamber pressure [Pa]
    c_star : float – characteristic velocity [m/s]
    dt : float – throat diameter [m]
    rc : float – throat radius of curvature [m]
    T_c : float – combustion temperature [K]
    gamma : float – ratio of specific heats
    mol_weight : float – mean molecular weight of exhaust [kg/kmol]
    mach : float – local Mach number
    T_w : float or None – gas-side wall temperature [K] (default: 0.5 * T_c)

    Returns
    -------
    dict with h_g [W/(m^2·K)], heat_flux [W/m^2], T_aw [K], sigma
    """
    if T_w is None:
        T_w = 0.5 * T_c

    g = gamma
    R_sp = specific_gas_constant(mol_weight)
    cp = g / (g - 1) * R_sp

    # Gas transport property estimates
    mu = 46.6e-10 * mol_weight**0.5 * T_c**0.6
    Pr = 4 * g / (9 * g - 5)

    # Recovery factor and adiabatic wall temperature
    r_factor = Pr**0.33
    M2 = mach**2
    T_ratio_denom = 1 + (g - 1) / 2 * M2
    T_aw = T_c * (1 + r_factor * (g - 1) / 2 * M2) / T_ratio_denom

    # Sigma correction
    sigma = bartz_sigma(T_c, T_w, g, mach)

    # Full Bartz equation
    h_g = ((0.026 / dt**0.2)
           * (mu**0.2 * cp / Pr**0.6)
           * (Pc / c_star)**0.8
           * (dt / rc)**0.1
           * sigma)

    q = h_g * (T_aw - T_w)

    return {"h_g": h_g, "heat_flux": q, "T_aw": T_aw, "sigma": sigma}


def coolant_htc(Re, Pr, D_h, k_coolant):
    """Dittus-Boelter coolant-side heat transfer coefficient.

    Nu = 0.023 * Re^0.8 * Pr^0.4

    Parameters
    ----------
    Re : float – Reynolds number
    Pr : float – Prandtl number of coolant
    D_h : float – hydraulic diameter [m]
    k_coolant : float – coolant thermal conductivity [W/(m·K)]

    Returns
    -------
    float – h_c [W/(m^2·K)]
    """
    Nu = 0.023 * Re**0.8 * Pr**0.4
    return Nu * k_coolant / D_h


def regen_cooling_analysis(x_stations, r_stations, throat_index,
                           Pc, c_star, T_c, gamma, mol_weight,
                           dt, n_channels, channel_width, channel_height,
                           wall_thickness, wall_k,
                           coolant_mdot, coolant_rho, coolant_cp,
                           coolant_mu, coolant_k, T_coolant_in,
                           rc=None):
    """Station-by-station regenerative cooling analysis.

    Marches from the nozzle exit toward the injector (counter-flow) computing
    gas-side heat flux (Bartz), wall conduction, and coolant-side convection.

    Parameters
    ----------
    x_stations : array – axial positions [m] (x=0 at injector)
    r_stations : array – inner wall radius at each station [m]
    throat_index : int – index of the throat station
    Pc : float – chamber pressure [Pa]
    c_star : float – characteristic velocity [m/s]
    T_c : float – combustion temperature [K]
    gamma : float – ratio of specific heats
    mol_weight : float – mean molecular weight [kg/kmol]
    dt : float – throat diameter [m]
    n_channels : int – number of coolant channels
    channel_width : float – channel width [m]
    channel_height : float – channel height [m] (radial depth)
    wall_thickness : float – hot-gas-side wall thickness [m]
    wall_k : float – wall thermal conductivity [W/(m·K)]
    coolant_mdot : float – total coolant mass flow rate [kg/s]
    coolant_rho : float – coolant density [kg/m^3]
    coolant_cp : float – coolant specific heat [J/(kg·K)]
    coolant_mu : float – coolant dynamic viscosity [Pa·s]
    coolant_k : float – coolant thermal conductivity [W/(m·K)]
    T_coolant_in : float – coolant inlet temperature [K]
    rc : float or None – throat radius of curvature [m] (default: 0.75 * dt)

    Returns
    -------
    dict with arrays and scalars for the thermal solution
    """
    if rc is None:
        rc = 0.75 * dt

    n = len(x_stations)
    r_throat = dt / 2

    # Channel hydraulic properties
    A_channel = channel_width * channel_height
    P_wet = 2 * (channel_width + channel_height)
    D_h = 4 * A_channel / P_wet
    v_coolant = coolant_mdot / (n_channels * coolant_rho * A_channel)
    Re_coolant = coolant_rho * v_coolant * D_h / coolant_mu
    Pr_coolant = coolant_mu * coolant_cp / coolant_k

    # Coolant-side HTC (constant properties assumption)
    h_c = coolant_htc(Re_coolant, Pr_coolant, D_h, coolant_k)

    # Output arrays
    heat_flux = np.zeros(n)
    h_g_arr = np.zeros(n)
    T_wall_gas = np.zeros(n)
    T_wall_coolant = np.zeros(n)
    T_coolant_arr = np.zeros(n)
    T_aw_arr = np.zeros(n)
    pressure_drop = np.zeros(n)

    # March from nozzle exit (last index) toward injector (index 0)
    T_cool = T_coolant_in
    cum_dp = 0.0
    total_heat = 0.0

    # Blasius friction factor for turbulent flow
    f_friction = 0.316 / Re_coolant**0.25 if Re_coolant > 0 else 0.0

    for i in range(n - 1, -1, -1):
        # Local area ratio
        area_ratio = (r_stations[i] / r_throat)**2
        area_ratio = max(area_ratio, 1.0)  # clamp to >= 1

        # Local Mach number
        supersonic = i >= throat_index
        if area_ratio < 1.001:
            local_mach = 1.0
        else:
            local_mach = mach_from_area_ratio(gamma, area_ratio, supersonic=supersonic)

        # Iterate Bartz (2 passes for sigma convergence)
        T_w_guess = 0.5 * T_c
        for _ in range(3):
            bartz = bartz_heat_flux(Pc, c_star, dt, rc, T_c, gamma,
                                    mol_weight, local_mach, T_w=T_w_guess)
            # Overall heat transfer: gas → wall → coolant
            # 1/U = 1/h_g + t_wall/k_wall + 1/h_c
            R_total = 1.0 / bartz["h_g"] + wall_thickness / wall_k + 1.0 / h_c
            U = 1.0 / R_total
            q = U * (bartz["T_aw"] - T_cool)
            q = max(q, 0.0)  # no negative heat flux

            # Update wall temperature guess
            T_w_guess = bartz["T_aw"] - q / bartz["h_g"] if bartz["h_g"] > 0 else T_w_guess

        # Store results
        heat_flux[i] = q
        h_g_arr[i] = bartz["h_g"]
        T_aw_arr[i] = bartz["T_aw"]
        T_wall_gas[i] = bartz["T_aw"] - q / bartz["h_g"] if bartz["h_g"] > 0 else T_c
        T_wall_coolant[i] = T_cool + q / h_c if h_c > 0 else T_cool
        T_coolant_arr[i] = T_cool

        # Heat absorbed over this segment
        if i < n - 1:
            dx = abs(x_stations[i + 1] - x_stations[i])
        elif i > 0:
            dx = abs(x_stations[i] - x_stations[i - 1])
        else:
            dx = 0.0

        dA = 2 * math.pi * r_stations[i] * dx
        dQ = q * dA
        total_heat += dQ

        # Update coolant temperature
        if coolant_mdot > 0 and coolant_cp > 0:
            T_cool += dQ / (coolant_mdot * coolant_cp)

        # Pressure drop (Darcy-Weisbach)
        dp = f_friction * (dx / D_h) * 0.5 * coolant_rho * v_coolant**2
        cum_dp += dp
        pressure_drop[i] = cum_dp

    return {
        "x": x_stations,
        "heat_flux": heat_flux,
        "h_g": h_g_arr,
        "T_wall_gas": T_wall_gas,
        "T_wall_coolant": T_wall_coolant,
        "T_coolant": T_coolant_arr,
        "T_aw": T_aw_arr,
        "pressure_drop": pressure_drop,
        "total_heat": total_heat,
        "coolant_velocity": v_coolant,
        "Re_coolant": Re_coolant,
    }
