"""First-order thermal / cooling estimates."""

import math

from lrd.utils import specific_gas_constant


def bartz_heat_flux(chamber_pressure, c_star, throat_diameter, T_c, gamma, mu, Pr, T_w=None):
    """Bartz correlation for convective heat transfer at the nozzle throat.

    q = h_g · (T_c - T_w)

    h_g from the simplified Bartz equation:
    h_g = (0.026 / d_t^0.2) · (μ^0.2 · cp / Pr^0.6) · (Pc / c*)^0.8 · (d_t / r_c)^0.1 · σ

    Uses a simplified form — adequate for first-pass design.

    Parameters
    ----------
    chamber_pressure : float – P_c [Pa]
    c_star : float – characteristic velocity [m/s]
    throat_diameter : float – [m]
    T_c : float – combustion temperature [K]
    gamma : float – ratio of specific heats
    mu : float – dynamic viscosity of combustion gas [Pa·s] (typical: 7e-5 to 1e-4)
    Pr : float – Prandtl number (typical: ~0.5 for combustion gases)
    T_w : float – wall temperature [K] (default: 0.5 * T_c)

    Returns
    -------
    dict with h_g [W/(m^2·K)] and q [W/m^2]
    """
    if T_w is None:
        T_w = 0.5 * T_c  # rough assumption

    d_t = throat_diameter
    g = gamma
    cp = g / (g - 1) * specific_gas_constant(23.0)  # approximate, using ~23 kg/kmol

    # Simplified Bartz
    h_g = (0.026 / d_t**0.2) * (mu**0.2 * cp / Pr**0.6) * (chamber_pressure / c_star)**0.8

    q = h_g * (T_c - T_w)

    return {"h_g": h_g, "heat_flux": q}


def film_cooling_fraction(heat_flux, coolant_cp, delta_T):
    """Estimate film cooling mass fraction needed to absorb a given heat flux.

    This is a very rough estimate — fraction of total mass flow diverted as film coolant.

    Parameters
    ----------
    heat_flux : float – [W/m^2]
    coolant_cp : float – specific heat of coolant [J/(kg·K)]
    delta_T : float – allowable coolant temperature rise [K]

    Returns
    -------
    float – approximate mass flux of coolant needed [kg/(m^2·s)]
    """
    return heat_flux / (coolant_cp * delta_T)


def regen_channel_sizing(heat_flux, coolant_mdot, coolant_cp, delta_T):
    """Basic regenerative cooling channel sizing.

    Parameters
    ----------
    heat_flux : float – local heat flux [W/m^2]
    coolant_mdot : float – coolant mass flow rate [kg/s]
    coolant_cp : float – coolant specific heat [J/(kg·K)]
    delta_T : float – allowable coolant temperature rise [K]

    Returns
    -------
    dict with required_area [m^2] — minimum channel wetted area to absorb heat
    """
    # Q = mdot · cp · ΔT  →  required area = Q / q
    Q = coolant_mdot * coolant_cp * delta_T  # total heat absorption capacity [W]
    required_area = Q / heat_flux if heat_flux > 0 else float("inf")
    return {"heat_absorption_capacity": Q, "required_area": required_area}
