"""Basic injector sizing functions."""

import math


def orifice_area(diameter):
    """Cross-sectional area of a single orifice [m^2]."""
    return math.pi / 4 * diameter**2


def mass_flow_per_orifice(cd, orifice_diameter, delta_p, density):
    """Mass flow through one orifice.

    mdot = Cd · A · sqrt(2 · ρ · ΔP)

    Parameters
    ----------
    cd : float – discharge coefficient (typically 0.6–0.8)
    orifice_diameter : float – [m]
    delta_p : float – pressure drop across injector [Pa]
    density : float – propellant density [kg/m^3]

    Returns
    -------
    float – mass flow per orifice [kg/s]
    """
    A = orifice_area(orifice_diameter)
    return cd * A * math.sqrt(2 * density * delta_p)


def orifice_count(mdot, cd, orifice_diameter, delta_p, density):
    """Number of orifices needed for a given total mass flow.

    Parameters
    ----------
    mdot : float – total mass flow for this circuit (ox or fuel) [kg/s]
    cd : float – discharge coefficient
    orifice_diameter : float – [m]
    delta_p : float – pressure drop [Pa]
    density : float – propellant density [kg/m^3]

    Returns
    -------
    int – number of orifices (rounded up)
    """
    mdot_per = mass_flow_per_orifice(cd, orifice_diameter, delta_p, density)
    return math.ceil(mdot / mdot_per)


def pressure_drop(mdot, cd, num_orifices, orifice_diameter, density):
    """Injector pressure drop for given geometry and flow.

    ΔP = (mdot / (N · Cd · A))^2 / (2 · ρ)

    Returns
    -------
    float – pressure drop [Pa]
    """
    A = orifice_area(orifice_diameter)
    return (mdot / (num_orifices * cd * A))**2 / (2 * density)


def orifice_diameter_from_flow(mdot, cd, num_orifices, delta_p, density):
    """Size orifice diameter for a target pressure drop.

    Returns
    -------
    float – orifice diameter [m]
    """
    A = mdot / (num_orifices * cd * math.sqrt(2 * density * delta_p))
    return math.sqrt(4 * A / math.pi)


def check_pressure_drop_ratio(delta_p, chamber_pressure, low=0.15, high=0.25):
    """Check if injector ΔP is in the recommended range (15–25% of Pc).

    Returns
    -------
    dict with ratio, in_range, and message
    """
    ratio = delta_p / chamber_pressure
    in_range = low <= ratio <= high
    if ratio < low:
        msg = f"ΔP/Pc = {ratio:.1%} — too low, risk of combustion instability"
    elif ratio > high:
        msg = f"ΔP/Pc = {ratio:.1%} — high, increases feed system pressure requirements"
    else:
        msg = f"ΔP/Pc = {ratio:.1%} — within recommended range"
    return {"ratio": ratio, "in_range": in_range, "message": msg}
