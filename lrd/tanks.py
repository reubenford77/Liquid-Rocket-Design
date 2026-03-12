"""Tank sizing and feed system functions."""

import math


def tank_volume(propellant_mass, density, ullage_fraction=0.05):
    """Tank volume including ullage.

    Parameters
    ----------
    propellant_mass : float – [kg]
    density : float – propellant density [kg/m^3]
    ullage_fraction : float – extra volume fraction (default 5%)

    Returns
    -------
    float – tank volume [m^3]
    """
    return (propellant_mass / density) * (1 + ullage_fraction)


def tank_mass_thin_wall(volume, pressure, material_yield, safety_factor=1.5, density_material=2700):
    """Estimate tank mass for a thin-wall spherical pressure vessel.

    Parameters
    ----------
    volume : float – internal volume [m^3]
    pressure : float – internal pressure [Pa]
    material_yield : float – yield strength [Pa] (e.g. Al 6061-T6: ~276 MPa)
    safety_factor : float
    density_material : float – wall material density [kg/m^3] (default: aluminum)

    Returns
    -------
    dict with keys: mass [kg], radius [m], wall_thickness [m]
    """
    # Sphere radius from volume
    r = (3 * volume / (4 * math.pi)) ** (1 / 3)
    # Thin-wall hoop stress: t = P·r / (2·σ)
    t = (pressure * r * safety_factor) / (2 * material_yield)
    # Shell mass
    mass = density_material * 4 * math.pi * r**2 * t
    return {"mass": mass, "radius": r, "wall_thickness": t}


def feed_pressure(chamber_pressure, injector_dp, line_losses=0.0):
    """Required tank / pump outlet pressure.

    Parameters
    ----------
    chamber_pressure : float – [Pa]
    injector_dp : float – injector pressure drop [Pa]
    line_losses : float – friction / minor losses in feed lines [Pa]

    Returns
    -------
    float – required feed pressure [Pa]
    """
    return chamber_pressure + injector_dp + line_losses


def blowdown_ratio(initial_pressure, final_pressure):
    """Blowdown pressure ratio for a pressure-fed system.

    Returns
    -------
    float – P_initial / P_final
    """
    return initial_pressure / final_pressure


def propellant_mass(mdot, burn_time):
    """Total propellant mass for a given flow rate and burn duration.

    Returns
    -------
    float – propellant mass [kg]
    """
    return mdot * burn_time
