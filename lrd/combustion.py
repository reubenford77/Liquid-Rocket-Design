"""Combustion chamber sizing functions."""

import math


def chamber_volume(throat_area, l_star):
    """Chamber volume from throat area and characteristic length L*.

    V_c = A_t · L*

    Parameters
    ----------
    throat_area : float – [m^2]
    l_star : float – characteristic length [m] (typical: 0.8–3.0 m depending on propellant)

    Returns
    -------
    float – chamber volume [m^3]
    """
    return throat_area * l_star


def chamber_dimensions(volume, contraction_ratio, throat_area):
    """Chamber diameter and cylindrical length from volume and geometry.

    Parameters
    ----------
    volume : float – chamber volume [m^3]
    contraction_ratio : float – A_c / A_t (typically 2–5)
    throat_area : float – [m^2]

    Returns
    -------
    dict with keys: diameter [m], length [m], area [m^2]
    """
    Ac = throat_area * contraction_ratio
    dc = math.sqrt(4 * Ac / math.pi)
    Lc = volume / Ac
    return {"diameter": dc, "length": Lc, "area": Ac}


def stay_time(l_star, c_star, chamber_pressure, combustion_temp, mol_weight):
    """Residence time estimate (how long gas stays in the chamber).

    t_s = L* · P_c / (c* · R · T_c)  — approximate

    A simpler approximation: t_s ≈ L* / c*

    Parameters
    ----------
    l_star : float – [m]
    c_star : float – [m/s]
    chamber_pressure : float – [Pa] (unused in simple approx, kept for interface)
    combustion_temp : float – [K] (unused in simple approx)
    mol_weight : float – [kg/kmol] (unused in simple approx)

    Returns
    -------
    float – approximate stay time [s]
    """
    return l_star / c_star


def l_star_typical(propellant_key):
    """Return a typical L* value for common propellant combinations.

    Returns
    -------
    float – L* [m]
    """
    values = {
        "N2O/ETHANOL": 1.0,  # m — TODO: refine with test data
    }
    key = propellant_key.upper().replace(" ", "")
    if key not in values:
        return 1.0  # safe default
    return values[key]
