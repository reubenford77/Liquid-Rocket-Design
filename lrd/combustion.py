"""Combustion chamber sizing functions."""

import math
import warnings


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


def chamber_dimensions(volume, contraction_ratio, throat_area, chamber_diameter=None):
    """Chamber diameter and cylindrical length from volume and geometry.

    Parameters
    ----------
    volume : float – chamber volume [m^3]
    contraction_ratio : float – A_c / A_t (typically 2–5), used when
        chamber_diameter is not provided.
    throat_area : float – [m^2]
    chamber_diameter : float or None – override chamber diameter [m].
        When provided, contraction_ratio is back-calculated.

    Returns
    -------
    dict with keys: diameter [m], length [m], area [m^2], contraction_ratio
    """
    if chamber_diameter is not None:
        Ac = math.pi / 4 * chamber_diameter ** 2
        dc = chamber_diameter
        if Ac < throat_area:
            raise ValueError(
                f"Chamber diameter {chamber_diameter} m gives chamber area "
                f"{Ac:.6f} m^2, which is smaller than throat area "
                f"{throat_area:.6f} m^2."
            )
        if contraction_ratio is not None:
            actual_cr = Ac / throat_area
            warnings.warn(
                f"Both chamber_diameter and contraction_ratio provided. "
                f"Using chamber_diameter={chamber_diameter:.4f} m "
                f"(implied CR={actual_cr:.2f})."
            )
    else:
        Ac = throat_area * contraction_ratio
        dc = math.sqrt(4 * Ac / math.pi)

    Lc = volume / Ac
    return {"diameter": dc, "length": Lc, "area": Ac, "contraction_ratio": Ac / throat_area}


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
