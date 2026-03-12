"""Isentropic nozzle design functions."""

import math
import numpy as np
from scipy.optimize import brentq

from lrd.utils import G0, specific_gas_constant


def throat_area(mdot, chamber_pressure, c_star):
    """Throat area from mass flow, chamber pressure, and characteristic velocity.

    A_t = mdot * c_star / P_c

    Parameters
    ----------
    mdot : float – total mass flow rate [kg/s]
    chamber_pressure : float – P_c [Pa]
    c_star : float – characteristic velocity [m/s]

    Returns
    -------
    float – throat area [m^2]
    """
    return mdot * c_star / chamber_pressure


def exit_area(throat_area, expansion_ratio):
    """Exit area from throat area and expansion ratio.

    Returns
    -------
    float – exit area [m^2]
    """
    return throat_area * expansion_ratio


def area_from_diameter(d):
    """Circle area from diameter [m]."""
    return math.pi / 4 * d**2


def diameter_from_area(a):
    """Diameter from circle area [m]."""
    return math.sqrt(4 * a / math.pi)


def _area_mach_relation(M, gamma):
    """(A/A*)^2 as a function of Mach number — used for root finding."""
    g = gamma
    t1 = 2 / (g + 1)
    t2 = 1 + (g - 1) / 2 * M**2
    exponent = (g + 1) / (g - 1)
    return (1 / M**2) * (t1 * t2) ** exponent


def expansion_ratio_from_pressure(gamma, pe_pc):
    """Expansion ratio ε = A_e/A_t for a given exit-to-chamber pressure ratio.

    Solves the isentropic area–Mach relation.

    Parameters
    ----------
    gamma : float – ratio of specific heats
    pe_pc : float – P_e / P_c (exit pressure / chamber pressure)

    Returns
    -------
    float – expansion ratio (A_e / A_t)
    """
    # First find exit Mach from pressure ratio
    Me = mach_from_pressure_ratio(gamma, pe_pc)
    return math.sqrt(_area_mach_relation(Me, gamma))


def mach_from_pressure_ratio(gamma, pe_pc):
    """Exit Mach number from pressure ratio P_e/P_c (isentropic)."""
    g = gamma
    return math.sqrt((2 / (g - 1)) * (pe_pc ** (-(g - 1) / g) - 1))


def exit_mach(gamma, expansion_ratio):
    """Supersonic Mach number for a given expansion ratio (iterative solve).

    Parameters
    ----------
    gamma : float
    expansion_ratio : float – A_e / A_t  (must be >= 1)

    Returns
    -------
    float – exit Mach number (supersonic solution)
    """
    eps_sq = expansion_ratio**2

    def f(M):
        return _area_mach_relation(M, gamma) - eps_sq

    # Supersonic branch: M > 1
    return brentq(f, 1.0001, 50.0)


def thrust_coefficient(gamma, expansion_ratio, pe_pc, pa_pc=0.0):
    """Thrust coefficient C_F.

    C_F = sqrt(2γ²/(γ-1) · (2/(γ+1))^((γ+1)/(γ-1)) · [1-(pe/pc)^((γ-1)/γ)])
          + ε·(pe/pc - pa/pc)

    Parameters
    ----------
    gamma : float
    expansion_ratio : float – ε = A_e / A_t
    pe_pc : float – exit pressure / chamber pressure
    pa_pc : float – ambient pressure / chamber pressure (0 for vacuum)

    Returns
    -------
    float – C_F (dimensionless)
    """
    g = gamma
    term1 = (2 * g**2) / (g - 1)
    term2 = (2 / (g + 1)) ** ((g + 1) / (g - 1))
    term3 = 1 - pe_pc ** ((g - 1) / g)
    cf_momentum = math.sqrt(term1 * term2 * term3)
    cf_pressure = expansion_ratio * (pe_pc - pa_pc)
    return cf_momentum + cf_pressure


def exit_velocity(gamma, T_c, mol_weight, expansion_ratio):
    """Exhaust velocity from isentropic expansion.

    Parameters
    ----------
    gamma : float
    T_c : float – chamber temperature [K]
    mol_weight : float – mean molecular weight of exhaust [kg/kmol]
    expansion_ratio : float

    Returns
    -------
    float – exit velocity [m/s]
    """
    g = gamma
    R = specific_gas_constant(mol_weight)
    Me = exit_mach(g, expansion_ratio)
    Te = T_c / (1 + (g - 1) / 2 * Me**2)
    return Me * math.sqrt(g * R * Te)


def nozzle_contour(throat_radius, expansion_ratio, half_angle_deg=15.0, n_points=100):
    """Simple conical nozzle contour (x, y) arrays for plotting.

    Parameters
    ----------
    throat_radius : float – [m]
    expansion_ratio : float
    half_angle_deg : float – cone half angle [degrees]
    n_points : int

    Returns
    -------
    x, y : numpy arrays [m] – axial and radial coordinates
    """
    r_t = throat_radius
    r_e = r_t * math.sqrt(expansion_ratio)
    alpha = math.radians(half_angle_deg)
    L = (r_e - r_t) / math.tan(alpha)

    x = np.linspace(0, L, n_points)
    y = r_t + x * math.tan(alpha)
    return x, y
