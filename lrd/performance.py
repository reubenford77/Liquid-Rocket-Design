"""Top-level engine performance calculations."""

import math

from lrd.utils import G0
from lrd import propellants, nozzle, combustion, injector, tanks


def thrust(cf, chamber_pressure, throat_area):
    """Thrust force.  F = C_F · P_c · A_t

    Returns float – thrust [N]
    """
    return cf * chamber_pressure * throat_area


def specific_impulse(c_star, cf):
    """Specific impulse.  Isp = c* · C_F / g0

    Returns float – Isp [s]
    """
    return c_star * cf / G0


def mass_flow_rate(thrust_n, isp):
    """Total mass flow rate.  mdot = F / (Isp · g0)

    Returns float – mdot [kg/s]
    """
    return thrust_n / (isp * G0)


def mixture_ratio_flows(mdot_total, of_ratio):
    """Split total mass flow into oxidizer and fuel.

    Returns dict with mdot_ox and mdot_fuel [kg/s]
    """
    mdot_ox = mdot_total * of_ratio / (1 + of_ratio)
    mdot_fuel = mdot_total / (1 + of_ratio)
    return {"mdot_ox": mdot_ox, "mdot_fuel": mdot_fuel}


def delta_v(isp, mass_ratio):
    """Tsiolkovsky rocket equation.  Δv = Isp · g0 · ln(mass_ratio)

    Parameters
    ----------
    isp : float – specific impulse [s]
    mass_ratio : float – m_initial / m_final

    Returns float – Δv [m/s]
    """
    return isp * G0 * math.log(mass_ratio)


def engine_summary(
    thrust_target,
    chamber_pressure,
    propellant_key,
    ambient_pressure=101325.0,
    burn_time=30.0,
    contraction_ratio=3.0,
    injector_dp_fraction=0.20,
    cd_injector=0.7,
    orifice_diameter_ox=0.002,
    orifice_diameter_fuel=0.0015,
):
    """Compute a complete engine design point.

    Parameters
    ----------
    thrust_target : float – desired thrust [N]
    chamber_pressure : float – P_c [Pa]
    propellant_key : str – e.g. "N2O/ETHANOL"
    ambient_pressure : float – [Pa] (sea level default)
    burn_time : float – [s]
    contraction_ratio : float – A_c / A_t
    injector_dp_fraction : float – ΔP_inj / P_c
    cd_injector : float – discharge coefficient
    orifice_diameter_ox : float – [m]
    orifice_diameter_fuel : float – [m]

    Returns
    -------
    dict – comprehensive engine design summary
    """
    # Propellant properties
    prop = propellants.get_propellant(propellant_key)
    c_star = prop["c_star"]
    gamma = prop["gamma"]
    T_c = prop["T_c"]
    of_ratio = prop["of_ratio"]

    # Pressure ratios
    pe_pc = ambient_pressure / chamber_pressure
    pa_pc = ambient_pressure / chamber_pressure

    # Nozzle
    eps = nozzle.expansion_ratio_from_pressure(gamma, pe_pc)
    cf = nozzle.thrust_coefficient(gamma, eps, pe_pc, pa_pc)
    isp = specific_impulse(c_star, cf)

    # Flow rates
    At = thrust_target / (cf * chamber_pressure)
    mdot = mass_flow_rate(thrust_target, isp)
    flows = mixture_ratio_flows(mdot, of_ratio)

    # Nozzle geometry
    dt = nozzle.diameter_from_area(At)
    Ae = nozzle.exit_area(At, eps)
    de = nozzle.diameter_from_area(Ae)

    # Chamber
    l_star = combustion.l_star_typical(propellant_key)
    Vc = combustion.chamber_volume(At, l_star)
    chamber = combustion.chamber_dimensions(Vc, contraction_ratio, At)

    # Injector
    dp_inj = injector_dp_fraction * chamber_pressure
    n_ox = injector.orifice_count(flows["mdot_ox"], cd_injector, orifice_diameter_ox, dp_inj, prop["ox_density"])
    n_fuel = injector.orifice_count(flows["mdot_fuel"], cd_injector, orifice_diameter_fuel, dp_inj, prop["fuel_density"])

    # Tanks
    m_ox = flows["mdot_ox"] * burn_time
    m_fuel = flows["mdot_fuel"] * burn_time
    v_ox = tanks.tank_volume(m_ox, prop["ox_density"])
    v_fuel = tanks.tank_volume(m_fuel, prop["fuel_density"])

    feed_p = tanks.feed_pressure(chamber_pressure, dp_inj)

    return {
        "propellant": prop["name"],
        "thrust": thrust_target,
        "chamber_pressure": chamber_pressure,
        "c_star": c_star,
        "gamma": gamma,
        "T_c": T_c,
        "of_ratio": of_ratio,
        "expansion_ratio": eps,
        "thrust_coefficient": cf,
        "isp": isp,
        "mdot_total": mdot,
        "mdot_ox": flows["mdot_ox"],
        "mdot_fuel": flows["mdot_fuel"],
        "throat_area": At,
        "throat_diameter": dt,
        "exit_area": Ae,
        "exit_diameter": de,
        "l_star": l_star,
        "chamber_volume": Vc,
        "chamber_diameter": chamber["diameter"],
        "chamber_length": chamber["length"],
        "injector_dp": dp_inj,
        "injector_orifices_ox": n_ox,
        "injector_orifices_fuel": n_fuel,
        "ox_mass": m_ox,
        "fuel_mass": m_fuel,
        "ox_tank_volume": v_ox,
        "fuel_tank_volume": v_fuel,
        "feed_pressure": feed_p,
        "burn_time": burn_time,
    }
