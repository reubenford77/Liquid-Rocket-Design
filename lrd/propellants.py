"""Propellant combination data from standard references (Sutton, Humble).

Values are representative at typical chamber pressures (~2-7 MPa).
For higher fidelity, use NASA CEA.
"""

# Each entry: c_star (m/s), gamma, T_c (K), of_ratio, mol_weight (kg/kmol)
_PROPELLANTS = {
    "LOX/RP-1": {
        "name": "LOX / RP-1 (Kerosene)",
        "c_star": 1774,       # m/s
        "gamma": 1.24,
        "T_c": 3571,          # K
        "of_ratio": 2.56,     # O/F by mass
        "mol_weight": 23.3,   # kg/kmol
        "ox_density": 1141,   # kg/m^3 (LOX)
        "fuel_density": 820,  # kg/m^3 (RP-1)
    },
    "LOX/LCH4": {
        "name": "LOX / Liquid Methane",
        "c_star": 1835,
        "gamma": 1.20,
        "T_c": 3526,
        "of_ratio": 3.20,
        "mol_weight": 21.4,
        "ox_density": 1141,
        "fuel_density": 422,
    },
    "LOX/LH2": {
        "name": "LOX / Liquid Hydrogen",
        "c_star": 2386,
        "gamma": 1.26,
        "T_c": 3250,
        "of_ratio": 5.00,
        "mol_weight": 10.0,
        "ox_density": 1141,
        "fuel_density": 71,
    },
    "N2O4/MMH": {
        "name": "NTO / Monomethylhydrazine",
        "c_star": 1724,
        "gamma": 1.25,
        "T_c": 3200,
        "of_ratio": 2.15,
        "mol_weight": 22.0,
        "ox_density": 1440,
        "fuel_density": 878,
    },
}


def list_propellants():
    """Return list of available propellant combination keys."""
    return list(_PROPELLANTS.keys())


def get_propellant(name):
    """Look up propellant properties by name.

    Parameters
    ----------
    name : str
        Key like "LOX/RP-1", "LOX/LCH4", "LOX/LH2", "N2O4/MMH".

    Returns
    -------
    dict with keys: name, c_star, gamma, T_c, of_ratio, mol_weight,
                    ox_density, fuel_density.
    """
    key = name.upper().replace(" ", "")
    if key not in _PROPELLANTS:
        available = ", ".join(_PROPELLANTS.keys())
        raise ValueError(f"Unknown propellant '{name}'. Available: {available}")
    return dict(_PROPELLANTS[key])  # return a copy
