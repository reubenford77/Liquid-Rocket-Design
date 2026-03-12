"""Propellant combination data.

Values are placeholders — update with real data or NASA CEA results.
All units are SI (metric).
"""

# Each entry: c_star (m/s), gamma, T_c (K), of_ratio, mol_weight (kg/kmol),
#             ox_density (kg/m^3), fuel_density (kg/m^3)
_PROPELLANTS = {
    "N2O/ETHANOL": {
        "name": "N2O / Ethanol",
        "c_star": 1419.0,          # m/s  
        "gamma": 1.198,          
        "T_c": 2662.6,             # K   
        "of_ratio": 4.5,        # O/F by mass
        "mol_weight": 26.01,      # kg/kmol
        "ox_density": 1220,   # kg/m^3 (liquid N2O at ~20°C, self-pressurising)
        "fuel_density": 789,  # kg/m^3 (ethanol)
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
        Key like "N2O/ETHANOL".

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
