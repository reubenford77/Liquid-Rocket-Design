"""Propellant combination data.

Values are placeholders — update with real data or NASA CEA results.
All units are SI (metric).
"""

# Each entry: c_star (m/s), gamma, T_c (K), of_ratio, mol_weight (kg/kmol),
#             ox_density (kg/m^3), fuel_density (kg/m^3)
_PROPELLANTS = {
    "N2O/ETHANOL": {
        "name": "N2O / Ethanol",
        "c_star": 0,          # m/s  — TODO: fill in
        "gamma": 0,           # — TODO: fill in
        "T_c": 0,             # K    — TODO: fill in
        "of_ratio": 0,        # O/F by mass — TODO: fill in
        "mol_weight": 0,      # kg/kmol — TODO: fill in
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
