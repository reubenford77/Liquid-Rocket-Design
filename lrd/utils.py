"""Shared constants and helpers. All units are SI (metric)."""

# --- Physical constants ---
G0 = 9.80665          # m/s^2, standard gravity
R_UNIVERSAL = 8314.46  # J/(kmol·K)

# --- Gas helpers ---

def specific_gas_constant(molecular_weight):
    """R_specific = R_universal / M  [J/(kg·K)]."""
    return R_UNIVERSAL / molecular_weight
