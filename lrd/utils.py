"""Shared constants and unit conversion helpers."""

import math

# --- Physical constants ---
G0 = 9.80665          # m/s^2, standard gravity
R_UNIVERSAL = 8314.46  # J/(kmol·K)

# --- Unit conversions ---

def psi_to_pa(psi):
    """Pounds per square inch → Pascals."""
    return psi * 6894.757

def pa_to_psi(pa):
    """Pascals → pounds per square inch."""
    return pa / 6894.757

def lb_to_n(lb):
    """Pounds-force → Newtons."""
    return lb * 4.44822

def n_to_lb(n):
    """Newtons → pounds-force."""
    return n / 4.44822

def in_to_m(inches):
    """Inches → meters."""
    return inches * 0.0254

def m_to_in(meters):
    """Meters → inches."""
    return meters / 0.0254

def lbm_to_kg(lbm):
    """Pounds-mass → kilograms."""
    return lbm * 0.453592

def kg_to_lbm(kg):
    """Kilograms → pounds-mass."""
    return kg / 0.453592

# --- Gas helpers ---

def specific_gas_constant(molecular_weight):
    """R_specific = R_universal / M  [J/(kg·K)]."""
    return R_UNIVERSAL / molecular_weight
