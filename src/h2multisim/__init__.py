"""
H2MultiSim™
===========
Open-source, dependency-free Python toolkit for multi-scale simulation of
underground hydrogen storage (UHS) across 0D, 1D, 2D, and 3D scales.

Author: Ziad Bennour (Subterranean Synergy / Curtin University)
License: GPL-3.0-or-later
"""

from h2multisim.models import (
    H2Reservoir0DModel,
    H2Reservoir1DModel,
    H2Reservoir2DModel,
    H2Reservoir3DModel,
    UHSAreal2D,
    UHSFieldScale3D,
    UHSRadial1D,
    UHSReservoir0D,
)
from h2multisim.thermo import (
    IsothermModel,
    PengRobinsonH2,
    PoreVolumeModel,
    adiabatic_temperature,
    coupled_temperature_update,
    cp_ideal_h2,
    gamma_real_h2,
    joule_thomson_coefficient,
    thermal_relaxation_update,
    viscosity_h2,
)
from h2multisim.units import (
    C_to_K,
    F_to_K,
    K_to_C,
    K_to_F,
    MPa_to_pa,
    bar_to_pa,
    kg_to_lb,
    kg_to_short_tons,
    kg_to_tonnes,
    lb_to_kg,
    m2_to_mD,
    m2_to_millidarcy,
    m3_to_cubic_feet,
    mD_to_m2,
    millidarcy_to_m2,
    pa_to_bar,
    pa_to_MPa,
    pa_to_psi,
    psi_to_pa,
    short_tons_to_kg,
    tonnes_to_kg,
)

__version__ = "1.0.0"
__all__ = [
    "__version__",
    # Main multi-scale simulators
    "UHSReservoir0D",
    "UHSRadial1D",
    "UHSAreal2D",
    "UHSFieldScale3D",
    "H2Reservoir0DModel",
    "H2Reservoir1DModel",
    "H2Reservoir2DModel",
    "H2Reservoir3DModel",
    # Thermodynamics
    "PengRobinsonH2",
    "PoreVolumeModel",
    "IsothermModel",
    "viscosity_h2",
    "gamma_real_h2",
    "joule_thomson_coefficient",
    "adiabatic_temperature",
    "thermal_relaxation_update",
    "coupled_temperature_update",
    "cp_ideal_h2",
    # Units
    "pa_to_bar",
    "bar_to_pa",
    "pa_to_MPa",
    "MPa_to_pa",
    "pa_to_psi",
    "psi_to_pa",
    "K_to_C",
    "C_to_K",
    "K_to_F",
    "F_to_K",
    "kg_to_tonnes",
    "tonnes_to_kg",
    "kg_to_short_tons",
    "short_tons_to_kg",
    "kg_to_lb",
    "lb_to_kg",
    "m3_to_cubic_feet",
    "millidarcy_to_m2",
    "m2_to_millidarcy",
    "mD_to_m2",
    "m2_to_mD",
]
