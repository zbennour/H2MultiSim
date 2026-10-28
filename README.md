# H2MultiSim™

**Multi-Scale Pure-Python Simulation Toolkit for Underground Hydrogen Storage (UHS)**

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
[![Python: 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/)
[![Lint: Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![Tests: 24/24 Passed](https://img.shields.io/badge/tests-24%2F24%20passed-brightgreen.svg)](tests/)

H2MultiSim™ is an open-source, dependency-free Python toolkit designed for researchers, subsurface reservoir engineers, and students to model underground hydrogen storage across four dimensional scales:

- **0D Tier**: Lumped material balance with strict BHP clamping and lognormal Monte Carlo uncertainty quantification (UQ).
- **1D Tier**: Radial finite-volume model on logarithmic grids with implicit Thomas algorithm solving near-wellbore transient drawdowns and buildups.
- **2D Tier**: Areal Cartesian grid with unconditionally stable Gauss-Seidel SOR ($\omega = 1.5$) modeling multi-well interference and permeability channeling.
- **3D Tier**: Field-scale Cartesian grid with red-black Gauss-Seidel SOR, layer stratigraphy, and buoyant gravity segregation ($+\rho_{\text{harm}} g \Delta z$).

---

## Key Technical Features

### 1. Robust Real-Gas Thermodynamics (`h2multisim.thermo`)
- **Closed-Form Analytical PR-EOS**: Solves the cubic Peng-Robinson equation of state using the Cardano/Viète trigonometric method. Guarantees the physical vapor root in a single pass without iterative Newton-Raphson failure or divergence.
- **Dynamic Viscosity**: Chapman-Enskog formulation calibrated against NIST and Hogeweg et al. (2022) benchmarks ($9.4\,\mu\text{Pa}\cdot\text{s}$ at $333\text{ K}, 82\text{ bar}$).
- **Coupled Thermal Model**: Real-gas isentropic adiabatic heating/cooling combined with exponential thermal relaxation toward geothermal equilibrium.
- **Thermodynamic Inversion**: Models the negative Joule-Thomson coefficient ($\mu_{JT} < 0$) across typical reservoir temperatures ($T > 202\text{ K}$).
- **Adsorption & Compressibility**: Dual-Site Langmuir (DSL) isotherms and rigid/elastic/breathing-MOF pore-volume models.

### 2. Multi-Scale Hierarchy

```
0D (Lumped)           1D (Radial)            2D (Areal)            3D (Field Scale)
┌─────────────────┐   ┌─────────────────┐    ┌─────────────────┐   ┌─────────────────┐
│ • Screening     │   │ • Near-wellbore │    │ • Multi-well    │   │ • Gravity plume │
│ • Monte Carlo   │──>│ • Wellbore skin │───>│ • Interference  │──>│ • Stratigraphy  │
│ • Thermodynamics│   │ • Transient P/T │    │ • Channel sweep │   │ • Buoyancy      │
└─────────────────┘   └─────────────────┘    └─────────────────┘   └─────────────────┘
```

---

## Benchmark Validation

H2MultiSim™ is validated against peer-reviewed experimental literature and NIST reference data:

| Metric | Condition | Reference Benchmark | H2MultiSim™ | Error |
| :--- | :--- | :--- | :--- | :--- |
| **Z-factor** | $293.15\text{ K}, 1\text{ bar}$ | $1.00060$ (NIST) | $1.00060$ | $< 0.01\%$ |
| **Z-factor** | $293.15\text{ K}, 100\text{ bar}$ | $1.05881$ (NIST) | $1.05814$ | $< 0.1\%$ |
| **Z-factor** | $293.15\text{ K}, 250\text{ bar}$ | $1.14283$ (NIST) | $1.14022$ | $< 0.25\%$ |
| **Gas Density** | $316\text{ K}, 80\text{ bar}$ | $5.98\text{ kg/m}^3$ (Okoroafor 2022) | $5.978\text{ kg/m}^3$ | $< 0.05\%$ |
| **Viscosity** | $293.15\text{ K}, 1\text{ atm}$ | $8.948\,\mu\text{Pa}\cdot\text{s}$ (NIST NTP) | $8.948\,\mu\text{Pa}\cdot\text{s}$ | $< 0.01\%$ |
| **Viscosity** | $333\text{ K}, 82\text{ bar}$ | $9.400\,\mu\text{Pa}\cdot\text{s}$ (Hogeweg 2022) | $9.385\,\mu\text{Pa}\cdot\text{s}$ | $< 0.2\%$ |
| **Thomas Solver**| $N=6$ Poisson system | Analytical solution | Residual $< 10^{-10}$ | Exact |

---

## Installation

H2MultiSim™ uses only the Python standard library for simulation and has zero external runtime dependencies.

```powershell
# Clone the repository
git clone https://github.com/zbennour/H2MultiSim.git
cd H2MultiSim

# Install package in editable mode (with uv or pip)
uv pip install -e .
# or: py -m pip install -e .
```

---

## Command-Line Interface (CLI)

Run simulations and export dual-unit CSV time series directly from PowerShell:

```powershell
# Run 0D Salt Cavern simulation (Teesside scenario)
py -m h2multisim run --scale 0d --scenario teesside

# Run 0D Depleted Reservoir with Monte Carlo UQ (500 realizations)
py -m h2multisim run --scale 0d --scenario okoroafor --mc --n-mc 500

# Run 1D Radial near-wellbore transient simulation
py -m h2multisim run --scale 1d --dt-days 0.1

# Run 2D Areal five-spot well interference simulation
py -m h2multisim run --scale 2d --scenario five_spot

# Run 2D Areal thief-zone channel simulation
py -m h2multisim run --scale 2d --scenario thief_zone

# Run 3D Field-scale gravity segregation simulation
py -m h2multisim run --scale 3d --scenario gravity_override

# Run complete 24-suite verification test suite
py -m h2multisim test
```

---

## Python API Usage

```python
from h2multisim import (
    UHSReservoir0D,
    build_teesside_cavern,
    schedule_cavern,
    ReservoirState,
    pa_to_bar,
)

# 1. Load canonical Teesside salt cavern benchmark
props = build_teesside_cavern()
schedule = schedule_cavern(props.reference_pressure)

# 2. Initialize simulator
model = UHSReservoir0D(props)
initial_mass = model.mass_from_state(props.reference_pressure, props.reference_temperature)
state = ReservoirState(
    pressure=props.reference_pressure,
    temperature=props.reference_temperature,
    mass=initial_mass,
)

# 3. Execute injection / withdrawal schedule
result = model.run_schedule(state, schedule, dt_s=86400.0)

print(f"Working Gas: {result.metrics['working_gas_kg'] / 1e3:.1f} tonnes")
print(f"Max Pressure: {pa_to_bar(result.metrics['max_pressure_Pa']):.1f} bar")
```

---

## Repository Structure

```
H2MultiSim/
├── pyproject.toml              # PEP 621 packaging and Ruff configuration
├── src/
│   └── h2multisim/
│       ├── __init__.py         # Public API and package version (1.0.0)
│       ├── __main__.py         # Module execution entrypoint
│       ├── cli.py              # Unified CLI (run, test)
│       ├── thermo.py           # PR-EOS, NIST viscosity, heat capacity, adsorption
│       ├── units.py            # Dual unit conversions & Cardano cubic solver
│       └── models/
│           ├── __init__.py     # Multi-scale exports
│           ├── zero_d.py       # 0D lumped material balance & Monte Carlo UQ
│           ├── one_d.py        # 1D radial finite-volume & Thomas solver
│           ├── two_d.py        # 2D areal SOR implicit multi-well model
│           └── three_d.py      # 3D field-scale model with gravity segregation
├── tests/
│   ├── __init__.py
│   └── test_validation.py     # 24 benchmark tests matching NIST & literature
└── README.md                   # Documentation and citation guide
```

---

## Citations & References

If you use H2MultiSim™ in academic research or industrial feasibility studies, please cite:

1. **Okoroafor, E. U., Saltzer, S., & Kovscek, A. R. (2022)**. *Assessing the potential for underground hydrogen storage in depleted gas fields*. International Journal of Hydrogen Energy, 47(78), 33781–33802.
2. **Hogeweg, S., Strobel, G., & Hagemann, B. (2022)**. *Benchmark study on numerical simulation of underground hydrogen storage in porous media*. Computational Geosciences, 26, 1367–1378.
3. **Peaceman, D. W. (1978)**. *Interpretation of well-block pressures in numerical reservoir simulation*. Society of Petroleum Engineers Journal, 18(3), 183–194.
4. **Peng, D.-Y., & Robinson, D. B. (1976)**. *A new two-constant equation of state*. Industrial & Engineering Chemistry Fundamentals, 15(1), 59–64.

---

## Authors & Acknowledgments

- **Ziad Bennour** ([@zbennour](https://github.com/zbennour)) — Architecture, multi-scale numerical engines (0D/1D/2D/3D), validation test suite, and package maintenance.
- **kwkushan** — Original author of the foundational V1.3 analytical Peng-Robinson EOS formulation, single-cell lumped model prototype, and benchmark scenario definitions.

---

## License

This project is licensed under the GNU General Public License v3.0 or later ([GPL-3.0-or-later](LICENSE)).
