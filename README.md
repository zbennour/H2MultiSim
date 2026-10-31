# H2MultiSim™

**Multi-Scale Pure-Python Simulation Toolkit for Underground Hydrogen Storage (UHS)**

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
[![Python: 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/)
[![Lint: Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![Tests: 86/86 Passed](https://img.shields.io/badge/tests-86%2F86%20passed-brightgreen.svg)](tests/)

H2MultiSim™ is an open-source, dependency-free Python toolkit designed for researchers, subsurface reservoir engineers, and students to model underground hydrogen storage across four dimensional scales:

- **0D Tier**: Lumped material balance with an exact mass-state step, active BHP/rate control, and Monte Carlo uncertainty quantification (UQ) with log-normal permeability and truncated-normal porosity, thickness, area and temperature.
- **1D Tier**: Radial finite-volume model on logarithmic grids with implicit Thomas algorithm solving near-wellbore transient drawdowns and buildups.
- **2D Tier**: Areal Cartesian grid with implicit finite volumes, a Jacobi-preconditioned conjugate-gradient pressure solver (`linear_solver="sor"` selects Gauss-Seidel SOR), Picard iteration, multi-well interference/permeability channeling, and optional binary H2 + cushion-gas (CH4/CO2/N2) mixing.
- **3D Tier**: Field-scale Cartesian grid with the same PCG pressure solver (red-black SOR optional), layer stratigraphy, buoyant gravity segregation ($+\rho_{\text{harm}} g \Delta z$), and optional binary H2 + cushion-gas mixing including H2 gravity override.

---

## Key Technical Features

### 1. Robust Real-Gas Thermodynamics (`h2multisim.thermo`)
- **Closed-Form Analytical PR-EOS & Peneloux Volume Translation**: Solves the cubic Peng-Robinson equation of state using the Cardano/Viète trigonometric method, guaranteeing the physical vapor root in a single pass without iterative Newton-Raphson divergence. Optional Peneloux volume translation ($c = -4.0 \times 10^{-6}\text{ m}^3/\text{mol}$) keeps the $Z$-factor error within $1.0\%$ across $1\text{–}350\text{ bar}$ at $280\text{–}450\text{ K}$ ($\le 1.25\%$ at $250\text{–}280\text{ K}$).
- **Dynamic Viscosity**: Chapman-Enskog dilute-gas formulation with a density-dependent residual fitted to Muzny et al. (2013) / CoolProp ($\mu = s_0\mu_0(T) + c_1\rho + c_2\rho^2$), matching CoolProp 8.0.0 within $0.5\%$ across $250\text{–}450\text{ K}$ and $1\text{–}350\text{ bar}$.
- **Thermal Model**: 0D tier solves a real-gas control-volume energy balance (Peng-Robinson enthalpy/internal energy, injected-gas enthalpy at `injection_temperature`, wall heat transfer via the relaxation time τ), matching CoolProp adiabatic cavern fill/discharge within 0.5–2 %. The 1D/2D/3D tiers use a per-cell real-gas isentropic update with exponential relaxation toward geothermal temperature (`thermal_model="lumped"` selects this in 0D).
- **Joule–Thomson**: Analytical real-gas throttling coefficient $\mu_{\text{JT}} = (\partial T/\partial P)_H$ with Peneloux covolume correction, reproducing negative $\mu_{\text{JT}}$ warming across the UHS operating envelope ($\sim 4.5\%$ vs CoolProp at $293.15\text{ K}, 100\text{ bar}$).
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

The thermodynamic core is a screening model.  The values below were measured
against CoolProp 8.0.0 (Leachman et al. normal-hydrogen reference EOS); they are
reproducible with `audit/bench_thermo.py`.  Do not use the cubic EOS for custody
transfer or regulatory reporting.

> [!TIP]
> For the complete quantitative benchmark report across all three tiers (**Tier 1: NIST/CoolProp 8.0.0**, **Tier 2: Al-Hussainy (1966) analytical pseudopressure**, and **Tier 3: OPM Flow & ECLIPSE cross-simulator decks**), see [BENCHMARKS.md](BENCHMARKS.md).
>
> What each tier does **not** model, and the validated accuracy envelope, are listed in [LIMITATIONS.md](LIMITATIONS.md).

| Metric | Condition | Reference (CoolProp) | H2MultiSim™ | Error |
| :--- | :--- | :--- | :--- | :--- |
| **Z-factor** (plain PR) | $293.15\text{ K}, 1\text{ bar}$ | $1.00059$ | $1.00027$ | $0.03\%$ |
| **Z-factor** (plain PR) | $293.15\text{ K}, 100\text{ bar}$ | $1.06081$ | $1.03537$ | $2.40\%$ |
| **Z-factor** (plain PR) | $293.15\text{ K}, 250\text{ bar}$ | $1.15751$ | $1.10938$ | $4.16\%$ |
| **Z-factor** (Peneloux) | $293.15\text{ K}, 250\text{ bar}$ | $1.15751$ | $1.15040$ | $0.61\%$ |
| **Gas Density** (plain PR)| $316\text{ K}, 80\text{ bar}$ | $5.868\text{ kg/m}^3$ | $5.978\text{ kg/m}^3$ | $1.86\%$ |
| **Viscosity** (NIST/Muzny)| $293.15\text{ K}, 1\text{ atm}$ | $8.797\,\mu\text{Pa}\cdot\text{s}$ | $8.746\,\mu\text{Pa}\cdot\text{s}$ | $0.58\%$ |
| **Viscosity** (NIST/Muzny)| $333\text{ K}, 82\text{ bar}$ | $9.676\,\mu\text{Pa}\cdot\text{s}$ | $9.603\,\mu\text{Pa}\cdot\text{s}$ | $0.75\%$ |
| **Viscosity** (NIST/Muzny)| $293.15\text{ K}, 250\text{ bar}$| $9.302\,\mu\text{Pa}\cdot\text{s}$ | $9.186\,\mu\text{Pa}\cdot\text{s}$ | $1.25\%$ |
| **Joule-Thomson** (Peneloux)| $293.15\text{ K}, 100\text{ bar}$| $-0.0347\text{ K/bar}$ | $-0.0331\text{ K/bar}$ | $4.5\%$ |
| **Thomas Solver**| $N=6$ Poisson system | Analytical solution | Residual $< 10^{-10}$ | Exact |

The uncorrected Peng–Robinson Z-factor error reaches ~4.2% at 250 bar. Enabling Peneloux volume
translation ($c=-4.0\times10^{-6}\ \mathrm{m^3/mol}$) keeps the maximum Z error within 1.0 % across
1–350 bar at 280–450 K (1.25 % at 250–280 K) and is supported across all models and JSON/YAML configs via `peneloux=True` or
`volume_translation=-4.0e-6` (see `AUDIT.md`, finding F-04). Dynamic viscosity uses a dense-fluid
residual fitted to Muzny et al. (2013) / CoolProp, within 0.5 % across 250–450 K and 1–350 bar.
`h2multisim verify` re-checks all of these figures against embedded CoolProp reference values.
Mass balance in the 1D/2D/3D solvers (including constant-pressure boundary influx in 1D) is closed to better than $10^{-6}$ relative error. The 0D tier reports step-average rates, so its mass balance holds by construction; solver health is reported by the `converged` / `unconverged_steps` metrics in every tier.

---

## Configuration-driven runs

Instead of the built-in scenarios, a full run can be described by a JSON (or YAML,
if PyYAML is installed) file:

```bash
h2multisim run --config examples/04_config_driven_2d.json
```

The file specifies `scale`, `reservoir`, `wells`, `schedule`, `dt_days` and `out_dir`;
every field maps onto the model dataclasses (see `src/h2multisim/config.py`). Pressures
accept either pascals or strings such as `"80 bar"`; permeability can be given in mD
via `permeability_mD`. This removes the need to edit Python source to study a new
reservoir.

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

# Export standard ECLIPSE 100/300 & OPM Flow (.DATA) input decks
py -m h2multisim export-deck --scale 1d --scenario okoroafor
py -m h2multisim export-deck --scale 2d --scenario five_spot
py -m h2multisim export-deck --scale 3d --scenario layered_reservoir
py -m h2multisim export-deck --config simulations/custom_reservoir.json

# Run complete 40-test verification test suite
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
├── pyproject.toml              # PEP 621 packaging, Pytest, and Ruff configuration
├── BENCHMARKS.md               # 3-Tier Industry Benchmark Report (NIST, Al-Hussainy, OPM Flow)
├── NUMERICAL_TOLERANCES.md     # Numerical stability, tolerances, and solver convergence limits
├── audit/
│   ├── bench_thermo.py         # CoolProp 8.0.0 vs PR-EOS thermodynamic benchmark
│   └── bench_analytical_1d.py  # Al-Hussainy real-gas pseudopressure radial line-source benchmark
├── src/
│   └── h2multisim/
│       ├── __init__.py         # Public API and package version (1.9.0)
│       ├── __main__.py         # Module execution entrypoint
│       ├── cli.py              # Unified CLI (run, export-deck, test)
│       ├── config.py           # JSON/YAML run configuration loader
│       ├── deck_export.py      # ECLIPSE 100/300 & OPM Flow (.DATA) input deck generator
│       ├── thermo.py           # PR-EOS, Peneloux, NIST/Muzny viscosity, heat capacity, adsorption
│       ├── units.py            # Dual unit conversions & Cardano cubic solver
│       └── models/
│           ├── __init__.py     # Multi-scale exports
│           ├── zero_d.py       # 0D lumped material balance & Monte Carlo UQ
│           ├── one_d.py        # 1D radial finite-volume & Thomas solver
│           ├── two_d.py        # 2D areal implicit multi-well model (PCG/SOR)
│           └── three_d.py      # 3D field-scale model with gravity segregation
├── tests/
│   ├── __init__.py
│   ├── test_char_2d.py         # 2D single-well golden trajectory characterization test
│   └── test_validation.py      # 39 benchmark tests matching NIST, analytical solutions & literature
└── README.md                   # Documentation and citation guide
```

---

## Citations & References

If you use H2MultiSim™ in academic research or industrial feasibility studies, please cite:

1. **Okoroafor, E. U., Saltzer, S., & Kovscek, A. R. (2022)**. *Assessing the potential for underground hydrogen storage in depleted gas fields*. International Journal of Hydrogen Energy, 47(78), 33781–33802.
2. **Hogeweg, S., Strobel, G., & Hagemann, B. (2022)**. *Benchmark study on numerical simulation of underground hydrogen storage in porous media*. Computational Geosciences, 26, 1367–1378.
3. **Al-Hussainy, R., Ramey, H. J., & Crawford, P. B. (1966)**. *The Flow of Real Gases Through Porous Media*. Journal of Petroleum Technology, 18(05), 624–636.
4. **Peaceman, D. W. (1978)**. *Interpretation of well-block pressures in numerical reservoir simulation*. Society of Petroleum Engineers Journal, 18(3), 183–194.
5. **Peng, D.-Y., & Robinson, D. B. (1976)**. *A new two-constant equation of state*. Industrial & Engineering Chemistry Fundamentals, 15(1), 59–64.

---

## Authors & Acknowledgments

- **Ziad Bennour** ([@zbennour](https://github.com/zbennour)); Architecture, multi-scale numerical engines (0D/1D/2D/3D), validation test suite, and package maintenance.
- **Kushan Sandunil** ([@kwkushan](https://github.com/kwkushan)); Original author of the foundational V1.3 analytical Peng-Robinson EOS formulation, single-cell lumped model prototype, and benchmark scenario definitions.

---

## License

This project is licensed under the GNU General Public License v3.0 or later ([GPL-3.0-or-later](LICENSE)).
