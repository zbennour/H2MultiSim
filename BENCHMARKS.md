# H2MultiSim™ Industry Benchmarking & Verification Suite

Comprehensive multi-tier benchmark validation of **H2MultiSim™** against peer-reviewed analytical reservoir engineering solutions, NIST reference thermodynamic models (CoolProp 8.0.0), and industry-standard numerical reservoir simulators (Open Porous Media / OPM Flow and Schlumberger ECLIPSE 100/300).

---

## 1. Benchmarking Architecture Overview

Underground Hydrogen Storage (UHS) poses distinct numerical and thermodynamic challenges that differ significantly from conventional natural gas ($CH_4$) and carbon capture ($CO_2$) storage:
1. **High Non-Ideality & Covolume Repulsion**: Pure $H_2$ exhibits $Z > 1.0$ across standard operating pressures ($50\text{–}250\text{ bar}$), requiring accurate covolume formulation.
2. **Low Molecular Weight & Viscosity**: With $M = 2.016\text{ g/mol}$ and dynamic viscosity $\mu \sim 0.009\text{–}0.015\text{ cP}$ ($\approx 1/10\text{th}$ of methane), $H_2$ possesses high mobility and diffusivity, driving steep near-wellbore pressure gradients.
3. **Severe Density Contrast & Gravity Segregation**: With $\rho \approx 5\text{–}20\text{ kg/m}^3$ at reservoir conditions, buoyancy forces drive rapid vertical migration and gravity override in porous aquifers and depleted reservoirs.

To guarantee industrial rigor while maintaining **zero external runtime dependencies**, H2MultiSim employs a three-tier benchmarking architecture:

```
┌────────────────────────────────────────────────────────────────────────┐
│                   H2MultiSim™ Verification Framework                   │
└────────────────────────────────────────────────────────────────────────┘
                                    │
    ┌───────────────────────────────┼───────────────────────────────┐
    ▼                               ▼                               ▼
┌──────────────────────┐ ┌──────────────────────┐ ┌──────────────────────┐
│  Tier 1: Thermo/PVT  │ │  Tier 2: Analytical  │ │  Tier 3: Simulators  │
│  NIST WebBook /      │ │  Al-Hussainy et al.  │ │  OPM Flow / ECLIPSE  │
│  CoolProp 8.0.0      │ │  Pseudopressure m(p) │ │  Input Deck Export   │
│  (Leachman 2009 EOS) │ │  Line-Source PDE Sol │ │  (.DATA Dry Gas)     │
└──────────────────────┘ └──────────────────────┘ └──────────────────────┘
```

---

## 2. Tier 1: Thermodynamic & Transport Verification (vs NIST / CoolProp 8.0.0)

H2MultiSim's analytical Peng-Robinson cubic equation of state (`PengRobinsonH2`) with Peneloux volume translation and a dense-fluid viscosity residual fitted to Muzny et al. (2013) was verified against **CoolProp 8.0.0**, which evaluates the international reference standard Helmholtz Energy Equation of State (HEOS) of Leachman et al. (2009) for normal hydrogen.

### Benchmark Results (Reproducible via `audit/bench_thermo.py`)

| Quantity | Temperature | Pressure | CoolProp 8.0.0 (Leachman) | H2MultiSim (PR + Peneloux) | Relative Error |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Z-factor** | $293.15\text{ K}$ | $1.0\text{ bar}$ | $1.00059$ | $1.00044$ | $0.02\%$ |
| **Z-factor** | $293.15\text{ K}$ | $50.0\text{ bar}$ | $1.02999$ | $1.02402$ | $0.58\%$ |
| **Z-factor** | $293.15\text{ K}$ | $100.0\text{ bar}$ | $1.06081$ | $1.05178$ | $0.85\%$ |
| **Z-factor** | $293.15\text{ K}$ | $200.0\text{ bar}$ | $1.12473$ | $1.11554$ | $0.82\%$ |
| **Z-factor** | $316.00\text{ K}$ | $80.0\text{ bar}$ | $1.04594$ | $1.03907$ | $0.66\%$ |
| **Z-factor** | $373.15\text{ K}$ | $100.0\text{ bar}$ | $1.05072$ | $1.04574$ | $0.47\%$ |
| **Density** | $293.15\text{ K}$ | $100.0\text{ bar}$ | $7.797\text{ kg/m}^3$ | $7.864\text{ kg/m}^3$ | $0.86\%$ |
| **Density** | $316.00\text{ K}$ | $80.0\text{ bar}$ | $5.868\text{ kg/m}^3$ | $5.907\text{ kg/m}^3$ | $0.66\%$ |
| **Viscosity $\mu$** | $293.15\text{ K}$ | $1.0\text{ bar}$ | $8.797\ \mu\text{Pa·s}$ | $8.824\ \mu\text{Pa·s}$ | $0.31\%$ |
| **Viscosity $\mu$** | $293.15\text{ K}$ | $100.0\text{ bar}$ | $8.922\ \mu\text{Pa·s}$ | $8.932\ \mu\text{Pa·s}$ | $0.10\%$ |
| **Viscosity $\mu$** | $373.15\text{ K}$ | $100.0\text{ bar}$ | $10.464\ \mu\text{Pa·s}$ | $10.474\ \mu\text{Pa·s}$ | $0.10\%$ |
| **Viscosity $\mu$** | $293.15\text{ K}$ | $300.0\text{ bar}$ | $9.462\ \mu\text{Pa·s}$ | $9.476\ \mu\text{Pa·s}$ | $0.15\%$ |
| **Heat capacity $c_p$** | $293.15\text{ K}$ | $100.0\text{ bar}$ | $14.535\text{ kJ/(kg·K)}$ | $14.567\text{ kJ/(kg·K)}$ | $0.22\%$ |
| **Joule–Thomson $\mu_{\text{JT}}$** | $293.15\text{ K}$ | $100.0\text{ bar}$ | $-0.0347\text{ K/bar}$ | $-0.0325\text{ K/bar}$ | $6.28\%$ |

Regenerated for v1.8.0. `h2multisim verify` re-checks 64 such points against embedded CoolProp references and exits non-zero if any exceeds its tolerance.

> [!NOTE]
> Across 280–450 K and 1–350 bar, H2MultiSim's Peneloux-translated Peng-Robinson EOS is within **1.0 % on $Z$-factor and density** (1.25 % at 250–280 K) and **0.5 % on dynamic viscosity** compared to the reference Leachman (2009) EOS.

---

## 3. Tier 2: Analytical Pseudopressure Verification (Al-Hussainy, Ramey & Crawford 1966)

To rigorously verify the numerical spatial discretisation, logarithmic grid refinement, and implicit Picard nonlinear solvers without simulator cross-dependency, H2MultiSim is benchmarked against the exact closed-form analytical line-source solution for transient compressible gas flow in porous media.

### Mathematical Formulation

The non-linear partial differential equation governing compressible radial gas flow is:
$$\frac{1}{r} \frac{\partial}{\partial r} \left( r \frac{p}{\mu(p) Z(p)} \frac{\partial p}{\partial r} \right) = \frac{\phi c_t}{k} \frac{p}{Z(p)} \frac{\partial p}{\partial t}$$

Following the classic Al-Hussainy, Ramey & Crawford (1966) transformation, the equation is linearized through the real-gas pseudopressure integral:
$$m(p) = 2 \int_0^p \frac{p}{\mu(p) Z(p)} dp$$

For a well producing or injecting at constant mass rate $q_{\text{mass}}$ from an infinite-acting homogeneous reservoir with initial pressure $p_i$, the exact analytical solution is:
$$m(r, t) = m(p_i) - \frac{q_{\text{mass}} R T}{\pi k h M} E_1\left( \frac{\phi \mu_i c_{ti} r^2}{4 k t} \right)$$
where:
- $E_1(x) = \int_x^\infty \frac{e^{-u}}{u} du$ is the exponential integral.
- $m(p)$ and its inverse $m^{-1}(m) \to p$ are computed with high-order adaptive Simpson quadrature coupled to H2MultiSim's analytical EOS and viscosity model.
- $E_1(x)$ is evaluated via a machine-precision pure-Python series expansion for $x < 1$ and a continuous fraction expansion for $x \ge 1$ (relative accuracy $< 10^{-14}$ vs `scipy.special.exp1`).

### Benchmark Scenario Specification
- **Initial Reservoir Pressure ($p_i$)**: $80.0\text{ bar}$ ($8.0\times 10^6\text{ Pa}$)
- **Reservoir Temperature ($T$)**: $316.0\text{ K}$ ($42.85^\circ\text{C}$)
- **Permeability ($k$)**: $100\text{ mD}$ ($9.869\times 10^{-14}\text{ m}^2$)
- **Net Pay Thickness ($h$)**: $50.0\text{ m}$
- **Porosity ($\phi$)**: $0.20$
- **Total Compressibility ($c_t$)**: $1.25\times 10^{-9}\text{ Pa}^{-1}$
- **Wellbore Radius ($r_w$)**: $0.10\text{ m}$
- **Reservoir Outer Radius ($r_e$)**: $3000.0\text{ m}$
- **Withdrawal Mass Rate ($q_{\text{mass}}$)**: $5.0\text{ kg/s}$ ($\approx 432\text{ tonnes/day}$)
- **Simulation Duration**: $5.0\text{ days}$

### Benchmark Verification Results

Executing `py audit/bench_analytical_1d.py` or running test `test_38_analytical_real_gas_pseudopressure`:

#### 1. Flowing Wellbore BHP vs Time
| Time (days) | Analytical $p_{\text{wf}}$ (bar) | H2MultiSim $p_{\text{wf}}$ (bar) | Absolute Error (bar) | Relative Error (%) |
| :--- | :--- | :--- | :--- | :--- |
| **0.50** | $57.942$ | $58.067$ | $+0.125$ | **0.216%** |
| **1.00** | $57.519$ | $57.579$ | $+0.060$ | **0.104%** |
| **1.50** | $57.269$ | $57.294$ | $+0.025$ | **0.043%** |
| **2.00** | $57.091$ | $57.093$ | $+0.001$ | **0.002%** |
| **3.00** | $56.840$ | $56.811$ | $-0.029$ | **0.052%** |
| **4.00** | $56.662$ | $56.611$ | $-0.051$ | **0.090%** |
| **5.00** | $56.524$ | $56.519$ | $-0.005$ | **0.010%** |

#### 2. Spatial Radial Pressure Profile at $t = 5.0\text{ days}$
| Radius $r$ (m) | Analytical $P(r)$ (bar) | H2MultiSim $P(r)$ (bar) | Error (bar) | Relative Error (%) |
| :--- | :--- | :--- | :--- | :--- |
| **0.114** (cell 1) | $56.861$ | $56.877$ | $+0.016$ | **0.028%** |
| **0.297** | $59.570$ | $59.574$ | $+0.004$ | **0.007%** |
| **1.018** | $62.827$ | $62.836$ | $+0.009$ | **0.014%** |
| **5.077** | $67.060$ | $67.078$ | $+0.018$ | **0.027%** |
| **25.32** | $71.603$ | $71.637$ | $+0.034$ | **0.048%** |
| **95.66** | $74.908$ | $74.964$ | $+0.056$ | **0.075%** |
| **212.8** | $76.845$ | $76.918$ | $+0.073$ | **0.095%** |
| **473.4** | $78.841$ | $78.926$ | $+0.085$ | **0.108%** |
| **1561.4** | $79.999$ | $80.000$ | $+0.001$ | **0.001%** |

- **Max Flowing BHP Error**: $\mathbf{0.228\%}$ at $t=0.5\text{ d}$, improving to $\mathbf{0.004\%}$ at $t=5.0\text{ d}$ (v1.2.0: 0.216 %).
- **Max Spatial Profile Error**: $\mathbf{0.121\%}$ across 4 orders of magnitude in radius ($0.11\text{–}1561\text{ m}$) (v1.2.0: 0.108 %).
- **Global Mass Balance Closure Error**: $\mathbf{1.66\times 10^{-13}}$.

> The tables above were captured with v1.2.0. From v1.3.1 the 1D solver applies the rock compressibility $C_t$ that the analytical solution already included in $c_t$, and uses the real-gas isentropic temperature update; the summary figures above are from `audit/bench_analytical_1d.py` on v1.3.1.

---

## 4. Tier 3: Numerical Cross-Simulator Deck Bridge (OPM Flow & ECLIPSE 100/300)

H2MultiSim provides a native, built-in deck export engine (`h2multisim.deck_export`) that converts any 1D, 2D, or 3D simulation model or scenario into fully valid ECLIPSE 100/300 and OPM Flow `.DATA` format.

### Model equivalence (v1.3.2)

OPM Flow has no single-phase gas mode (decks declaring only `GAS` stop with *"No suitable configuration found"*), so the decks are **gas–water** with immobile connate water:

- `SWAT`/`SGAS` (or `EQUIL` in 3D) place water at the solver's $S_{wr}$; `SGWFN` tabulates the solver's own `RelativePermeabilityModel` $k_{rg}(S_g)$, $k_{rw}(S_w)$ with zero capillary pressure; `PORO`/`PERM*` are the unscaled rock values.
- `ROCK` carries the pore compressibility $C_t$ the 1D/2D/3D solvers apply.
- 3D layers and perforations are written top-down (ECLIPSE `K = 1` is the top; H2MultiSim layer 0 is the bottom), with `PERMZ = PERMX` (the solver is isotropic).
- The 1D radial model is written as an $N\times1\times1$ strip whose `PORV` and `TRANX` (EDIT section) are the exact log-grid pore volumes and radial transmissibilities, with an explicit `COMPDAT` connection factor using the solver's equivalent radius. OPM builds a cylindrical grid from `RADIAL` but does not connect the well to it.
- Producer rate caps go to the `WCONPROD` GRAT item; a `SUMMARY` section requests `FPR`, `FGIP`, `FGIT`, `FGPT`, `WBHP`, `WGIR`, `WGPR`.
- The decks are isothermal; H2MultiSim's adiabatic/relaxation thermal model has no deck counterpart.

### Cross-simulator results (OPM Flow 2026.04, H2MultiSim v1.10.0, reproducible via `audit/bench_opm.py`)

Each canonical scenario is exported, run in OPM Flow with its time step capped at H2MultiSim's $\Delta t$, and compared at every OPM report time. "dGIP" is the error in the change of gas in place (OPM: cumulative injection − production), normalised by the largest inventory swing; "P" is the error in pore-volume-average pressure (`FPR`, uniform Cartesian grids only).

| Scenario | H2MultiSim isothermal: max dGIP err | max P err | H2MultiSim thermal: max dGIP err | max P err |
| :--- | ---: | ---: | ---: | ---: |
| 1D Okoroafor (radial, 40 cells) | **0.40 %** | – | 0.48 % | – |
| 2D single well (15×15) | **0.71 %** | **0.010 %** | 0.74 % | 0.009 % |
| 2D five-spot (15×15) | **1.18 %** | **0.052 %** | 1.25 % | 0.048 % |
| 3D gravity override (8×8×5) | **0.19 %** | **0.001 %** | 0.21 % | 0.002 % |
| 3D layered reservoir (8×8×5) | **0.49 %** | **0.032 %** | 0.54 % | 0.027 % |

With matched physics (isothermal) the two simulators agree to ≤ 1.2 % in stored-gas change and ≤ 0.05 % in average pressure. The "thermal" columns use H2MultiSim's default thermal model. Since v1.10.0 that model includes the rock and connate-water heat capacity, so porous reservoirs stay within ~1 K of geothermal and also agree with the isothermal decks to ≤ 1.25 %. Before v1.10.0 the gas-only thermal model heated porous reservoirs by tens of kelvin and differed by 3–13 %, which was a modelling error rather than real physics. With OPM's default adaptive stepping (steps up to tens of days), the isothermal gap grows to 1–9 %, which is OPM's own time-discretisation error on these fast transients.

### Exporting Input Decks via CLI

Use the `export-deck` subcommand:

```powershell
# 1D Radial Wellbore Model (exported as an exact PORV/TRANX strip)
py -m h2multisim export-deck --scale 1d --scenario okoroafor --out-file output/1d_okoroafor.DATA

# 2D Multi-Well Areal Model (Five-Spot Pattern)
py -m h2multisim export-deck --scale 2d --scenario five_spot --out-file output/2d_five_spot.DATA

# 3D Field-Scale Stratigraphic Reservoir
py -m h2multisim export-deck --scale 3d --scenario layered_reservoir --out-file output/3d_layered.DATA

# Custom JSON/YAML Model Configuration
py -m h2multisim export-deck --config simulations/custom_reservoir.json
```

### Running the Exported Decks in OPM Flow

OPM Flow is the premier open-source reservoir simulator maintained by the Open Porous Media initiative. To execute an exported H2MultiSim deck in OPM Flow:

#### Native Linux / WSL2:
```bash
sudo apt install opm-simulators
flow output/1d_okoroafor.DATA
```

#### Docker:
```bash
docker run --rm -v ${PWD}:/models openporousmedia/opm:latest flow /models/output/1d_okoroafor.DATA
```

The resulting `.EGRID`, `.UNRST`, and `.UNSMRY` output files can be visualised directly in ResInsight or compared against H2MultiSim's exported `.csv` timeseries.

---

## 5. Automated Verification Test Suite

All benchmark criteria are embedded directly into H2MultiSim's automated test suite:

| Test ID | Method | Verified Standard | Benchmark Criteria |
| :--- | :--- | :--- | :--- |
| `test_01` | `test_01_eos_zfactor_vs_nist` | NIST WebBook | Plain PR $Z(P)$ error $< 1.0\%$ |
| `test_02` | `test_02_eos_density_okoroafor` | Okoroafor (2022) | Density error $< 1.5\%$ at $40\text{–}100\text{ bar}$ |
| `test_05` | `test_05_viscosity_benchmarks` | Muzny / NIST (2013) | Dense gas viscosity error $< 2.0\%$ |
| `test_07` | `test_07_joule_thomson` | CoolProp 8.0.0 | Negative $\mu_{\text{JT}}$ within $< 5.0\%$ |
| `test_32` | `test_32_peneloux_translation` | CoolProp 8.0.0 | Peneloux $Z(P)$ error $\le 1.0\%$ ($1\text{–}350\text{ bar}$, 280–450 K) |
| `test_38` | `test_38_analytical_real_gas_pseudopressure` | Al-Hussainy et al. (1966) | Flowing BHP error $< 0.50\%$, spatial profile $< 0.25\%$ |
| `test_39` | `test_39_eclipse_deck_export` | OPM Flow / ECLIPSE | Monotonic $B_g$, valid keyword structures across 1D/2D/3D |

To execute all 40 automated tests:
```powershell
uv run pytest -v
```
