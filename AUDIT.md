# H2MultiSim™ v1.0.0 — Independent Technical Audit & Industrial Benchmark

**Repository:** `zbennour/H2MultiSim` @ `7b10e35` (main)
**Audit date:** 2026-10-28
**Environment:** Python 3.12.3 · CoolProp 8.0.0 (HEOS / Leachman 2009 reference EOS) · NumPy
**Scope:** all 10 package modules, 24 tests (31 after the v1.0.1 patch), CLI, packaging,
documentation; 0D/1D/2D/3D flow models; thermodynamics.

---

## 0. Verdict

> **Not ready for quantitative use.** The architecture, packaging and thermodynamics are a solid
> screening-tool foundation, but the 1D/2D/3D flow solvers are **not mass-conservative**, and the
> reported throughput is disconnected from the mass that actually enters the model by **12× to 889×**
> (measured). The flagship 0D Teesside case reports **3.11 Mt injected and zero net cycles** while the
> gas inventory actually **falls by 58 %**. The “24/24 validation tests” pass because most are
> sign/bound checks, not validation. The README benchmark table does not agree with NIST/HEOS.

| Area | Grade | Comment |
|---|---|---|
| Architecture / packaging | **B+** | Clean module split, zero-dep wheel, typed, ruff-clean |
| Thermodynamics | **C+** | Correct structure; 2–5 % EOS error, pressure-independent μ, wrong J-T |
| 0D model | **C−** | Physically sensible, but clamp + reporting breaks mass balance |
| 1D/2D/3D flow | **F** | Wrong accumulation term; 1 213× too stiff; no mass balance |
| Well model | **F** | Cell-pressure clamp instead of Peaceman; rate cap applied post-hoc |
| Validation suite | **D** | 10/24 tests are tautologies or bounds; no conservation test |
| Documentation | **C−** | Good prose; benchmark claims are materially wrong |

**Recommendation:** issue a `v1.0.1` correctness release (P0 items in §4) before any external
citation, and do not advertise “validated against NIST / Hogeweg” until §5 is satisfied.

---

## 0.1 v1.0.1 correctness patch — status

The P0 patch is implemented on `fix/v1.0.1-correctness` (31 tests).  It replaces the
`φ²Ctρ` accumulation with the exact mass-conservative `Vp·(ρ−ρ_old)/dt`, solves the
well constraint inside the pressure system (active-set BHP/rate, no cell clamp),
adds Picard iterations with a consistent temperature update, fixes the 3D gravity sign,
makes the 0D step an exact scalar mass solve, corrects the Teesside 40–120 bar envelope,
and fixes the five-spot well modes plus the Python 3.9/CLI issues.

### Before / after

| Case | v1.0.0 reported | v1.0.0 actual ΔM | v1.0.1 reported | v1.0.1 ΔM | v1.0.1 MB |
|---|---:|---:|---:|---:|---:|
| 0D Teesside, 180 cycles | 0 t net (3.11 Mt “injected”) | −133.9 t (58 %) | +33.8 t net (53 928.9 t in / 53 895.1 t out) | +33.8 t | 1.5 × 10⁻¹² |
| 0D Okoroafor | −53 837.6 t | −34 377.3 t (+36.2 %) | −32 375.5 t | −32 375.5 t | 3.0 × 10⁻¹⁵ |
| 1D 90 d injection | 34.8 t | 30 919 t | 29 199.9 t | 29 199.9 t | 5.2 × 10⁻¹⁰ |
| 1D const BHP 60 bar, 10 d | 14.23 t | −22 958 t | 6 556.4 t (ref 6 575.0, −0.28 %) | −6 556.4 t | 7.6 × 10⁻¹⁰ |
| 2D 15×15 single well, full cycle | — | — | −406.6 t net | −406.6 t | 1.9 × 10⁻⁸ |
| 3D 8×8×5 single well, 30 d inject | 2 692.9 t | 32 500 t | 3 623.7 t | 3 623.7 t | 5.6 × 10⁻⁸ |

The 2D five-spot scenario is now physically meaningful (four corner injectors + one
central producer; v1.0.0 injected or withdrew all five wells together).  v1.0.1 gives
10 308.6 t injected, 3 995.0 t produced, MB 4.4 × 10⁻⁷.

### Convergence after the patch

| Case | Metric | v1.0.0 | v1.0.1 |
|---|---|---:|---:|
| 3D timestep, dt 3 → 1 d | working-gas spread | 4.35× | **−0.31 %** |
| 2D five-spot, dt 0.5 → 0.25 d | working-gas spread | 0.1 % (clamp-dominated) | **+1.4 %** |
| 2D grid, 10×10 → 30×30 | working gas | clamp-dominated | **−2.8 %** |
| 2D 40×40, default → 10× tighter SOR | working gas | 56 950 → 42 089 t (−26 %) | **same to 4 s.f. (0.00 %)** |
| 1D grid, N 10 → 80 | injected | pinned by the clamp | 29 155 → 29 128 t (0.09 %) |
| 3D hydrostatic, bottom − top | pressure | **−47 Pa** (top heavy) | **+4 682 Pa** |

The 40×40 SOR row compares convergence *within* each version; the absolute working gas
changed because the storage term (i.e. the physics) changed between versions.  The v1.0.1
default (`tol_Pa = 10⁻⁴`, 2 000 sweeps, 20 Picard iterations) matches a 10× tighter run to
four significant figures and closes mass to 1.3 × 10⁻⁷.  The 2D five-spot still shows a
~1.4 % timestep sensitivity between 0.5 and 0.25 d, and the 2D areal grid ~2.8 % between
10×10 and 30×30; those are P1 convergence items, not mass-balance errors.

### Remaining P1 (not in v1.0.1)

- **EOS accuracy (F-04):** Z error up to 4.2 % at 250 bar; Peneloux `c = −4.0 × 10⁻⁶ m³/mol`
  is verified to bring it below 1 %.  The README table is now honest about the 4.2 %.
- **Thermal model (F-12):** J-T error 37–148 %, sign wrong at 1 bar; the bulk isentropic
  update over-cools/over-heats the inventory.  Needs a real-gas enthalpy balance.
- **Hydrostatic initialisation (F-13b):** the model still starts from uniform pressure; the
  transient is now correct, but a hydrostatic initial state would remove the startup kick.
- **Tests (F-24):** the tautological tests were tightened or replaced.  The suite now has 31
  tests, including cross-tier mass balance, deliverability, gravity, and five-spot modes.

---

## 0.2 v1.1.0 release — status

The P1 thermophysical, hydrostatic, solver-validation, and configuration items above are now
resolved on `main` (`v1.1.0`, 38 tests + 29 subtests):

- **Peneloux Volume Translation (F-04, F-11):** Opt-in `peneloux=True` / `volume_translation = -4.0e-6`
  across all four model tiers and `config.py` reduces $Z$-factor error to $<0.85\%$ across
  1–350 bar without distorting cubic departure derivatives (`cp_real_h2`, `gamma_real_h2`,
  `fugacity_coefficient`).
- **Muzny Dense-Fluid Viscosity (F-12, F-21):** `viscosity_h2` now includes the Muzny et al.
  (2013) density residual ($\mu_0(T) + c_1\rho + c_2\rho^2$), matching CoolProp 8.0.0 within
  $<1.5\%$ across 250–450 K and 1–350 bar (mean $1.11\%$), and is used in 0D `_productivity_index`
  as well as 1D/2D/3D transmissibilities and well indices.
- **Analytical Joule-Thomson & Real-Gas $C_p$ (F-13, F-23):** Added `cp_real_h2` and
  `cp_specific_real_h2` ($1.81\%$ mean error vs CoolProp) and analytical $\mu_{\text{JT}}$ with
  Peneloux covolume translation ($0.74\%$ mean error vs CoolProp, negative across 1–350 bar).
- **3D Discrete Hydrostatic Initialisation (F-05, F-13b):** `UHSFieldScale3D._hydrostatic_init`
  initialises the vertical column around mid-depth $z_{\text{ref}} = 0.5 L_z$ using fixed-point
  iteration on discrete harmonic density $\rho_h$, so $\text{mean}(P) = P_{\text{init}}$ and
  shut-in vertical flux is zero to machine precision.
- **Solver & Config Hardening (F-16, F-17, F-24, F-25, F-26):** Added `max_picard_residual_Pa`
  and `converged` metrics to 1D/2D/3D results, post-solve `q_max` rate-cap checks, unique well
  name and grid-bounds validation, full 0D Monte Carlo schedule/property propagation, JSON/YAML
  config support across all scales, `CITATION.cff`, `CONTRIBUTING.md`, `NUMERICAL_TOLERANCES.md`,
  and `pytest-cov` in CI.

---

## 1. Method

1. **Static review** of every source file, test, CI workflow and metadata.
2. **Dynamic testing**: full test suite, CLI, examples, wheel build + clean-venv install.
3. **Thermophysical benchmark** against CoolProp HEOS, which implements the NIST reference EOS for
   normal hydrogen (Leachman et al., 2009), over 20 (T, P) points spanning 273–700 K, 1–350 bar.
4. **Flow benchmark** against:
   - exact material balance (inventory from the same EOS),
   - an independently written **mass-conservative implicit finite-volume radial simulator** on the
     same grid, same EOS, same well index (the “industrial” reference for this reduced problem),
   - timestep/grid refinement,
   - analytic hydrostatics.
5. **Reproduction harness:** `audit/bench_thermo.py`, `audit/bench_flow.py`, `audit/bench_1d_bhp.py`,
   `audit/proto_fix.py` (run with the package on `PYTHONPATH`).

A commercial-simulator cross-check (OPM Flow / Eclipse / GEM) was **not** possible in this
environment; it is the recommended next step (§5, L2).

---

## 2. Benchmark results

### 2.1 Thermophysics vs NIST (CoolProp HEOS)

20 points, 273–700 K, 1–350 bar:

| Property | Mean abs. error | Max abs. error | Where the max occurs |
|---|---:|---:|---|
| Z-factor | **2.41 %** | **4.68 %** | 293 K, 350 bar |
| Density | **2.49 %** | **4.91 %** | 293 K, 350 bar |
| cp° ideal gas (Shomate) | 0.05 % | 0.10 % | 1000 K |
| cp real gas | *not implemented* | −1.2 … −2.3 % | uses cp° only |
| Viscosity | **3.19 %** | **9.21 %** | 293 K, 350 bar |
| γ = cp/cv | 1.49 % | 2.84 % | 350 K, 250 bar |
| Joule–Thomson | **37–148 %** | sign wrong at 1 bar | 316 K, 80 bar: −0.0055 vs −0.0366 K/bar |

Representative points:

| T [K] | P [bar] | Z model | Z NIST | error |
|---:|---:|---:|---:|---:|
| 293.15 | 80 | 1.0272 | 1.0484 | −2.03 % |
| 293.15 | 100 | 1.0354 | 1.0608 | −2.40 % |
| 293.15 | 250 | 1.1094 | 1.1575 | −4.16 % |
| 316.0 | 80 | 1.0269 | 1.0459 | −1.82 % |
| 350.0 | 250 | 1.0974 | 1.1358 | −3.38 % |

**Consequence:** gas-in-place and working-gas capacities are systematically **2–5 % high** at
reservoir pressures (≈ +2.1 % at 80 bar, +2.8 % at 120 bar). That is outside the ±1 % usually
required for material balance and capacity booking. The README claims “Z < 0.25 %”; measured error is
**17× larger**.

**Cheap fix (verified):** a Peneloux volume translation with `c = −4.0 × 10⁻⁶ m³/mol`
(`Z_corr = Z_PR − cP/RT`) reduces the error over the same envelope to **mean 0.43 %, max 1.04 %**.
A better long-term fix is the Leachman (2009) Helmholtz EOS (≈0.1 %, still pure Python) or an
optional CoolProp/REFPROP backend.

### 2.2 The repository’s “NIST” table is not NIST

`thermo.py:417–428` and `tests/test_validation.py:76–87`:

| T [K] | P [bar] | repo “NIST” | CoolProp/NIST | repo value error |
|---:|---:|---:|---:|---:|
| 293.15 | 100 | 1.05881 | 1.06081 | −0.19 % |
| 293.15 | 200 | 1.11522 | 1.12473 | −0.85 % |
| 293.15 | 250 | 1.14283 | 1.15751 | **−1.27 %** |
| 316.00 | 80 | 1.02389 | 1.04594 | **−2.11 %** |
| 350.00 | 250 | 1.11665 | 1.13579 | −1.69 % |

The test tolerance is **3.0–3.5 %**, so the cubic EOS passes *despite* a −4.2 % error: at 316 K /
80 bar the model (−1.82 %) and the wrong reference (−2.11 %) agree with each other while both are
wrong. This is a calibration-to-error trap. The README table reproduces the same wrong values
(e.g. “Z 250 bar = 1.14283, error < 0.25 %”).

### 2.3 Material balance — every tier

Inventory is integrated with the same PR-EOS the code uses; `ΔM` is the true mass change,
`net_q` is the model’s reported injection minus withdrawal.

| Tier / case | Reported net | Actual ΔM | Relative error |
|---|---:|---:|---:|
| 0D Teesside (full schedule) | 0 t | **−133.9 t** | undefined (net 0) |
| 0D Okoroafor porous | −53 837.6 t | −34 377.3 t | **+36.2 %** |
| 1D Okoroafor, 90 d injection | 34.8 t | **30 919 t** | **889× under-report** |
| 1D full cycle | −0.7 t | −28 150 t | **4 × 10⁶ %** |
| 2D five-spot, 30 d injection | 119.4 t | **78 205 t** | **654× under-report** |
| 2D full cycle | −163.4 t | −96 930 t | **59 213 %** |
| 3D single-well, 30 d injection | 2 692.9 t | **32 500 t** | **12× under-report** |
| 3D full cycle | −1 348.6 t | −14 180 t | **951 %** |

Commercial simulators target **10⁻⁶–10⁻⁴** relative material-balance error. The gap is **8–12 orders
of magnitude**.

The 1D injection result is the clearest: the model raises a 500 m-radius reservoir from 80 to
**108 bar**, which requires **30 919 t** of hydrogen (31 501 t for a uniform 108 bar), yet reports
**34.8 t** injected — while the 0D tank model correctly says the same pore volume at 108 bar holds
≈ 125 400 t (93 900 t initial + 31 500 t). The three tiers therefore disagree by ~10³ on the same
physics.

### 2.4 Root cause: the accumulation term is a liquid term, counted twice

`one_d.py:216`, `two_d.py:207`, `three_d.py:277`:

```python
B[i] = phi * Ct * self.Vp[i] * rho_i / dt      # Vp already contains phi!
```

with `Vp = phi * V_cell` (`one_d.py:149–152`, `two_d.py:139`, `three_d.py:152`). Per bulk volume
this is

```
φ² · Ct · ρ · ∂P/∂t      instead of      φ · ρ · cg · ∂P/∂t
```

At 80 bar / 316 K (φ = 0.2, Ct = 5 × 10⁻¹⁰ Pa⁻¹, cg = 1.213 × 10⁻⁷ Pa⁻¹):

```
physical storage  φ·cg  = 2.43 × 10⁻⁸ Pa⁻¹
code storage      φ²·Ct = 2.00 × 10⁻¹¹ Pa⁻¹
ratio                         1 213×
```

The pressure diffusivity `η = k/(μ φ² Ct)` is therefore **1 213× too large**; the reservoir
equilibrates instantly to the well pressure. This is a *formulation* error, independent of the
well-model bug, and it invalidates BHP, deliverability, pressure-interference and working-gas
numbers in 1D/2D/3D.

`Ct` is a rock/water compressibility (0.5 × 10⁻⁴ bar⁻¹); a gas at 80 bar has
`cg ≈ 1.2 × 10⁻³ bar⁻¹`, i.e. **24× larger even before the φ² double count**.

### 2.5 Root cause: the well model is a clamp, not a boundary condition

`one_d.py:237–264`, `two_d.py:244–286`, `three_d.py:334–371`:

1. The linear solve inserts a Peaceman-like source `WI·(P_bhp − P_cell)`.
2. The resulting **cell pressure is then clamped to the BHP target** (`P_new = min(P_new, bhp)`).
3. The reported rate is recomputed from the *clamped* pressure:
   `q = WI · max(P_bhp − P_cell_clamped, 0)`, which is **0** whenever the clamp binds.
4. `q_max` is applied only to this reported rate, never to the solve.

So the mass that the solve actually added is **not** the mass that is reported. In the shipped 1D
run, the first step injects ~30 900 t, the clamp then pins the well cell at 108 bar, and every
subsequent step reports `q = 0`; the run reports 34.8 t total. The 2D/3D “BHP” is likewise the
clamped cell pressure and is **equal to the target by construction** — the reported BHP is not a
flowing bottom-hole pressure at all (no `BHP = P_cell − q/WI`).

**Fix:** choose one well constraint per well and enforce it in the solve (active set):
- rate-specified → Neumann source, mass balance closes by construction;
- BHP-specified → solve for `q` (scalar Newton) or add the Peaceman relation
  `q = WI·(P_cell − P_wf)` with `P_wf` fixed and no cell clamp.

### 2.6 Reference comparison, 1D, constant BHP = 60 bar, 10 days

| Model | Cumulative produced | ΔM | MB closure | P₀ | avg P |
|---|---:|---:|---:|---:|---:|
| Shipped v1.0.0 | 14.23 t | −22 958 t | −161 235 % | 60.00 bar | 60.00 bar |
| Prototype fix (EOS storage + no clamp) | 6 512.4 t | −6 561.8 t | **−0.76 %** | 60.20 bar | 67.97 bar |
| Mass-conservative reference | 6 578.4 t | −6 578.4 t | ~10⁻⁷ % | 60.15 bar | 67.94 bar |

The one-line storage-term replacement plus removal of the cell clamp reproduces the independent
reference to **1 %** — confirming the diagnosis and giving a concrete patch (`audit/proto_fix.py`).
The shipped model under-produces by **462×** and collapses the entire 500 m reservoir to the BHP
target instead of a near-well drawdown.

### 2.7 Timestep and grid convergence

| Case | Metric | Coarse | Fine | Spread |
|---|---|---:|---:|---:|
| 0D Teesside | WG, dt 2 → 0.25 d | 359.1 t | 347.7 t | 3.3 % |
| 3D single well | WG, dt 3 → 0.5 d | 31 256 t | 136 095 t | **4.35×** |
| 2D 15×15 five-spot | WG, dt 1 → 0.25 d | 175 026 t | 175 187 t | 0.1 % (clamp-dominated) |
| 2D 40×40, 30 d injection | WG, 80 → 3000 SOR iters | 56 950 t | 42 089 t | **−26 %** (still drifting at 400 iters: 62 794 t) |
| 1D radial | BHP at end, N = 10 → 80 | 55.00 bar | 55.00 bar | 0 (pinned by clamp) |

- The 3D working gas is **not converging**; 4.4× between dt = 3 and 0.5 d. A time step must be
  halved until the metric changes by <1 %; this model is not close.
- The default 2D `max_iter=80` (cells: `two_d.py:122`) silently returns a non-converged solution on
  a 40×40 grid; increasing to 3 000 iterations changes the answer by 35 %. There is no residual,
  iteration count or `converged` flag in the result.
- The 1D grid study is meaningless because the well cell is clipped to the BHP target; `N` cannot be
  assessed.

### 2.8 Storage capacity, 70 000 m³ cavern at 293.15 K

| P [bar] | ρ NIST | ρ model | error | inventory NIST | inventory model |
|---:|---:|---:|---:|---:|---:|
| 40 | 3.231 | 3.268 | +1.15 % | 226.2 t | 228.8 t |
| 80 | 6.311 | 6.442 | +2.07 % | 441.8 t | 450.9 t |
| 120 | 9.246 | 9.506 | +2.81 % | 647.2 t | 665.4 t |

Working gas 40–120 bar: **421.1 t (NIST)** vs 436.7 t (PR), **+3.7 %**. The shipped Teesside
scenario reports **356.1 t** because it actually runs **15–80 bar**, not the documented 40–120 bar:
`schedule_cavern` (`zero_d.py:370–384`) computes `P_max = min(2·P_ref, P_ref+80) = 80 bar` and
`P_min = max(0.35·P_ref, P_ref−25) = 15 bar` for `P_ref = 40 bar`.

The 70 000 m³ cavern has `q_max = 200 kg/s`; that is ~720 t/h, which would fill the entire cavern in
under an hour. Every 0.5-day step therefore saturates the BHP clamp.

### 2.9 Gravity in 3D has the wrong sign

Shut-in, uniform initial pressure, 6×6×5 grid (`audit/bench_flow.py`):

```
column bottom → top [bar]: 79.9764  79.9882  80.0000  80.0118  80.0236
```

Hydrostatics requires `P_bottom − P_top = ρgΔz > 0`; the model gives `P_top > P_bottom` by 47 Pa.
The RHS terms at `three_d.py:329–332` are

```python
rhs += TzD * (P_new[k-1] + gD)     # below
rhs += TzU * (P_new[k+1] - gU)     # above   ← sign inverted
```

So “buoyancy” pushes hydrogen **downward**. The unit test (`test_20`) only asserts
`P_bottom ≥ P_top − 500 Pa`, which passes with either sign — it cannot detect this. Also, the model
initialises a **uniform** pressure, not hydrostatic equilibrium, so it starts with a spurious
2.4 kPa perturbation.

### 2.10 Thermal model

| Run | Temperature range | Comment |
|---|---:|---|
| 0D Teesside | **209.2 – 397.1 K** | −84 K to +104 K around geothermal, in 12 h cycles |
| 2D five-spot | 259.7 – 344.1 K | 56 K cooling of the whole field |
| 3D single well | 308.9 – 321.6 K | mild |

The isentropic update is applied to the **entire inventory** on every pressure step
(`thermo.py:282–293`); for 15→80 bar it predicts 481 K. Real cavern gas stays within a few tens of
K of geothermal because (i) only the injected stream is compressed, (ii) the injected gas arrives at
compressor temperature, and (iii) the rock heat capacity buffers the cavern. The relaxation constant
τ = 3 d (`zero_d.py:122–124`) is a magic number; it should be derived from `UA/(m·cp)` with the
cavern wall area and rock conductivity, and the energy balance should use real-gas enthalpy, not a
constant γ.

The J-T coefficient itself is wrong by 37–148 % (§2.1), including a **positive** value at 293 K /
1 bar where the true sign is negative. Any wellbore-thermal or throttling calculation built on it
will be wrong.

### 2.11 Scenario semantics

- **“Five-spot”** (`two_d.py:469–481`): all five wells are driven by the *same* `ctrl.mode`, so all
  four corner wells inject and then all five produce. A true five-spot has four injectors and one
  producer simultaneously. There is no saturation or tracer, so “sweep” and “thief zone” cannot be
  represented at all — the single-phase pressure field is the only output.
- **“Working gas”** = `m_max − m_min` over the whole schedule, including the cushion-fill ramp. In
  Teesside this includes the 15 bar state, giving 79 % “efficiency” against a 40 bar design; the
  design working-gas fraction at 40–120 bar is ≈ 65 %.
- **“Cushion gas”** = `m_min`, i.e. 94.9 t at 15 bar, not the gas retained at the minimum operating
  pressure (226 t at 40 bar).
- **“cycled_gas_t”** (`two_d.py:405`, `three_d.py:538`) is set to `total_injected_t`, which is not
  cycled gas.
- Units are inconsistent across tiers: `max_pressure_Pa` (0D), `max_P_bar` (2D/3D), `mass_t`
  (2D/3D) vs `masses` in kg (0D/1D), BHP in Pa (1D) vs bar (2D/3D). This leaks into CSV headers and
  the CLI.

### 2.12 Performance

| Run | Time | Per cell-step |
|---|---:|---:|
| 2D 15×15, 188 steps, default | 1.4 s | ~0.03 ms |
| 3D 8×8×5, 35 steps, default | 1.8 s | ~0.16 ms |
| 2D 40×40, 30 steps, tight (3000 it) | 68.9 s | ~1.4 ms |

Extrapolating a 50×50×10 field (25 000 cells, 365 steps) at the tight setting: **several hours per
run** in pure Python. Fine for teaching; not field-scale. The EOS is re-evaluated per link per
iteration with no caching.

---

## 3. Findings register

Severity: **P0** = wrong results / false claims; **P1** = accuracy or API defect; **P2** = polish.

| ID | Sev | Location | Finding |
|---|---|---|---|
| F-01 | P0 | `one_d.py:216`, `two_d.py:207`, `three_d.py:277` | Gas storage term is `φ²Ctρ`; should be `φρcg`. **1 213× too stiff.** |
| F-02 | P0 | `one_d.py:248–264`, `two_d.py:250–286`, `three_d.py:344–371` | Cell pressure clamped to BHP; rate reported after clamp; `q_max` not in solve. Mass balance broken by 12–889×. |
| F-03 | P0 | `zero_d.py:199–238`, `265–307` | 0D linearises `dm/dP`, clamps, recomputes mass, but accumulates `q·dt` as throughput. Teesside reports 0 net while losing 134 t (58 %). |
| F-04 | P0 | `thermo.py:417–428`, `tests:76–87` | “NIST” table is off by up to 2.1 %; tolerances 3–3.5 % hide a 4.2 % model error. README claims <0.25 %. |
| F-05 | P0 | `three_d.py:329–332` | Gravity sign inverted; hydrostatic column wrong by 2×ρgΔz; test cannot detect it. |
| F-06 | P0 | `zero_d.py:314–331`, `370–384` | Teesside docstring says 40–120 bar; code runs 15–80 bar; `q_max` 200 kg/s is ~13× realistic. |
| F-07 | P0 | `two_d.py:469–499` | Five-spot wells all share one mode; no injector/producer split; no tracer/saturation. |
| F-08 | P0 | `tests/test_validation.py` | 10/24 tests are tautological/bounds (see Appendix A); no mass-balance, convergence or analytical test. |
| F-09 | P0 | `zero_d.py:516`, `one_d.py:411` | `zip(..., strict=False)` requires Python ≥3.10; project claims ≥3.9; CI 3.9 never exercises CSV export. |
| F-10 | P0 | `cli.py:219–232` | `h2multisim test` discovers `"tests"` relative to CWD → `ImportError` outside the repo (verified). |
| F-11 | P1 | `thermo.py:68–70` | PR EOS with constant κ = 0.024 has 2–5 % density error above 150 bar. Add Peneloux or Leachman EOS. |
| F-12 | P1 | `thermo.py:242–257` | Viscosity is pressure-independent (dilute-gas Chapman–Enskog); −9 % at 350 bar, −1.4 % at 80 bar. |
| F-13 | P1 | `thermo.py:261–278` | J-T coefficient wrong by 37–148 %, wrong sign at 1 bar; uses ideal Cp. |
| F-14 | P1 | `thermo.py:282–323` | Whole-inventory isentropic update; τ heuristic; no rock heat capacity; no enthalpy. |
| F-15 | P1 | `one_d.py:170–177`, `two_d.py:180–187`, `three_d.py:215–224` | Well index uses Cartesian `0.2Δr` on a radial log grid, clamped to `ln=0.1`; ad hoc, grid-dependent; no skin in 2D/3D. |
| F-15b | P1 | `three_d.py:357–371` | `q_max` is applied **per perforation**, so a well can deliver `nz × q_max`; “BHP” is the mid-perforation cell pressure, not a flowing BHP. |
| F-16 | P1 | `two_d.py:122`, `three_d.py:130` | Fixed 80 SOR iterations; silent non-convergence (35 % error at 40×40); no residual/flag. |
| F-17 | P1 | `zero_d.py:420–481` | MC ignores `pore_model`, `isotherm`, `m_sorbent`; Gaussian/truncated instead of LHS; no CI/sensitivity; CLI uses dt=2 d vs 1 d deterministic. |
| F-18 | P1 | all tiers | Result units/semantics inconsistent (Pa vs bar, kg vs t, K vs °C); metric names differ. |
| F-19 | P1 | `thermo.py:238`, `257`, `73`, `94`, `123` | Silent clamps (`gamma` to [1.10,1.80], μ to [6,25] µPa·s, T to 50 K, P to 1 Pa) hide invalid inputs. |
| F-20 | P1 | `cli.py:66, 163, 199` | Unknown `--scenario` silently falls through to a default (typos produce the wrong model). |
| F-21 | P1 | `zero_d.py:164–174` | PI uses fixed μ = 9 µPa·s and ignores `Z`, ρ(P) in the drawdown; no pseudopressure. |
| F-22 | P1 | repo root | 9 legacy scripts (~8 500 lines, filenames with spaces/parentheses) untested, unlinted, undocumented; README structure omits them. |
| F-23 | P1 | `thermo.py:200–238` | Real-gas cp not exposed/used; γ clamped; no speed of sound/enthalpy for an energy balance. |
| F-24 | P2 | `pyproject.toml` | `license = {text=...}` deprecated (PEP 639); no `[project.urls]`, no `py.typed`, no lockfile, no citation file. |
| F-25 | P2 | `.github/workflows/ci.yml` | No mypy, coverage, wheel-build/install smoke, or minimum-version runtime test. |
| F-26 | P2 | `two_d.py:244` | `next(...)` linear well scan inside the inner SOR loop (O(N·W·iters)). |
| F-27 | P2 | `one_d.py:118–120` | Thomas pivot clamp to 1e-30 returns silent garbage instead of raising. |
| F-28 | P2 | CLI/API | `print` instead of `logging`; no `--json`; no input validation. |
| F-29 | P2 | docs | No governing equations, assumptions, validity envelope, or benchmark provenance; “™” in source strings. |
| F-30 | P2 | metrics | No energy metrics (LHV, GWh, round-trip efficiency), which is what UHS studies actually report. |

**Status in v1.0.1:** F-01, F-02, F-03, F-05, F-06, F-07, F-08, F-09, F-10, F-16 and F-20
are fixed (§0.1).  F-04 is mitigated by honest references and tolerances, but the cubic EOS
itself is unchanged.  F-15 is fixed for 1D only (2D/3D still use the Cartesian `0.2Δ` well
index and no skin); F-02/F-15b’s 3D rate cap and BHP reporting are fixed.  The remaining P1/P2
items are unchanged.

---

## 4. Proposed fixes

> The P0 subsections below are the design record for the patch.  The shipped implementation
> uses the full Picard mass-balance loop rather than the one-line `cg` replacement shown in
> §4.1; compare with `one_d.py`, `two_d.py`, `three_d.py`, `zero_d.py` on the branch.

### 4.1 P0 core: mass-conservative accumulation (verified)

Replace the liquid storage term in 1D/2D/3D:

```python
# before (one_d.py:216)
B[i] = phi * Ct * self.Vp[i] * rho_i / dt

# after — gas storage, no double phi
dP = 1.0e3
cg = (self.eos.density(P[i] + dP, T[i])
      - self.eos.density(max(P[i] - dP, 1.0), T[i])) / (2.0 * dP)   # kg/m³/Pa
B[i] = self.Vp[i] * cg / dt
```

Better still, solve the mass equation with a Picard/Newton loop
`Vp_i (ρ_i^{n+1} − ρ_i^n)/Δt = Σ T_ij (P_j − P_i) + q_i` and iterate to `|ΔP| < 1 Pa`.
`audit/proto_fix.py` implements the one-line version; it reproduces the independent reference
within 1 % and closes mass balance to −0.76 % (vs −161 235 % shipped) on the 1D BHP case.

For 2D/3D the same replacement applies; the well cap must also be enforced (§4.2), otherwise the
field drains to the 1 Pa floor when the unlimited source exceeds `q_max`.

### 4.2 P0 well control

- Never clamp a cell pressure to a BHP target.
- For each well choose **one** constraint:
  - *rate control*: source term `q` (mass balance closes exactly);
  - *BHP control*: solve the Peaceman relation `q = WI (P_cell − P_wf)` and either fix `q` or
    Newton-iterate on `q` until `P_wf = target`; switch with an active set when `q_max` is reached.
- Report `BHP = P_cell − q/WI` (sign per mode), not `P_cell`.
- Use the correct radial index `WI = 2πkhρ/(μ (ln(r_c/r_w) − 1/2 + s))`, not `0.2Δr`; add skin and
  k-anisotropy for 2D/3D.

### 4.3 P0 exact 0D step

Use mass as the state and solve for `P` (or `q`) implicitly:

```python
def step(self, state, period, dt):
    def residual(P_new):
        T_new = coupled_temperature_update(state.T, state.P, P_new, ...)
        q = self._mass_rate(ReservoirState(P_new, T_new, 0.0), period)  # rate at new P
        return self.mass_from_state(P_new, T_new) - state.mass - q * dt
    P_new = brentq(residual, 1e3, 1e9)
    ...
```

No clamp, no linearisation, mass closes to machine precision. Sub-step to BHP targets rather than
overshooting in one 1-day step.

### 4.4 P0 metrics from state, not from `q·dt`

Add to every result:

```python
metrics["mass_balance_error"] = (m_final - m_initial - (inj - wdr)) / max(abs(inj - wdr), 1.0)
```

and derive `total_injected` / `total_withdrawn` from the state changes (or verify the two agree to
<1e-6). Define:

- `working_gas` = inventory at design `P_max` minus inventory at design `P_min`, evaluated on the
  last stable cycle;
- `cushion_gas` = inventory at design `P_min`;
- `working_gas_fraction` = WG / (WG + cushion), not “storage efficiency”;
- `cycled_gas` = Σ|q|dt, separate from `total_injected`;
- SI internally, with a unit-aware view for output.

### 4.5 P0 scenario corrections

- Teesside: explicit `P_min = 40 bar`, `P_max = 120 bar`, `q_max = 5–20 kg/s`, and a ramp to the
  target; keep the 70 000 m³ volume.
- Five-spot: per-well schedules (`INJ_*` inject while `CTRE` produces), or rename to
  “multi-well pressure interference”. Add a tracer/saturation model before using “sweep”.
- 3D gravity: fix the sign (`+g` for the upper link, `−g` for the lower), initialise hydrostatically,
  and add a test that asserts the analytic `ρgΔz` within 1 %.

### 4.6 P0 validation suite

- Move reference data to `tests/data/*.json` with provenance (source, conditions, uncertainty).
- Replace the 3–3.5 % Z tolerance with 0.5 % against true NIST values (and fail if the reference
  itself disagrees with CoolProp when available).
- Add: material-balance tests at 1e-8 for every tier; analytical Theis/line-source BHP test;
  hydrostatic test; timestep halving <1 %; grid refinement; MMS convergence for 1D/2D.
- Delete or rewrite the tautologies in Appendix A.

### 4.7 P0 packaging

- Either fix `zip(..., strict=False)` (drop the argument) or set `requires-python = ">=3.10"`;
  add a 3.9/3.10 job that calls `export_timeseries_csv` and `export_radial_timeseries_csv`.
- `h2multisim test`: ship `tests/` inside the package (`src/h2multisim/tests/`) and discover via
  `importlib.resources`, or invoke `pytest` on the installed package; return non-zero on collection
  errors.
- Validate `--scenario` with `choices` / a registry; error on typos.

### 4.8 P1 thermodynamics

1. Add a **Peneloux translation** `c = −4.0e-6 m³/mol` (0.43 % mean, 1.04 % max) as an immediate
   patch; expose it as `volume_translation` and document the fitted envelope.
2. Or implement **Leachman (2009) normal-hydrogen HEOS** (pure Python, ~0.1 %); add
   **Muzny (2013)** viscosity (pressure-dependent, ~1–2 %).
3. Use real-gas `cp`, `cv`, `h`, `s`, `w` (speed of sound) from the departure functions; remove the
   γ clamp; compute J-T from `(T ∂V/∂T − V)/cp` with the real cp.
4. Replace the lumped thermal update with an energy balance:
   `dU/dt = Σ q_in h_in − Σ q_out h_out + UA(T_rock − T) − P dV/dt`, with `UA` from cavern/reservoir
   geometry and rock conductivity; apply compression heating to the injected stream, not the whole
   inventory.
5. Add a `ThermoBackend` protocol (`ideal`, `PR`, `PR-Peneloux`, `table`, `coolprop`) so users can
   plug REFPROP without changing the models.

### 4.9 P1 numerics

- Replace fixed-iteration GS-SOR with a preconditioned Krylov solver (optional SciPy) or
  line-SOR/ADI with optimal ω (`ω_opt = 2/(1+sin(π/N))`); return `iterations`, `residual`,
  `converged`.
- Adaptive time stepping with a target ΔP per step and event handling for BHP/rate switches;
  sub-step to targets instead of overshooting.
- Cache ρ, μ, cg per cell per Newton iteration; vectorise with an optional NumPy backend.
- Add `mass_balance_error` and a `verify()` method to every model.

### 4.10 P1/P2 polish

- Unify `SchedulePeriod` / `WellControl` / `Schedule2D` / `Schedule3D` into one control type with
  per-well entries; unify results and metric names (SI internally).
- Add input validation (`k > 0`, `0 < φ ≤ 1`, `T > 0`, `P > 0`) with actionable errors; add a
  `strict=True` mode that raises instead of clamping.
- Add `logging`, `--json`, `--quiet`, and a `benchmark`/`verify` CLI command that prints the NIST
  table (from tabulated data, optionally recomputed with CoolProp).
- Energy metrics: LHV = 33.33 kWh/kg, GWh, round-trip efficiency, levelised cost hooks.
- Move the 9 legacy scripts to `legacy/` (or delete), exclude the directory in ruff, and document
  the canonical entry points; remove `key_improvements_in_v1.3.txt` from the root.
- Packaging: PEP 639 license string, `[project.urls]`, `py.typed`, `CITATION.cff`, `CHANGELOG.md`,
  `CONTRIBUTING.md`, `SECURITY.md`, issue templates, `uv.lock`, pre-commit.
- CI: `mypy` (at least `--check-untyped-defs`), coverage, `python -m build` + install + smoke run,
  a job on the minimum Python that exercises CSV export.
- Docs: add an “Assumptions & limitations” page (single-phase gas, no brine, no rel-perm, no
  capillarity, no dissolution, no hysteresis, no geochemistry), the governing equations, the
  validity envelope, and a validation page with provenance for every number.
- Correct the README benchmark table and remove “validated against Hogeweg” until a two-phase
  comparison exists (Hogeweg 2022 is a two-phase H₂-brine benchmark).

---

## 5. Industrial validation ladder

| Level | Test | Acceptance |
|---|---|---|
| L0 | Properties vs NIST/REFPROP tables (Z, ρ, cp, cv, w, μ, J-T) | Z, ρ ≤ 0.5 %; μ ≤ 2 %; cp, γ ≤ 1 % |
| L1 | 0D tank vs analytic isothermal/adiabatic; exact mass closure | MB error ≤ 1e-10; P ≤ 0.1 % |
| L2 | 1D radial vs Theis/line-source and vs OPM Flow/Eclipse (single-phase gas) | BHP ≤ 2 %; MB ≤ 1e-6 |
| L3 | 2D/3D vs SPE 1/2/9/10 adapted to gas; grid/timestep refinement | 2nd-order convergence; metrics <1 % |
| L4 | Two-phase H₂-brine vs Hogeweg et al. (2022), Feldmann et al. (2016), Lysyy et al. (2022) | pressure/saturation within published spread |
| L5 | Field history match (pressure, produced gas) | within gauge/measurement uncertainty |

Key references: Leachman et al. (2009); Muzny et al. (2013); Peneloux et al. (1982);
Peaceman (1978, 1983); Al-Hussainy et al. (1966); Aziz & Settari (1979); Hogeweg et al. (2022);
Okoroafor et al. (2022); Feldmann et al. (2016).

---

## 6. Roadmap

**v1.0.1 — correctness (1–2 weeks)**
F-01…F-10: storage term, well control, 0D exact step, metrics, gravity sign, Teesside/five-spot,
validation suite, 3.9/CLI. Re-run `audit/` and require MB < 1e-6 on all tiers.

**v1.1 — accuracy & usability (1–2 months)**
Peneloux/Leachman EOS, Muzny viscosity, real cp/J-T, energy balance, adaptive dt, SOR → Krylov,
LHS + bootstrap CI, unified API/units, JSON/logging, docs & CI hardening.

**v2.0 — industrial scope (3–6 months)**
Two-phase H₂-brine with rel-perm/capillarity/hysteresis/dissolution, optional NumPy/SciPy backend,
SPE/OPM benchmark suite, field-scale performance, validation paper.

---

## Appendix A — Test-by-test assessment

| Test | Verdict |
|---|---|
| 01 Z vs “NIST” | Fails to validate: wrong references, 3–3.5 % tolerance hides 4.2 % error |
| 02 density vs Okoroafor | ±5 % band; passes both 5.98 and 5.87 (can’t discriminate) |
| 03 density monotone | trivial |
| 04 Z bounds | trivial (0.05–5) |
| 05 viscosity | 5 %/7 % tolerances; pressure independence untested |
| 06 gamma bounds | **tautology** — function clamps to the same [1.10, 1.80] |
| 07 J-T sign | sign only; magnitude is 85 % wrong |
| 08 thermal | checks a formula against >50 K; no physics |
| 09 unit round-trip | fine |
| 10 Thomas | 6×6 only; no random/ill-conditioned cases |
| 11 BHP clamp | tests the clamp, not conservation |
| 12 rigid volume | fine |
| 13 pressure monotone | trivial |
| 14 mass conservation | **20 % tolerance** — not a conservation test (0D only) |
| 15 MC reproducibility | fine |
| 16 1D buildup | BHP is clamped to 108 bar → assertion is guaranteed |
| 17 1D drawdown | `BHP ≤ avg + 1 kPa` — trivially true (BHP is cell 0) |
| 18 5-spot vs single | arbitrary 0.5× threshold |
| 19 thief-zone CV | **tautology** (`cv ≥ 0`) |
| 20 gravity | passes with the inverted sign (47 Pa < 500 Pa tolerance) |
| 21 gravity override | arbitrary ordering of two wrong answers |
| 22 layered vs homogeneous | arbitrary ordering |
| 23 cross-scale EOS | **tautology** — same object compared to itself 3× |
| 24 cross-scale WG | only `> 0` |

---

## Appendix B — Reproduction

```bash
python -m venv .venv && . .venv/bin/activate
pip install -e . CoolProp numpy pytest
python -m pytest -q
python audit/bench_thermo.py      # NIST/CoolProp tables
python audit/bench_flow.py        # mass balance, dt/grid, reference
python audit/bench_1d_bhp.py      # deliverability vs reference
python audit/proto_fix.py         # verified 1D storage-term patch
```

Raw CLI outputs used above: `python -m h2multisim run --scale {0d,1d,2d,3d} ...`.

---

## Appendix C — What is genuinely good

- Clear four-tier structure and a coherent public API (`UHSReservoir0D`, `UHSRadial1D`,
  `UHSAreal2D`, `UHSFieldScale3D`).
- Closed-form Cardano/Viète cubic solver avoids Newton divergence; the gas-root selection is
  sensible.
- Zero runtime dependencies; wheel builds and installs cleanly; ruff and format pass; mypy is clean.
- Good docstrings and unit-conversion layer; the CSV dual-unit export is a nice touch.
- The 0D model’s `dm/dP` idea is sound — it just needs to be made exact and applied to the other
  tiers.
- The README’s multi-scale story is compelling; with the P0 fixes it can become a credible
  screening tool and teaching package.
