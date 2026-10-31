# Changelog

## v1.3.1 — 2026-10-31

Fixes from the multi-lens audit of v1.3.0.

### Thermodynamics
- Fixed the sign of the Peng-Robinson residual heat capacity and added the exact $C_p - C_v$ term (`cp_real_h2`, `gamma_real_h2`). Real-gas $C_p$ now rises with pressure and matches CoolProp within 0.3 % (was −3.6 % at 300 bar); $\gamma$ error fell from +3.6 % to ≤ 2 %.
- `adiabatic_temperature` now uses the real-gas isentropic exponent $n = R(Z + T\,\partial Z/\partial T)/C_p$ instead of the ideal-gas $(\gamma-1)/\gamma$; isentropic $\Delta T$ error fell from 4–7 % to ≤ 1.6 % vs CoolProp.
- `joule_thomson_coefficient` no longer replaces a mixture EOS with pure H2 (the 80/20 H2/CH4 sign was wrong). The pure-H2 path's effective JT translation ($-8\times10^{-6}$ m³/mol) is now an explicit, documented empirical constant.
- Brooks-Corey `gas_saturation_from_capillary_pressure` honours `sgc`, so it inverts `capillary_pressure_Pa` exactly.

### Solvers
- 1D/2D/3D solvers now apply pore (rock) compressibility `Ct` ($V_p = V_{p0}e^{C_t(P-P_{init})}$), consistent with the 0D tier and the exported `ROCK` keyword. `Ct` was previously ignored.
- 2D/3D well index uses the Peaceman equivalent radius $0.14\sqrt{\Delta x^2+\Delta y^2}$ (was $0.28\sqrt{\Delta x\Delta y}$, ~6 % low productivity on the default grid).
- Well modes are validated (`inject`, `withdraw`, `shut_in`); a typo such as `"injection"` previously made a 1D well flow backwards.
- `run_monte_carlo` passes `swr`, `sgc` and `rel_perm` to each realisation (immobile water was ignored, overstating gas in place by 43 % at `swr = 0.3`).
- 1D mass balance includes constant-pressure outer-boundary influx (`boundary_influx_kg`); 1D temperature updates use the local mixture composition.
- 2D/3D viscosity uses the solver's EOS; new `sor_unconverged_solves` and `active_set_unconverged_steps` metrics feed `converged`.
- 0D reports `converged` / `unconverged_steps` from the root bracket; its inventory-based `mass_balance_error` is documented as a bookkeeping identity.

### Deck export
- `WCONPROD` writes the rate cap to GRAT (it was written to the THP item, so producer caps were lost).
- 3D decks write layers and perforations top-down (layered permeability and the gravity-override bottom perforations were flipped) and use isotropic `PERMZ` like the solver.
- `PORO`/`PERM*` include immobile-water factors; 1D `COMPDAT` carries the skin factor; deck headers document the remaining isothermal difference.

### Config & CLI
- Pressure strings accept `bar`, `bara`, `MPa`, `kPa`, `Pa`, `psi`, `psia`.
- 0D config pore-model `V0` defaults to the gas-accessible pore volume (honours `swr`).
- `--out-dir` and `--dt-days` override a config file's values instead of being silently replaced.

### Tests
- Added `tests/test_audit_fixes.py` (16 regression tests against CoolProp references and ECLIPSE item order); re-captured the 2D golden values.

## v1.2.0 — 2026-10-30

### Benchmarks & Verification
- Added **Tier 2 Analytical Pseudopressure Benchmark** (`audit/bench_analytical_1d.py` and `test_38`):
  - Solves the Al-Hussainy, Ramey & Crawford (1966) infinite-acting radial line source PDE via machine-precision pure-Python exponential integral $E_1(x)$ (series expansion for $x < 1$, continued fraction for $x \ge 1$) and adaptive Simpson quadrature real-gas pseudopressure $m(p) = 2 \int_0^p [p/(\mu Z)] dp$.
  - Verified `UHSRadial1D` against analytical solution: max flowing BHP error $0.216\%$ at $0.5\text{ d}$, improving to $0.010\%$ at $5.0\text{ d}$; max spatial profile error $0.108\%$ across 4 orders of magnitude in radius; mass conservation error $9.93 \times 10^{-13}$.
- Added comprehensive **3-Tier Industry Benchmark Report** (`BENCHMARKS.md`) documenting quantitative comparison against NIST / CoolProp 8.0.0 (HEOS / Leachman 2009 reference EOS), Al-Hussainy analytical solutions, and OPM Flow.

### Industry Simulator Bridge (ECLIPSE & OPM Flow)
- Added `h2multisim.deck_export` module exporting standard ECLIPSE 100/300 and OPM Flow `.DATA` input decks:
  - `generate_pvdg_table`: generates gas PVT table (`PVDG`: $P$, $B_g$, $\mu$) using pure-Python Peneloux PR-EOS and Muzny dense-fluid viscosity.
  - `export_eclipse_1d_radial`: generates native `RADIAL` logarithmic grid decks (`INRAD`, `DRV`, `DTHETAV`, `DZV`, `PORO`, `PERMR`, `PERMTHT`, `PERMZ`, `PVDG`, `PRESSURE`, `WELSPECS`, `COMPDAT`, `WCONINJE`, `WCONPROD`).
  - `export_eclipse_2d`: generates Cartesian areal decks (`DXV`, `DYV`, `DZV`, `PORO`, `PERMX`, `PERMY`, `PERMZ`, `PVDG`, `PRESSURE`, multi-well controls).
  - `export_eclipse_3d`: generates 3D stratigraphic field-scale decks with `EQUIL` hydrostatic datum matching mid-depth $z_{\text{ref}} = 0.5 L_z$ and layer perforations.
- Added `export-deck` CLI subcommand (`h2multisim export-deck --scale ... --scenario ... --out-file ...`).
- Added `test_39_eclipse_deck_export` in `tests/test_validation.py` bringing the test suite to 40 passing tests.
- Configured `[tool.pytest.ini_options]` with `pythonpath = ["src"]` in `pyproject.toml`.

## v1.1.0 — 2026-10-30

### Thermodynamics & Accuracy
- Integrated Peneloux volume translation across all model tiers (0D, 1D, 2D, 3D) and config loaders (`peneloux=True` or `volume_translation=-4.0e-6`), reducing peak PR-EOS Z-factor error from 4.16% to <0.85% across 1–350 bar.
- Fixed Peneloux volume-translation double-counting in `PengRobinsonH2.fugacity_coefficient` and evaluated cubic pressure derivatives at the untranslated molar volume $V_{\text{PR}} = V_m + c$ in `gamma_real_h2` and `cp_real_h2`.
- Added Muzny et al. (2013) / NIST density-dependent dense-fluid residual correction to `viscosity_h2`, reducing high-pressure viscosity error from 9.2% down to <1.5% across 250–450 K and 1–350 bar.
- Added `cp_real_h2` and `cp_specific_real_h2` evaluating real-gas isobaric heat capacity via Peng-Robinson enthalpy departure.
- Resolved analytical Joule-Thomson coefficient formulation via $T(\partial V/\partial T)_P - V = (R T^2/P)(\partial Z/\partial T) + c$ (using Peneloux covolume translation when `eos.volume_translation == 0.0`), restoring negative $\mu_{\text{JT}}$ across the entire 1–350 bar operating envelope.

### Solvers, Hydrostatics & Well Controls
- Rebuilt 3D `_hydrostatic_init` around the exact mid-depth datum ($z_{\text{ref}} = 0.5 L_z$) using fixed-point iteration on discrete harmonic density $\rho_h(P_k, P_{k+1})$ so $\text{mean}(P) = P_{\text{init}}$ and initial shut-in vertical flux is zero to machine precision.
- Added post-solve `q_max` rate-cap enforcement across 1D, 2D, and 3D active-set solvers and passed precomputed cell densities directly to `viscosity_h2`.
- Added upfront validation of unique well names, grid bounds, perforation intervals, and `k_field` dimensions in `UHSAreal2D` and `UHSFieldScale3D`.
- Added per-well operating mode overrides (`well_modes` and `mode_for`) to `Schedule3D`.
- Replaced hardcoded $9.0\,\mu\text{Pa}\cdot\text{s}$ viscosity in `UHSReservoir0D._productivity_index` with `viscosity_h2(T, P, rho=rho, eos=self.eos)`.
- Exposed `max_picard_residual_Pa` and `converged` metrics across 1D, 2D, and 3D results, and fixed float-accumulation snapshot boundary checks.

### Scenarios, Config, CLI & Testing
- Added `schedule_rough` in 0D representing Rough depleted gas field operating schedule (270d cushion gas fill, 3 annual cycles) and `rough_radial` in 1D (150 bar, 369 K).
- Added `compute_sweep_efficiency` in 2D supporting volumetric sweep fraction evaluation on both 2D grids and flat cell arrays.
- Updated `run_monte_carlo` and `cli.py` so `--scenario rough --mc` and `--config ... --mc` execute the active schedule and preserve `peneloux`, `pore_model`, `isotherm`, and `m_sorbent`.
- Extended `config.py` to support 0D `pore_model`/`isotherm`, 2D `k_field`/`k_field_mD`, 3D `k_layers`/`k_layers_mD`, and 3D `well_modes`.
- Added `NUMERICAL_TOLERANCES.md`, `CONTRIBUTING.md`, `CITATION.cff`, `tests/test_char_2d.py` characterization test, and expanded the test suite to 38 tests (with 29 subtests) and `pytest-cov` CI reporting.

## v1.0.1 — 2026-10-28

Correctness release. No public class names were removed; several result dictionaries
gained keys, and some default solver tolerances changed. See `AUDIT.md` for the full
audit and the before/after benchmark.

### Mass conservation
- Replaced the `φ²·Ct·ρ` accumulation in 1D/2D/3D with the exact mass-conservative
  `Vp·(ρ − ρ_old)/dt` term and Picard iteration.
- 0D `step()` now solves `m(P,T(P)) − m_old − q(P,T(P))·dt = 0` by a safeguarded
  Newton/bisection root; the reported rate is the actual inventory change.
- Added `mass_balance_error`, `initial_mass_kg/t`, `net_stored_kg/t`, and
  `cycled_gas_kg/t` to all result metrics. Measured closure is now `10⁻¹²`–`10⁻⁷`
  (was 12×–889× wrong).

### Wells
- Removed the post-solve cell-pressure clamp. BHP-controlled wells are solved inside
  the pressure system; `q_max` is enforced by an active set that switches the well to
  rate control and re-solves.
- Added a no-backflow constraint (an injection well shuts in if the reservoir exceeds
  its BHP target, and vice versa).
- 3D wells now use the total well rate cap and distribute a rate-controlled well over
  its perforations through a single implied BHP; reported BHP is no longer the cell
  pressure.
- 1D uses the innermost radial-cell radius as the Peaceman equivalent radius instead of
  `0.2·Δr`.

### Physics and scenarios
- Fixed the 3D gravity sign: the hydrostatic column is now bottom-heavy
  (`+4.68 kPa` over 100 m instead of `−47 Pa`).
- Teesside cavern now runs the documented 40–120 bar envelope (was 15–80 bar) with a
  15 kg/s well limit (was 200 kg/s).
- Five-spot and thief-zone scenarios now use four corner injectors and one central
  producer instead of applying one mode to all five wells.
- 0D schedule steps respect the requested period duration when `dt` is longer than the
  period (previously a 0.5 d period advanced by a full 1 d).

### Packaging / CLI / tests
- Version bumped to 1.0.1.
- Removed `zip(strict=...)`, restoring Python 3.9 compatibility.
- `h2multisim test` now locates the suite from a source checkout and prints a clear
  error (exit 2) when run from an installed wheel without tests.
- Unknown `--scenario` values are rejected per scale (exit 2).
- Tightened the NIST/CoolProp Z, density, viscosity, and γ tests to real reference
  values; added cross-tier mass-balance, constant-BHP deliverability, 3D gravity,
  five-spot well-mode, and no-backflow tests (31 tests total).
- `ruff` and `mypy` clean; `audit/` and Markdown are excluded from formatting.

### Known limitations (P1)
- Cubic PR-EOS Z error up to ~4.2 % at 250 bar; the Peneloux translation is verified
  but deferred.
- Joule–Thomson coefficient is 37–148 % off and has the wrong sign at 1 bar.
- Thermal model applies the bulk isentropic update to the whole inventory.
- Initial pressure is uniform, not hydrostatic.
