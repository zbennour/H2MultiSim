# Changelog

## v1.1.0 — Unreleased

### Thermodynamics & Accuracy
- Integrated Peneloux volume translation across all model tiers (0D, 1D, 2D, 3D) and config loaders (`peneloux=True` or `volume_translation=-4.0e-6`), reducing peak PR-EOS Z-factor error from 4.16% to <0.85% across 1–350 bar.
- Added Muzny et al. (2013) / NIST density-dependent dense-fluid residual correction to `viscosity_h2`, reducing high-pressure viscosity error from 9.2% down to <1.5% across 250–450 K and 1–350 bar.
- Added `cp_real_h2` and `cp_specific_real_h2` evaluating real-gas isobaric heat capacity via Peng-Robinson enthalpy departure.
- Resolved analytical Joule-Thomson coefficient formulation via $T(\partial V/\partial T)_P - V = (R T^2/P)(\partial Z/\partial T) + c$, eliminating numerical cancellation and restoring exact physical negative sign and magnitude.
- Expanded validation test suite to 33 test suites (29 subtests).

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
