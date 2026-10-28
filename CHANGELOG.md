# Changelog

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
