# Audit Remediation Notes — v1.3.0 → v1.11.0 (October 2026)

This document records what changed after the multi-lens audit of v1.3.0 and **why**, release by
release, together with how each change was verified. The per-release details are in
[CHANGELOG.md](CHANGELOG.md); the validated accuracy envelope and remaining gaps are in
[LIMITATIONS.md](LIMITATIONS.md); cross-simulator and analytical benchmarks are in
[BENCHMARKS.md](BENCHMARKS.md).

## Guiding principles

1. **Measure before claiming.** Every accuracy statement is checked against an independent
   reference: CoolProp 8.0.0 (Leachman 2009 normal-H2 EOS and reference EOS for CH4/CO2/N2), the
   Al-Hussainy analytical solution, or OPM Flow. Several long-standing claims turned out to be
   wrong and were corrected rather than kept.
2. **Fix the physics, not the test.** When a reference disagreed, the root cause was found first;
   test reference values were regenerated only when the model *should* change (and the provenance
   is recorded next to the value).
3. **Keep the simple path identical.** New capabilities (cushion-gas transport, PCG, adaptive
   stepping) leave single-gas / fixed-step results bit-identical, and the previous behaviour is
   available through an explicit option (`linear_solver="sor"`, `adaptive=False`,
   `rock_heat_capacity=0`, `thermal_model="lumped"`).
4. **Every fix ships with a regression test**, and CI checks the things that broke silently
   before (OPM accepting the decks, examples running, `h2multisim verify`, mypy).

## Release-by-release summary

### v1.3.1 — audit fixes (PR #19)
| Change | Rationale |
|---|---|
| Real-gas Cp sign fixed and exact Cp−Cv term | The residual heat capacity had the wrong sign: Cp *fell* with pressure. Cp is now within 0.3 % of CoolProp. |
| Real-gas isentropic exponent in `adiabatic_temperature` | The ideal-gas (γ−1)/γ exponent was 4–7 % off for isentropic ΔT; now within 1.6 %. |
| JT honours mixture EOS | A mixture EOS was silently replaced by pure H2, giving the wrong sign for H2/CH4. |
| Rock compressibility `Ct` used by 1D/2D/3D | `Ct` was ignored by the solvers but exported to decks, so deck comparisons were not like-for-like. |
| Peaceman radius 0.14·√(Δx²+Δy²) | The previous radius made well productivity ~6 % low. |
| Well-mode validation | A typo such as `"injection"` made a 1D well flow backwards. |
| Monte Carlo passes `swr`/`sgc`/`rel_perm` | Immobile water was dropped, overstating gas in place by 43 % at `swr=0.3`. |
| Deck: WCONPROD GRAT item, 3D layer order, isotropic PERMZ | Producer rate caps were written to the THP item, and 3D decks were upside down (solver layer 0 = bottom, ECLIPSE K=1 = top). |
| Convergence metrics in every tier | Non-convergence was silent. |

### v1.3.2 — decks that actually run in OPM Flow (PR #20)
- OPM Flow has no single-phase gas mode, so every exported deck had always been rejected. Decks
  are now gas–water with immobile connate water (SGWFN from the solver's own rel-perm model).
- OPM builds a `RADIAL` grid but never connects the well to it, so the 1D model is exported as a
  Cartesian strip with the exact radial PORV/TRANX.
- `audit/bench_opm.py` runs every scenario in both simulators. Isothermal agreement: ≤ 1.2 % in
  stored-gas change, ≤ 0.05 % in pressure. The apparent 4–9 % gaps seen first were OPM's own
  large adaptive steps and a day-1 baseline in the comparison script, not model errors.
- CH4/CO2/N2 Peneloux shifts and viscosities refitted: CO2 viscosity error 62 % → 3 %.

### v1.3.3 — input validation, CI package and OPM jobs (PR #21)
- Non-physical inputs (k ≤ 0, φ ∉ (0,1], T ≤ 0, …) are rejected with the field name.
- CI builds and installs the wheel and runs the CLI, and installs OPM Flow to run every exported
  deck. That job later caught a real deck regression (v1.8.0).

### v1.4.0 — 0D energy balance (PR #22)
- The lumped thermal model treated filling a cavern as isentropic compression of the whole
  inventory, overstating heating by 55 % (+110 K vs +70.8 K CoolProp for 40→120 bar). The 0D
  tier now solves m·u balance with injected-gas enthalpy and wall heat transfer: +70.5 K
  (fill) and −79.4 K vs −81.2 K (discharge).
- The 0D Newton solver bisects when it stalls at the saturated-rate kink near a BHP target.

### v1.5.0 — PCG pressure solver, LIMITATIONS.md (PR #23)
- Jacobi-preconditioned CG with a true-residual convergence test replaced point-SOR as default.
  Combined with caching the per-step isentropic exponent: 2D 40×40 57 s → 31 s, 3D 54 s → 21 s,
  results unchanged to 1e-8 bar.
- LIMITATIONS.md states what is not modelled and the measured accuracy envelope.

### v1.6.0 / v1.7.0 — energy metrics; 2D and 3D cushion-gas mixing (PRs #24, #25)
- GWh (HHV/LHV) metrics in every tier, valued at the actual stream compositions.
- Binary H2 + cushion-gas (CH4/CO2/N2) transport in 2D and 3D, with gravity override in 3D.
  The pressure–composition coupling is stiff (H2 ~14× lighter than N2), so Anderson acceleration
  was added. Verified: mass balance ~1e-8, exact H2 species balance, symmetric plume, and a
  lock-exchange test showing H2 tonguing along the top and N2 slumping along the bottom.
- `min_withdrawn_h2_purity` started from the initial fraction (stuck at 0 for a pure cushion).

### v1.8.0 — `verify`, JSON output, viscosity refit (PR #26)
- `h2multisim verify` checks 64 reference points against embedded CoolProp values. On its first
  run it showed the documented "< 1.5 %" H2 viscosity was up to 3 % off (refitted: ≤ 0.5 %) and
  that the "< 0.85 %" Z claim held only on sparse points (measured: ≤ 1.0 % at 280–450 K,
  ≤ 1.25 % at 250–280 K; docs corrected).
- PVDG viscosities were rounded to 5 decimals, so OPM rejected the deck after the refit. The
  OPM CI job caught this, and values are now written at full precision.
- The BENCHMARKS Tier 1 table had a JT reference ~5× too large; it was regenerated.

### v1.9.0 — 0D adaptive sub-stepping (PR #27)
- With 1-day steps the Teesside cycle's 0.4-day withdrawals were one step each, so the peak
  temperature came out 12.5 K low. Sub-stepping to ≤ 10 bar / ≤ 5 K gives 86.5 °C vs 86.2 °C
  for a 0.02-day reference. `adaptive=False` reproduces fixed stepping exactly.

### v1.10.0 — rock heat capacity; typing and tooling; examples (PR #28)
- **Most consequential physics fix.** Porous-reservoir temperatures ignored the rock: compression
  heated only the gas, so the Rough case went from 96 to 138 °C. With local thermal equilibrium the
  rock (~2.2 MJ/m³K, ~100× the gas) absorbs the work and porous stores stay within ~1 K of
  geothermal, as observed in the field. Working gas rises 3–6 %. This also closed the 3–13 %
  "thermal" gap against isothermal OPM decks (now ≤ 1.25 %), which had wrongly been read as
  real physics. Caverns are unchanged.
- `py.typed`, mypy in CI, pre-commit, SECURITY.md; cushion-gas examples; every example runs in CI.

### v1.11.0 — real-gas entropy and compression energy (PR #29)
- `entropy_molar` and an exact isentropic solver (outlet T within ~1 K of CoolProp; the one-point
  estimate was 1–2.4 % off in work), intercooled multi-stage `compression_work_mass`, and 0D
  `compression_energy_GWh` (Teesside ~1.5 %, Okoroafor ~2.1 % of injected LHV).
- A "storage efficiency" metric was prototyped and dropped: it exceeded 100 % for schedules
  ending below the starting inventory, and without H2 losses it only restates compression cost.

## Behaviour changes to be aware of

| Since | Change | Restore old behaviour with |
|---|---|---|
| 1.4.0 | 0D thermal model is an energy balance | `thermal_model="lumped"` |
| 1.5.0 | 2D/3D pressure solver is PCG | `linear_solver="sor"` |
| 1.8.0 | H2 viscosity refit (~1–3 % higher at high pressure) | — (old fit was inaccurate) |
| 1.9.0 | 0D adaptive sub-stepping | `adaptive=False` / `--fixed-dt` |
| 1.10.0 | Rock heat capacity in porous media | `rock_heat_capacity=0` |

Reference values that moved with these changes were regenerated with recorded provenance: the
2D golden values and the 1D constant-BHP reference (6 496.9 t, from `audit/bench_1d_bhp.py`).

## Open items

- **Two-phase H2–brine flow** (saturation transport, capillarity, hysteresis, dissolution) is not
  modelled. Published two-phase benchmarks (e.g. Hogeweg et al. 2022) cannot be reproduced yet.
- **Heat advection and injected-gas temperature in 1D/2D/3D**: near-well thermal fronts are not
  represented.
- **Releases**: v1.3.1–v1.11.0 exist as merged PRs and version bumps; no git tags or GitHub
  releases have been created.
- CO2 density near its critical point (T < 310 K, P > 70 bar) is outside the validated envelope.
