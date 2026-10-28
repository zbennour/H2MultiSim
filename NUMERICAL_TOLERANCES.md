# Numerical tolerances (H2MultiSim)

Screening-accuracy tolerances used by the validation suite (`tests/test_validation.py`).
Peng–Robinson cubic EOS is a screening model; values below are relative error vs
NIST/CoolProp 8.0.0 (Leachman et al. 2009 normal-H₂ reference EOS).

| Property | Condition | Reference | Tolerance |
|---|---|---|---|
| Z-factor | 293 K, 1–250 bar | CoolProp | 0.2%–5.0% (see cases) |
| Density | 80 bar, 316 K | 5.868 kg/m³ | ±0.15 kg/m³ (~2.5%) |

Solver convergence (Picard / Newton) uses relative residual ≤ 1e-6 on pressure
and saturation; mass balance is checked to ≤ 1e-8 relative per timestep.

See `AUDIT.md` (F-04) for the planned Peneloux volume-translation improvement.
