# Numerical Tolerances (H2MultiSim)

Screening-accuracy tolerances used by the validation suite (`tests/test_validation.py`).
Values below are relative or absolute error bounds vs NIST/CoolProp 8.0.0 (Leachman et al. 2009 normal-H₂ reference EOS).

| Property | Condition | Reference | Tolerance |
|---|---|---|---|
| Z-factor (plain PR) | 293.15–400 K, 1–250 bar | CoolProp 8.0.0 | 0.2%–5.0% (`test_01`) |
| Z-factor (Peneloux, $c = -4.0\times 10^{-6}\text{ m}^3/\text{mol}$) | 293.15–400 K, 1–250 bar | CoolProp 8.0.0 | $\le 1.5\%$ (`test_32`, max $\sim 0.85\%$) |
| Density (plain PR) | 80 bar, 316 K | 5.868 kg/m³ | $\pm 0.15\text{ kg/m}^3$ ($\sim 2.5\%$, `test_02`) |
| Dynamic viscosity (Chapman-Enskog + Muzny) | 293.15–350 K, 1–250 bar | CoolProp 8.0.0 | $< 2.0\%$ (`test_05`, max $\sim 1.25\%$) |
| Heat capacity $C_p$ (PR departure) | 350 K, 250 bar | 14,769 J/(kg·K) | $< 4.0\%$ (`test_33`) |
| Isentropic exponent $\gamma$ (PR departure) | 250–600 K, 50–300 bar | CoolProp 8.0.0 | $< 8.0\%$ (`test_06`) |
| Joule-Thomson $\mu_{\text{JT}}$ (Peneloux) | 293.15 K, 100 bar | $-3.47\times 10^{-7}\text{ K/Pa}$ | $\pm 1.5\times 10^{-7}\text{ K/Pa}$ (`test_33`) |

## Solver & Mass-Balance Tolerances

- **0D Exact Mass-State Root Solver:** Bisection/secant pressure tolerance $\le 10^{-3}\text{ Pa}$; schedule mass balance relative error $< 10^{-9}$.
- **1D Implicit Radial Solver:** Picard pressure tolerance `picard_tol_Pa = 0.1 Pa`; Thomas tridiagonal residual $< 10^{-10}$; schedule mass balance relative error $< 10^{-6}$.
- **2D / 3D Implicit Finite-Volume Solvers:** Inner Gauss-Seidel SOR tolerance `tol_Pa = 1e-4 Pa`, outer Picard pressure tolerance `picard_tol_Pa = 1e-4 Pa`; schedule mass balance relative error $< 10^{-6}$.
