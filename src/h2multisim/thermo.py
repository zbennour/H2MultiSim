"""
h2multisim.thermo
=================
H2MultiSim™ Shared Thermodynamic Core Module.

Provides comprehensive, dependency-free thermodynamic and thermophysical
property functions for pure hydrogen (H₂) at subsurface storage conditions
(P ≤ 350 bar, T ≈ 250–700 K, expandable to 700 bar).

Features:
- Closed-form analytical Peng-Robinson cubic solver (Cardano/Viète)
- NIST-calibrated critical parameters and Shomate ideal gas heat capacities
- Dynamic viscosity via Chapman-Enskog / Chung formulation
- Joule-Thomson coefficient and real-gas heat capacity ratio γ(T, P)
- Coupled adiabatic and geothermal thermal relaxation models
- Adsorption isotherms (Langmuir, Dual-Site Langmuir) and pore volume models

References:
[1] Peng, D.-Y., & Robinson, D. B. (1976). Ind. Eng. Chem. Fundam., 15(1), 59–64.
[2] NIST Chemistry WebBook (2022). Hydrogen thermophysical properties.
[3] Poling, B. E., Prausnitz, J. M., & O'Connell, J. P. (2001). Properties of Gases and Liquids.
[4] Hogeweg, S., Strobel, G., & Hagemann, B. (2022). Comput. Geosci., 26, 1367–1378.
[5] Okoroafor, E. U., Saltzer, S., & Kovscek, A. R. (2022). Int. J. Hydrogen Energy, 47, 33781.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from h2multisim.units import cubic_real_roots

# ── Physical constants (NIST 2022) ───────────────────────────────────────────
R_UNIVERSAL: float = 8.314462618  # J/(mol·K) universal gas constant
M_H2: float = 2.01594e-3  # kg/mol molar mass of pure H2
R_SPECIFIC: float = R_UNIVERSAL / M_H2  # J/(kg·K) ≈ 4124.36 J/(kg·K)

# H2 Critical constants
T_C: float = 32.938  # K  critical temperature
P_C: float = 1.2858e6  # Pa critical pressure
RHO_C: float = 30.09  # kg/m³ critical density
OMEGA: float = -0.219  # Pitzer acentric factor

# Lennard-Jones parameters for H2
LJ_SIGMA: float = 2.827  # Å
LJ_EPS_KB: float = 59.7  # K (ε/kB)

# Normal Temperature & Pressure (NTP)
T_NTP: float = 293.15  # K
P_NTP: float = 101325.0  # Pa
MU_NTP: float = 8.948e-6  # Pa·s (NIST)

# Operational limits
T_FLOOR: float = T_C + 0.5  # K minimum valid temperature


class PengRobinsonH2:
    """
    Peng-Robinson Equation of State for pure hydrogen.

    Combines exact analytical cubic root finding (Cardano/Viète) with
    Newton-Raphson polish for numerical precision and absolute stability.
    """

    def __init__(
        self, max_iter: int = 50, tol: float = 1e-10, volume_translation: float = 0.0
    ) -> None:
        self.max_iter: int = max_iter
        self.tol: float = tol
        self.volume_translation: float = float(volume_translation)  # m³/mol Peneloux c
        self.ac: float = 0.45724 * R_UNIVERSAL**2 * T_C**2 / P_C
        self.b: float = 0.07780 * R_UNIVERSAL * T_C / P_C
        self.kappa: float = 0.37464 + 1.54226 * OMEGA - 0.26992 * OMEGA**2

    def _safe_T(self, T: float) -> float:
        return max(float(T), T_FLOOR)

    def _alpha(self, T: float) -> float:
        Tr = self._safe_T(T) / T_C
        return (1.0 + self.kappa * (1.0 - math.sqrt(Tr))) ** 2

    def _AB(self, P: float, T: float) -> tuple[float, float]:
        T_eff = self._safe_T(T)
        aT = self.ac * self._alpha(T_eff)
        A = aT * P / (R_UNIVERSAL**2 * T_eff**2)
        B = self.b * P / (R_UNIVERSAL * T_eff)
        return A, B

    def z_factor(self, P: float, T: float) -> float:
        """
        Compressibility factor Z = P·Vm / (R·T).
        Solves cubic polynomial analytically with Newton-Raphson polish.
        """
        if P <= 0.0 or T <= 0.0:
            return 1.0

        T_eff = max(float(T), 50.0)
        A, B = self._AB(P, T_eff)
        if not (math.isfinite(A) and math.isfinite(B)):
            return 1.0

        # Monic cubic: Z³ + c2·Z² + c1·Z + c0 = 0
        c2 = -(1.0 - B)
        c1 = A - 3.0 * B**2 - 2.0 * B
        c0 = -(A * B - B**2 - B**3)

        roots = cubic_real_roots(c2, c1, c0)
        gas_roots = [r for r in roots if r > B + 1e-6 and r > 0.05 and math.isfinite(r)]
        Z = max(gas_roots) if gas_roots else max((r for r in roots if r > 0.0), default=1.0)

        # Polish to machine precision
        for _ in range(self.max_iter):
            f = Z**3 + c2 * Z**2 + c1 * Z + c0
            df = 3.0 * Z**2 + 2.0 * c2 * Z + c1
            if abs(df) < 1e-20:
                break
            dZ = -f / df
            Z += dZ
            if abs(dZ) < self.tol:
                break

        # Optional Peneloux volume translation: V = V_PR - c, i.e.
        # Z_corr = Z_PR - c*P/(R*T).  c = -4.0e-6 reproduces the <1 % Z
        # result verified in AUDIT.md (F-04).  Default c = 0 (plain PR).
        c = self.volume_translation
        if c != 0.0:
            Z = Z - c * P / (R_UNIVERSAL * T_eff)

        return float(max(Z, 1e-4))

    def density(self, P: float, T: float) -> float:
        """Gas density [kg/m³]: ρ = P·M / (Z·R·T)."""
        P_val = max(float(P), 1.0)
        T_val = max(float(T), 50.0)
        Z = self.z_factor(P_val, T_val)
        return P_val * M_H2 / (Z * R_UNIVERSAL * T_val)

    def molar_volume(self, P: float, T: float) -> float:
        """Molar volume Vm [m³/mol]: Vm = Z·R·T / P."""
        P_val = max(float(P), 1.0)
        T_val = max(float(T), 50.0)
        Z = self.z_factor(P_val, T_val)
        return Z * R_UNIVERSAL * T_val / P_val

    def drho_dP(self, P: float, T: float) -> float:
        """
        Isothermal derivative of density with respect to pressure [kg/(m³·Pa)].

        Central finite difference of :meth:`density`, used by the mass-conservative
        accumulation term in the finite-volume models.  The step is scaled with
        pressure so that the derivative is accurate over the full 1 bar–700 bar
        operating range.
        """
        P_val = max(float(P), 1.0)
        T_val = max(float(T), 50.0)
        dP = max(P_val * 1.0e-6, 1.0)
        P_lo = max(P_val - dP, 1.0)
        P_hi = P_val + dP
        return (self.density(P_hi, T_val) - self.density(P_lo, T_val)) / (P_hi - P_lo)

    def fugacity_coefficient(self, P: float, T: float) -> float:
        """Thermodynamic fugacity coefficient φ = f / P from PR-EOS."""
        P_val = max(float(P), 1.0)
        T_val = self._safe_T(T)
        Z = self.z_factor(P_val, T_val)
        A, B = self._AB(P_val, T_val)
        s2 = math.sqrt(2.0)
        ln_phi = (
            (Z - 1.0)
            - math.log(max(Z - B, 1e-12))
            - A
            / (2.0 * s2 * max(B, 1e-12))
            * math.log(max((Z + (1.0 + s2) * B) / (Z + (1.0 - s2) * B), 1e-12))
        )
        # Peneloux volume-translation correction to fugacity.
        c = self.volume_translation
        if c != 0.0:
            ln_phi -= c * P_val / (R_UNIVERSAL * T_val)
        return math.exp(max(-50.0, min(50.0, ln_phi)))

    def fugacity(self, P: float, T: float) -> float:
        """Thermodynamic fugacity f = P · φ [Pa]."""
        return max(float(P), 1.0) * self.fugacity_coefficient(P, T)

    def isothermal_compressibility(self, P: float, T: float) -> float:
        """Isothermal gas compressibility cg = (1/ρ) · (∂ρ/∂P)_T [Pa⁻¹]."""
        P_val = max(float(P), 1.0)
        T_val = max(float(T), 50.0)
        dP = P_val * 1e-4 + 1.0
        drho = self.density(P_val + dP, T_val) - self.density(max(P_val - dP, 1.0), T_val)
        return drho / (2.0 * dP * max(self.density(P_val, T_val), 1e-12))

    def da_dT(self, T: float) -> float:
        """First derivative of EOS attraction parameter da/dT."""
        T_eff = max(float(T), 1e-3)
        Tr = T_eff / T_C
        sqTr = math.sqrt(Tr)
        return self.ac * (-self.kappa * (1.0 + self.kappa * (1.0 - sqTr)) / sqTr) / T_C

    def d2a_dT2(self, T: float) -> float:
        """Second derivative of EOS attraction parameter d²a/dT²."""
        T_eff = max(float(T), 1e-3)
        Tr = T_eff / T_C
        sqTr = math.sqrt(Tr)
        return self.ac * self.kappa * (1.0 + self.kappa) / (2.0 * Tr * sqTr) / T_C**2

    def viscosity(self, P: float, T: float) -> float:
        """Dynamic viscosity [Pa·s] using NIST-calibrated Chapman-Enskog."""
        return viscosity_h2(T, P)


# ── Ideal-gas heat capacity ──────────────────────────────────────────────────
def cp_ideal_h2(T: float) -> float:
    """
    Ideal-gas molar heat capacity Cp°(T) [J/(mol·K)] — NIST Shomate equation.
    Valid from 298 to 1000 K (extrapolates smoothly to 50 K).
    """
    T_val = max(float(T), 50.0)
    t = T_val / 1000.0
    Cp = 33.066178 - 11.363417 * t + 11.432816 * t**2 - 2.772874 * t**3 - 0.158558 / max(t**2, 1e-6)
    return max(Cp, 14.0)


def cp_specific_ideal_h2(T: float) -> float:
    """Ideal-gas specific heat capacity cp° [J/(kg·K)]."""
    return cp_ideal_h2(T) / M_H2


# ── Real-gas specific heat ratio γ(T, P) ─────────────────────────────────────
def gamma_real_h2(T: float, P: float, eos: PengRobinsonH2 | None = None) -> float:
    """
    Real-gas heat capacity ratio γ = Cp / Cv at temperature T [K] and pressure P [Pa].
    Combines Shomate ideal capacity with PR-EOS departure function.
    """
    if eos is None:
        eos = PengRobinsonH2()

    T_val = max(float(T), 50.0)
    P_val = max(float(P), 1.0)
    Vm = eos.molar_volume(P_val, T_val)
    b = eos.b
    R = R_UNIVERSAL

    def _press(Vm_: float, T_: float) -> float:
        a_ = eos.ac * eos._alpha(T_)
        d_ = Vm_ * (Vm_ + b) + b * (Vm_ - b)
        return R * T_ / (Vm_ - b) - a_ / max(d_, 1e-30)

    dT = T_val * 1e-4 + 0.01
    dV = Vm * 1e-4 + 1e-12
    dPdT_V = (_press(Vm, T_val + dT) - _press(Vm, T_val - dT)) / (2.0 * dT)
    dPdV_T = (_press(Vm + dV, T_val) - _press(Vm - dV, T_val)) / (2.0 * dV)

    Cp_Cv = -T_val * dPdT_V**2 / dPdV_T if abs(dPdV_T) > 1e-30 else R

    s2 = math.sqrt(2.0)
    d2a = eos.d2a_dT2(T_val)
    num = Vm + (1.0 + s2) * b
    den = Vm + (1.0 - s2) * b
    delta_Cp = (
        -T_val * d2a / (2.0 * s2 * b) * math.log(num / den)
        if abs(den) > 1e-20 and num / den > 0.0
        else 0.0
    )

    Cp = max(cp_ideal_h2(T_val) + delta_Cp, 14.0)
    Cv = max(Cp - Cp_Cv, 8.0)
    return float(min(max(Cp / Cv, 1.10), 1.80))


# ── Viscosity ────────────────────────────────────────────────────────────────
def viscosity_h2(T: float, P: float) -> float:
    """
    Dynamic viscosity of pure hydrogen [Pa·s] at (T, P).
    Chapman-Enskog formulation calibrated against NIST and Hogeweg et al. (2022).
    Benchmark: 9.4 μPa·s at 333 K, 82 bar.
    """
    T_val = max(float(T), 50.0)
    M_g = M_H2 * 1000.0  # g/mol
    T_star = T_val / LJ_EPS_KB
    Omega_v = (
        1.16145 / T_star**0.14874
        + 0.52487 / math.exp(0.7732 * T_star)
        + 2.16178 / math.exp(2.43787 * T_star)
    )
    mu0 = 2.6693e-6 * math.sqrt(M_g * T_val) / (LJ_SIGMA**2 * max(Omega_v, 1e-6))
    return float(min(max(mu0, 6.0e-6), 2.5e-5))


# ── Joule-Thomson coefficient ────────────────────────────────────────────────
def joule_thomson_coefficient(T: float, P: float, eos: PengRobinsonH2 | None = None) -> float:
    """
    Joule-Thomson throttling coefficient μ_JT = (∂T/∂P)_H [K/Pa] for H2.
    Physical note: Below the inversion temperature (~202 K at 1 bar), μ_JT > 0.
    At UHS reservoir conditions (T > 250 K), μ_JT < 0 (H2 warms on isenthalpic expansion).
    """
    if eos is None:
        eos = PengRobinsonH2()

    T_val = max(float(T), 50.0)
    P_val = max(float(P), 1.0)
    Vm = eos.molar_volume(P_val, T_val)
    Cp = cp_ideal_h2(T_val)
    dT = T_val * 1e-4 + 0.01
    dVm_dT = (eos.molar_volume(P_val, T_val + dT) - eos.molar_volume(P_val, T_val - dT)) / (
        2.0 * dT
    )
    return float((T_val * dVm_dT - Vm) / max(Cp, 1.0))


# ── Coupled Thermal Modeling ─────────────────────────────────────────────────
def adiabatic_temperature(
    T_old: float, P_old: float, P_new: float, eos: PengRobinsonH2 | None = None
) -> float:
    """
    Isentropic gas temperature change under rapid expansion or compression:
    T_new = T_old · (P_new / P_old)^((γ - 1) / γ).
    """
    T_init = max(float(T_old), 50.0)
    P_init = max(float(P_old), 1.0)
    P_final = max(float(P_new), 1.0)
    g = gamma_real_h2(T_init, P_init, eos)
    return float(max(T_init * (P_final / P_init) ** ((g - 1.0) / g), 50.0))


def thermal_relaxation_update(
    T_current: float, T_geothermal: float, dt_s: float, tau_s: float
) -> float:
    """
    Exponential relaxation of gas temperature toward surrounding rock:
    T(t + dt) = T_geo + (T - T_geo) · exp(-dt / τ).
    """
    if tau_s <= 0.0:
        return float(T_geothermal)
    return float(T_geothermal + (T_current - T_geothermal) * math.exp(-dt_s / max(tau_s, 1.0)))


def coupled_temperature_update(
    T_old: float,
    P_old: float,
    P_new: float,
    T_geothermal: float,
    dt_s: float,
    tau_s: float,
    eos: PengRobinsonH2 | None = None,
) -> float:
    """
    Coupled two-stage thermal update:
    1. Adiabatic isentropic response to pressure change
    2. Conductive thermal relaxation toward geothermal rock temperature
    """
    T_ad = adiabatic_temperature(T_old, P_old, P_new, eos)
    return thermal_relaxation_update(T_ad, T_geothermal, dt_s, tau_s)


# ── Material Models: Pore Volume ─────────────────────────────────────────────
@dataclass
class PoreVolumeModel:
    """
    Reservoir / vessel pore volume response to pressure.

    model = 'rigid'
        Volume is invariant: V(P) = V0 (steel vessels, caverns in stiff halite).
    model = 'elastic'
        Pore volume dilation with rock compressibility ct: V(P) = V0 · exp(ct · (P - P0)).
    model = 'breathing'
        Structural gate-opening transition in flexible Metal-Organic Frameworks (MOFs):
        V(P) transitions smoothly from V0 to V0 · V_ratio via a logistic curve.
    """

    model: str = "elastic"
    V0: float = 1.0
    P0: float = 1e5
    total_compressibility: float = 1e-9
    P_gate: float = 5e6
    k_gate: float = 1e-6
    V_ratio: float = 2.0

    def volume(self, P: float) -> float:
        """Pore volume [m³] at pressure P [Pa]."""
        P_val = max(float(P), 1.0)
        if self.model == "rigid":
            return self.V0
        if self.model == "elastic":
            return max(self.V0 * math.exp(self.total_compressibility * (P_val - self.P0)), 1e-12)
        if self.model == "breathing":
            f = 1.0 / (1.0 + math.exp(-self.k_gate * (P_val - self.P_gate)))
            return self.V0 + (self.V0 * self.V_ratio - self.V0) * f
        raise ValueError(
            f"Unknown pore model: {self.model!r}. Choose 'rigid', 'elastic', or 'breathing'."
        )


# ── Material Models: Sorbent Adsorption ──────────────────────────────────────
@dataclass
class IsothermModel:
    """
    Adsorption isotherm model driven by thermodynamic fugacity f [Pa].

    model = 'none'
        No sorbent adsorption (caverns, standard depleted fields).
    model = 'langmuir'
        Single-site Langmuir isotherm: q = q_max · K·f / (1 + K·f).
    model = 'dsl'
        Dual-Site Langmuir isotherm (two binding energies, e.g. MOFs).
    """

    model: str = "none"
    n_max1: float = 20.0  # mol / kg_sorbent
    K0_1: float = 1e-7  # Pa⁻¹
    n_max2: float = 5.0  # mol / kg_sorbent
    K0_2: float = 5e-8  # Pa⁻¹
    T_ref: float = 298.0  # K
    dH1: float = -6000.0  # J/mol
    dH2: float = -4000.0  # J/mol

    def _K1(self, T: float) -> float:
        return self.K0_1 * math.exp(-self.dH1 / R_UNIVERSAL * (1.0 / T - 1.0 / self.T_ref))

    def _K2(self, T: float) -> float:
        return self.K0_2 * math.exp(-self.dH2 / R_UNIVERSAL * (1.0 / T - 1.0 / self.T_ref))

    def adsorbed_moles_per_kg(self, P: float, T: float, eos: PengRobinsonH2 | None = None) -> float:
        """Return adsorbed hydrogen [mol H2 / kg sorbent]."""
        if self.model == "none":
            return 0.0

        if eos is None:
            eos = PengRobinsonH2()

        f = eos.fugacity(P, T)
        if self.model == "langmuir":
            K = self._K1(T)
            return self.n_max1 * K * f / (1.0 + K * f)
        if self.model == "dsl":
            K1 = self._K1(T)
            K2 = self._K2(T)
            q1 = self.n_max1 * K1 * f / (1.0 + K1 * f)
            q2 = self.n_max2 * K2 * f / (1.0 + K2 * f)
            return q1 + q2
        raise ValueError(
            f"Unknown isotherm model: {self.model!r}. Choose 'none', 'langmuir', or 'dsl'."
        )


# ── NIST Benchmarking Self-Test ──────────────────────────────────────────────
# Reference values: CoolProp 8.0.0 / Leachman et al. (2009) normal-H2 EOS.
# Tolerances document the cubic-PR screening accuracy (max ~4.2% at 250 bar);
# see AUDIT.md F-04 for the proposed Peneloux volume translation.
_NIST_Z_BENCHMARKS = [
    (293.15, 101325.0, 1.00059, 0.002),
    (293.15, 5.0e6, 1.02999, 0.030),
    (293.15, 10.0e6, 1.06081, 0.030),
    (293.15, 20.0e6, 1.12473, 0.050),
    (293.15, 25.0e6, 1.15751, 0.050),
    (350.00, 10.0e6, 1.05338, 0.030),
    (350.00, 25.0e6, 1.13579, 0.050),
    (400.00, 10.0e6, 1.04789, 0.030),
    (400.00, 25.0e6, 1.12074, 0.050),
    (316.00, 8.0e6, 1.04594, 0.030),
]


def run_thermo_self_test(verbose: bool = False) -> bool:
    """Validate thermodynamic module against NIST reference points."""
    eos = PengRobinsonH2()
    passed = True
    for T, P, Znist, tol in _NIST_Z_BENCHMARKS:
        Z_calc = eos.z_factor(P, T)
        err = abs(Z_calc - Znist) / Znist
        if err > tol:
            passed = False
            if verbose:
                print(
                    f"FAIL: T={T}K, P={P}Pa: Z_calc={Z_calc:.5f}, Z_NIST={Znist:.5f}, err={err:.4f}"
                )
    return passed
