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
T_FLOOR: float = 50.0  # K minimum valid temperature for EOS & Shomate correlations


class PengRobinsonH2:
    """
    Peng-Robinson Equation of State for pure hydrogen.

    Combines exact analytical cubic root finding (Cardano/Viète) with
    Newton-Raphson polish for numerical precision and absolute stability.
    """

    def __init__(
        self,
        max_iter: int = 50,
        tol: float = 1e-10,
        volume_translation: float = 0.0,
        peneloux: bool = False,
    ) -> None:
        self.max_iter: int = max_iter
        self.tol: float = tol
        if peneloux and volume_translation == 0.0:
            volume_translation = -4.0e-6
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
        c = self.volume_translation
        Z_pr = Z + c * P_val / (R_UNIVERSAL * T_val)
        A, B = self._AB(P_val, T_val)
        s2 = math.sqrt(2.0)
        ln_phi = (
            (Z_pr - 1.0)
            - math.log(max(Z_pr - B, 1e-12))
            - A
            / (2.0 * s2 * max(B, 1e-12))
            * math.log(max((Z_pr + (1.0 + s2) * B) / (Z_pr + (1.0 - s2) * B), 1e-12))
        )
        # Peneloux volume-translation correction to fugacity (applied once).
        if c != 0.0:
            ln_phi -= c * P_val / (R_UNIVERSAL * T_val)
        return math.exp(max(-50.0, min(50.0, ln_phi)))

    def fugacity(self, P: float, T: float) -> float:
        """Thermodynamic fugacity f = P · φ [Pa]."""
        return max(float(P), 1.0) * self.fugacity_coefficient(P, T)

    def isothermal_compressibility(self, P: float, T: float) -> float:
        """Isothermal gas compressibility cg = (1/ρ) · (∂ρ/∂P)_T [Pa⁻¹]."""
        P_val = max(float(P), 1.0)
        T_val = self._safe_T(T)
        dP = P_val * 1e-4 + 1.0
        drho = self.density(P_val + dP, T_val) - self.density(max(P_val - dP, 1.0), T_val)
        return drho / (2.0 * dP * max(self.density(P_val, T_val), 1e-12))

    def da_dT(self, T: float) -> float:
        """First derivative of EOS attraction parameter da/dT."""
        T_eff = self._safe_T(T)
        Tr = T_eff / T_C
        sqTr = math.sqrt(Tr)
        return self.ac * (-self.kappa * (1.0 + self.kappa * (1.0 - sqTr)) / sqTr) / T_C

    def d2a_dT2(self, T: float) -> float:
        """Second derivative of EOS attraction parameter d²a/dT²."""
        T_eff = self._safe_T(T)
        Tr = T_eff / T_C
        sqTr = math.sqrt(Tr)
        return self.ac * self.kappa * (1.0 + self.kappa) / (2.0 * Tr * sqTr) / T_C**2

    def viscosity(self, P: float, T: float) -> float:
        """Dynamic viscosity [Pa·s] using NIST-calibrated Chapman-Enskog + Muzny."""
        return viscosity_h2(T, P, eos=self)

    def cp_ideal(self, T: float) -> float:
        """Ideal-gas molar heat capacity Cp°(T) [J/(mol·K)]."""
        return cp_ideal_h2(T)


# ── Ideal-gas heat capacity ──────────────────────────────────────────────────
def cp_ideal_h2(T: float) -> float:
    """
    Ideal-gas molar heat capacity Cp°(T) [J/(mol·K)] — NIST Shomate equation.
    Valid from 298 to 1000 K; acceptable down to ~200 K, unreliable below ~100 K.
    """
    T_val = max(float(T), 50.0)
    t = T_val / 1000.0
    Cp = 33.066178 - 11.363417 * t + 11.432816 * t**2 - 2.772874 * t**3 - 0.158558 / max(t**2, 1e-6)
    return max(Cp, 14.0)


def cp_specific_ideal_h2(T: float) -> float:
    """Ideal-gas specific heat capacity cp° [J/(kg·K)]."""
    return cp_ideal_h2(T) / M_H2


# ── Real-gas heat capacities ─────────────────────────────────────────────────
def _real_heat_capacities(
    T: float, P: float, eos: PengRobinsonH2 | PengRobinsonMixture
) -> tuple[float, float]:
    """
    Real-gas molar (Cp, Cv) [J/(mol·K)] from the Peng-Robinson departure functions.

        Cv = Cv° + Cv_res,   Cv_res = T·a''/(2√2·b) · ln[(V + (1+√2)b) / (V + (1−√2)b)]
        Cp = Cv − T·(∂P/∂T)²_V / (∂P/∂V)_T

    Departures are evaluated at the untranslated PR molar volume V_PR = Vm + c;
    a constant Peneloux shift leaves both capacities unchanged.
    """
    T_val = max(float(T), 50.0)
    P_val = max(float(P), 1.0)
    Vm = eos.molar_volume(P_val, T_val) + eos.volume_translation
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
    Cp_minus_Cv = -T_val * dPdT_V**2 / dPdV_T if dPdV_T < -1e-30 else R

    s2 = math.sqrt(2.0)
    num = Vm + (1.0 + s2) * b
    den = Vm + (1.0 - s2) * b
    Cv_res = (
        T_val * eos.d2a_dT2(T_val) / (2.0 * s2 * b) * math.log(num / den)
        if abs(den) > 1e-20 and num / den > 0.0
        else 0.0
    )

    cp0 = eos.cp_ideal(T_val) if hasattr(eos, "cp_ideal") else cp_ideal_h2(T_val)
    Cv = max(cp0 - R + Cv_res, 8.0)
    Cp = max(Cv + Cp_minus_Cv, 14.0)
    return Cp, Cv


def gamma_real_h2(
    T: float, P: float, eos: PengRobinsonH2 | PengRobinsonMixture | None = None
) -> float:
    """
    Real-gas heat capacity ratio γ = Cp / Cv at temperature T [K] and pressure P [Pa].
    Combines Shomate ideal capacity with PR-EOS departure functions.
    """
    if eos is None:
        eos = PengRobinsonH2()
    Cp, Cv = _real_heat_capacities(T, P, eos)
    return float(min(max(Cp / Cv, 1.10), 1.80))


def cp_real_h2(
    T: float, P: float, eos: PengRobinsonH2 | PengRobinsonMixture | None = None
) -> float:
    """
    Real-gas isobaric molar heat capacity Cp [J/(mol·K)] at (T, P).

    Combines the NIST Shomate ideal-gas capacity with the Peng-Robinson
    residual Cv and the exact Cp − Cv relation.
    """
    if eos is None:
        eos = PengRobinsonH2()
    return float(_real_heat_capacities(T, P, eos)[0])


def cp_specific_real_h2(
    T: float, P: float, eos: PengRobinsonH2 | PengRobinsonMixture | None = None
) -> float:
    """Real-gas isobaric specific heat capacity cp [J/(kg·K)] at (T, P)."""
    M = eos.molar_mass() if isinstance(eos, PengRobinsonMixture) else M_H2
    return cp_real_h2(T, P, eos) / M


# ── Viscosity ────────────────────────────────────────────────────────────────
# Dense-fluid residual viscosity parameters for normal hydrogen (Muzny et al., 2013 / NIST).
# Δμ(ρ) = c1 * ρ + c2 * ρ² [Pa·s] with ρ in kg/m³
_MU_RESIDUAL_C1 = 3.464e-9  # Pa·s / (kg/m³)
_MU_RESIDUAL_C2 = 1.083e-9  # Pa·s / (kg/m³)²


def viscosity_h2(
    T: float,
    P: float,
    rho: float | None = None,
    eos: PengRobinsonH2 | PengRobinsonMixture | None = None,
) -> float:
    """
    Dynamic viscosity of pure hydrogen (or gas mixture if PengRobinsonMixture is provided) [Pa·s].

    Combines the dilute-gas Chapman-Enskog collision-integral formulation with the
    Muzny et al. (2013) / NIST density-dependent dense fluid residual correction:
        μ(T, P) = μ₀(T) + c₁·ρ + c₂·ρ²
    Matches NIST reference data (CoolProp HEOS) within <1.5% across 250–450 K and 1–350 bar.
    """
    if isinstance(eos, PengRobinsonMixture) and not eos.is_pure_h2():
        return eos.viscosity(P, T, rho=rho)
    T_val = max(float(T), 50.0)
    P_val = max(float(P), 1.0)
    M_g = M_H2 * 1000.0  # g/mol
    T_star = T_val / LJ_EPS_KB
    Omega_v = (
        1.16145 / T_star**0.14874
        + 0.52487 / math.exp(0.7732 * T_star)
        + 2.16178 / math.exp(2.43787 * T_star)
    )
    mu0 = 2.6693e-6 * math.sqrt(M_g * T_val) / (LJ_SIGMA**2 * max(Omega_v, 1e-6))

    if rho is None:
        calc_eos = eos if eos is not None else PengRobinsonH2()
        rho_val = calc_eos.density(P_val, T_val)
    else:
        rho_val = max(float(rho), 0.0)

    delta_mu = _MU_RESIDUAL_C1 * rho_val + _MU_RESIDUAL_C2 * (rho_val**2)
    mu = mu0 + delta_mu
    return float(min(max(mu, 6.0e-6), 3.0e-5))


# ── Joule-Thomson coefficient ────────────────────────────────────────────────
_JT_EFFECTIVE_C_H2: float = -8.0e-6  # m³/mol empirical JT translation for pure H2


def joule_thomson_coefficient(
    T: float, P: float, eos: PengRobinsonH2 | PengRobinsonMixture | None = None
) -> float:
    """
    Joule-Thomson throttling coefficient μ_JT = (∂T/∂P)_H [K/Pa].

    From the thermodynamic identity μ_JT = [T·(∂Vm/∂T)_P − Vm] / Cp with
    Vm = Z·R·T/P, where Z is the (volume-translated) compressibility factor
    returned by ``eos.z_factor``:

        μ_JT = (R·T² / P) · (∂Z/∂T)_P / Cp

    Pure H2 (``PengRobinsonH2`` or ``eos=None``): plain PR badly underpredicts
    the H2 second virial coefficient, so an *effective* JT translation
    ``c_JT = -8.0e-6`` m³/mol (twice the volumetric Peneloux value) is applied,
    i.e. μ_JT·Cp = (R·T²/P)·(∂Z_PR/∂T)_P + c_JT.  This is an empirical fit to
    CoolProp/Leachman (within ±6 % at ~100 bar; -46 %…+23 % over 10-300 bar at
    293-350 K), not a consequence of the volumetric translation.  Mixture EOS
    instances are evaluated as given, with their own composition and
    translation; cubic-PR mixture JT values are screening estimates only.

    Physical note: Below the inversion temperature (~202 K at 1 bar), μ_JT > 0.
    At UHS reservoir conditions (T > 250 K), μ_JT < 0 for H2-rich gas.
    """
    T_val = max(float(T), 50.0)
    P_val = max(float(P), 1.0)
    dT = max(T_val * 1e-4, 0.01)

    if isinstance(eos, PengRobinsonMixture):
        Z_hi = eos.z_factor(P_val, T_val + dT)
        Z_lo = eos.z_factor(P_val, T_val - dT)
        num = (R_UNIVERSAL * T_val**2 / P_val) * (Z_hi - Z_lo) / (2.0 * dT)
        Cp = cp_real_h2(T_val, P_val, eos)
        return float(num / max(Cp, 1.0))

    plain = PengRobinsonH2()
    Z_hi = plain.z_factor(P_val, T_val + dT)
    Z_lo = plain.z_factor(P_val, T_val - dT)
    num = (R_UNIVERSAL * T_val**2 / P_val) * (Z_hi - Z_lo) / (2.0 * dT) + _JT_EFFECTIVE_C_H2
    Cp = cp_real_h2(T_val, P_val, plain)
    return float(num / max(Cp, 1.0))


# ── Coupled Thermal Modeling ─────────────────────────────────────────────────
def adiabatic_temperature(
    T_old: float,
    P_old: float,
    P_new: float,
    eos: PengRobinsonH2 | PengRobinsonMixture | None = None,
) -> float:
    """
    Isentropic gas temperature change under rapid expansion or compression:
    T_new = T_old · (P_new / P_old)^n, with the real-gas isentropic exponent

        n = (∂lnT/∂lnP)_S = P·(∂Vm/∂T)_P / Cp = R·(Z + T·(∂Z/∂T)_P) / Cp

    evaluated at the initial state (n → (γ−1)/γ for an ideal gas).
    """
    if eos is None:
        eos = PengRobinsonH2()
    T_init = max(float(T_old), 50.0)
    P_init = max(float(P_old), 1.0)
    P_final = max(float(P_new), 1.0)
    if P_final == P_init:
        return T_init
    dT = max(T_init * 1e-4, 0.01)
    Z = eos.z_factor(P_init, T_init)
    dZ_dT = (eos.z_factor(P_init, T_init + dT) - eos.z_factor(P_init, T_init - dT)) / (2.0 * dT)
    Cp = cp_real_h2(T_init, P_init, eos)
    n = R_UNIVERSAL * (Z + T_init * dZ_dT) / max(Cp, 1.0)
    n = min(max(n, 0.05), 0.5)
    return float(max(T_init * (P_final / P_init) ** n, 50.0))


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
    eos: PengRobinsonH2 | PengRobinsonMixture | None = None,
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


# ── Multicomponent Gas Mixture Thermodynamics (H2 / CH4 / CO2 / N2) ──────────
M_AIR: float = 28.9647e-3  # kg/mol dry air molar mass
T_STD_METRIC: float = 288.15  # K (15 °C metric standard temperature)
P_STD_METRIC: float = 101325.0  # Pa (1.01325 bar metric standard pressure)
V_MOLAR_STD: float = R_UNIVERSAL * T_STD_METRIC / P_STD_METRIC  # ~0.0236445 sm³/mol


@dataclass(frozen=True)
class GasComponentSpec:
    """Thermophysical specification for a pure subsurface gas component."""

    name: str
    molar_mass: float  # kg/mol
    T_c: float  # K critical temperature
    P_c: float  # Pa critical pressure
    omega: float  # Pitzer acentric factor
    c_peneloux: float  # m³/mol default Peneloux volume translation
    lj_sigma: float  # Å Lennard-Jones collision diameter
    lj_eps_kb: float  # K Lennard-Jones well depth (ε/kB)
    hhv_molar: float  # J/mol Higher Heating Value at 298.15 K
    lhv_molar: float  # J/mol Lower Heating Value at 298.15 K
    cp_ideal_298: float  # J/(mol·K) reference ideal-gas molar heat capacity
    mu_res_c1: float  # Pa·s / (kg/m³) linear dense-fluid viscosity residual
    mu_res_c2: float  # Pa·s / (kg/m³)² quadratic dense-fluid viscosity residual

    @property
    def ac(self) -> float:
        return 0.45724 * R_UNIVERSAL**2 * self.T_c**2 / self.P_c

    @property
    def b(self) -> float:
        return 0.07780 * R_UNIVERSAL * self.T_c / self.P_c

    @property
    def kappa(self) -> float:
        return 0.37464 + 1.54226 * self.omega - 0.26992 * self.omega**2

    def alpha(self, T: float) -> float:
        Tr = max(float(T), T_FLOOR) / self.T_c
        return (1.0 + self.kappa * (1.0 - math.sqrt(Tr))) ** 2

    def a_T(self, T: float) -> float:
        return self.ac * self.alpha(T)


# CH4/CO2/N2 c_peneloux and mu_res_* are fitted to CoolProp 8.0.0 (300-370 K,
# 20-250 bar; CO2 excludes T < 310 K above 70 bar): c minimises the maximum
# density error (CH4 1.2 %, CO2 4.2 %, N2 0.2 %), viscosity residuals are least
# squares on that density (max error CH4 1.0 %, CO2 3.2 %, N2 2.1 %).
GAS_COMPONENTS: dict[str, GasComponentSpec] = {
    "H2": GasComponentSpec(
        name="H2",
        molar_mass=M_H2,
        T_c=T_C,
        P_c=P_C,
        omega=OMEGA,
        c_peneloux=-4.0e-6,
        lj_sigma=LJ_SIGMA,
        lj_eps_kb=LJ_EPS_KB,
        hhv_molar=285830.0,
        lhv_molar=241820.0,
        cp_ideal_298=28.84,
        mu_res_c1=_MU_RESIDUAL_C1,
        mu_res_c2=_MU_RESIDUAL_C2,
    ),
    "CH4": GasComponentSpec(
        name="CH4",
        molar_mass=16.0425e-3,
        T_c=190.564,
        P_c=4.5992e6,
        omega=0.01142,
        c_peneloux=-3.3e-6,
        lj_sigma=3.758,
        lj_eps_kb=148.6,
        hhv_molar=890580.0,
        lhv_molar=802340.0,
        cp_ideal_298=35.69,
        mu_res_c1=2.0627e-8,
        mu_res_c2=1.9903e-10,
    ),
    "CO2": GasComponentSpec(
        name="CO2",
        molar_mass=44.0095e-3,
        T_c=304.128,
        P_c=7.3773e6,
        omega=0.22394,
        c_peneloux=2.1e-6,
        lj_sigma=3.941,
        lj_eps_kb=195.2,
        hhv_molar=0.0,
        lhv_molar=0.0,
        cp_ideal_298=37.11,
        mu_res_c1=2.2277e-10,
        mu_res_c2=8.4037e-11,
    ),
    "N2": GasComponentSpec(
        name="N2",
        molar_mass=28.0134e-3,
        T_c=126.192,
        P_c=3.3958e6,
        omega=0.03720,
        c_peneloux=-3.9e-6,
        lj_sigma=3.798,
        lj_eps_kb=71.4,
        hhv_molar=0.0,
        lhv_molar=0.0,
        cp_ideal_298=29.12,
        mu_res_c1=1.6912e-8,
        mu_res_c2=4.4211e-11,
    ),
}

# Symmetric binary interaction parameters k_ij (Qian et al. 2013; Hassanpouryouzband et al. 2020)
BINARY_K_IJ: dict[tuple[str, str], float] = {
    ("H2", "CH4"): 0.0156,
    ("CH4", "H2"): 0.0156,
    ("H2", "CO2"): -0.0500,
    ("CO2", "H2"): -0.0500,
    ("H2", "N2"): 0.1030,
    ("N2", "H2"): 0.1030,
    ("CH4", "CO2"): 0.0919,
    ("CO2", "CH4"): 0.0919,
    ("CH4", "N2"): 0.0311,
    ("N2", "CH4"): 0.0311,
    ("CO2", "N2"): -0.0170,
    ("N2", "CO2"): -0.0170,
}


def normalize_composition(comp: dict[str, float] | None) -> dict[str, float]:
    """Normalize a mole-fraction dictionary over GAS_COMPONENTS to sum to 1.0."""
    if not comp:
        return {"H2": 1.0}
    cleaned: dict[str, float] = {}
    total = 0.0
    for k, v in comp.items():
        key = k.strip().upper()
        if key not in GAS_COMPONENTS:
            raise ValueError(f"Unknown gas component {k!r}. Supported: {', '.join(GAS_COMPONENTS)}")
        val = max(float(v), 0.0)
        if val > 0.0:
            cleaned[key] = cleaned.get(key, 0.0) + val
            total += val
    if total <= 0.0:
        return {"H2": 1.0}
    return {k: v / total for k, v in cleaned.items()}


class PengRobinsonMixture:
    """
    Multicomponent Peng-Robinson Equation of State for subsurface gas mixtures
    (H2, CH4, CO2, N2) using van der Waals one-fluid mixing rules, Peneloux
    molar volume translation, and Wilke + dense-fluid residual viscosity.

    Fully duck-type compatible with :class:`PengRobinsonH2` so any 0D, 1D, 2D,
    or 3D solver can accept a ``PengRobinsonMixture`` instance directly.
    """

    def __init__(
        self,
        composition: dict[str, float] | None = None,
        peneloux: bool = False,
        volume_translation: float | None = None,
        k_ij_overrides: dict[tuple[str, str], float] | None = None,
        max_iter: int = 50,
        tol: float = 1e-10,
    ) -> None:
        self.max_iter = max_iter
        self.tol = tol
        self.peneloux = bool(peneloux)
        self._explicit_vt = volume_translation
        self._k_ij = dict(BINARY_K_IJ)
        if k_ij_overrides:
            for (c1, c2), val in k_ij_overrides.items():
                u1, u2 = c1.upper(), c2.upper()
                self._k_ij[(u1, u2)] = float(val)
                self._k_ij[(u2, u1)] = float(val)
        self.ac: float = 1.0
        self.composition: dict[str, float] = {"H2": 1.0}
        self.set_composition(composition)

    def set_composition(self, composition: dict[str, float] | None) -> None:
        """Update active mole-fraction composition (normalized to sum to 1.0)."""
        self.composition = normalize_composition(composition)

    def is_pure_h2(self, composition: dict[str, float] | None = None) -> bool:
        comp = normalize_composition(composition) if composition is not None else self.composition
        return len(comp) == 1 and comp.get("H2", 0.0) >= 1.0 - 1e-12

    def _get_comp(self, composition: dict[str, float] | None) -> dict[str, float]:
        return normalize_composition(composition) if composition is not None else self.composition

    @property
    def b(self) -> float:
        """Mixture covolume b_mix [m³/mol] for current composition."""
        return self._b_mix(self.composition)

    @property
    def volume_translation(self) -> float:
        """Mixture Peneloux volume translation c_mix [m³/mol]."""
        return self._c_mix(self.composition)

    def _b_mix(self, comp: dict[str, float]) -> float:
        return sum(y * GAS_COMPONENTS[sp].b for sp, y in comp.items())

    def _c_mix(self, comp: dict[str, float]) -> float:
        if self._explicit_vt is not None:
            return float(self._explicit_vt)
        if not self.peneloux:
            return 0.0
        return sum(y * GAS_COMPONENTS[sp].c_peneloux for sp, y in comp.items())

    def _a_mix(self, T: float, comp: dict[str, float]) -> float:
        T_eff = max(float(T), T_FLOOR)
        species = list(comp.items())
        a_pure = {sp: GAS_COMPONENTS[sp].a_T(T_eff) for sp, _ in species}
        a_val = 0.0
        for sp_i, y_i in species:
            for sp_j, y_j in species:
                k_ij = 0.0 if sp_i == sp_j else self._k_ij.get((sp_i, sp_j), 0.0)
                a_val += y_i * y_j * math.sqrt(a_pure[sp_i] * a_pure[sp_j]) * (1.0 - k_ij)
        return a_val

    def _alpha(self, T: float) -> float:
        """Compatibility hook for departure functions (ac * _alpha(T) == a_mix(T))."""
        return self._a_mix(T, self.composition)

    def d2a_dT2(self, T: float, composition: dict[str, float] | None = None) -> float:
        """Second temperature derivative of mixture attraction parameter d²a_mix/dT²."""
        comp = self._get_comp(composition)
        if len(comp) == 1 and "H2" in comp:
            T_eff = max(float(T), T_FLOOR)
            Tr = T_eff / T_C
            sqTr = math.sqrt(Tr)
            spec = GAS_COMPONENTS["H2"]
            return spec.ac * spec.kappa * (1.0 + spec.kappa) / (2.0 * Tr * sqTr) / T_C**2
        T_eff = max(float(T), T_FLOOR)
        dT = max(T_eff * 1e-4, 0.01)
        a_hi = self._a_mix(T_eff + dT, comp)
        a_0 = self._a_mix(T_eff, comp)
        a_lo = self._a_mix(max(T_eff - dT, T_FLOOR), comp)
        return (a_hi - 2.0 * a_0 + a_lo) / (dT**2)

    def molar_mass(self, composition: dict[str, float] | None = None) -> float:
        """Mixture molar mass M_mix [kg/mol]."""
        comp = self._get_comp(composition)
        return sum(y * GAS_COMPONENTS[sp].molar_mass for sp, y in comp.items())

    def mass_fractions(self, composition: dict[str, float] | None = None) -> dict[str, float]:
        """Convert mole fractions to mass fractions w_i [-]."""
        comp = self._get_comp(composition)
        M_mix = self.molar_mass(comp)
        return {sp: y * GAS_COMPONENTS[sp].molar_mass / M_mix for sp, y in comp.items()}

    @staticmethod
    def mole_fractions_from_mass(mass_comp: dict[str, float]) -> dict[str, float]:
        """Convert mass fractions (or component masses [kg]) to normalized mole fractions y_i [-]."""
        moles: dict[str, float] = {}
        for sp, m_val in mass_comp.items():
            key = sp.strip().upper()
            if key not in GAS_COMPONENTS:
                raise ValueError(f"Unknown gas component {sp!r}")
            if m_val > 0.0:
                moles[key] = moles.get(key, 0.0) + float(m_val) / GAS_COMPONENTS[key].molar_mass
        return normalize_composition(moles)

    def cp_ideal(self, T: float, composition: dict[str, float] | None = None) -> float:
        """Mixture ideal-gas molar heat capacity Cp°_mix(T) [J/(mol·K)]."""
        comp = self._get_comp(composition)
        T_val = max(float(T), T_FLOOR)
        cp_val = 0.0
        for sp, y in comp.items():
            if sp == "H2":
                cp_val += y * cp_ideal_h2(T_val)
            else:
                spec = GAS_COMPONENTS[sp]
                # Smooth temperature scaling around 298.15 K for CH4/CO2/N2
                cp_sp = spec.cp_ideal_298 * (1.0 + 0.00045 * (T_val - 298.15))
                cp_val += y * max(cp_sp, 20.0)
        return cp_val

    def z_factor(self, P: float, T: float, composition: dict[str, float] | None = None) -> float:
        """Mixture compressibility factor Z = P·Vm / (R·T)."""
        if P <= 0.0 or T <= 0.0:
            return 1.0
        comp = self._get_comp(composition)
        T_eff = max(float(T), T_FLOOR)
        a_mix = self._a_mix(T_eff, comp)
        b_mix = self._b_mix(comp)
        A = a_mix * P / (R_UNIVERSAL**2 * T_eff**2)
        B = b_mix * P / (R_UNIVERSAL * T_eff)
        if not (math.isfinite(A) and math.isfinite(B)):
            return 1.0

        c2 = -(1.0 - B)
        c1 = A - 3.0 * B**2 - 2.0 * B
        c0 = -(A * B - B**2 - B**3)

        roots = cubic_real_roots(c2, c1, c0)
        gas_roots = [r for r in roots if r > B + 1e-6 and r > 0.05 and math.isfinite(r)]
        Z = max(gas_roots) if gas_roots else max((r for r in roots if r > 0.0), default=1.0)

        for _ in range(self.max_iter):
            f = Z**3 + c2 * Z**2 + c1 * Z + c0
            df = 3.0 * Z**2 + 2.0 * c2 * Z + c1
            if abs(df) < 1e-20:
                break
            dZ = -f / df
            Z += dZ
            if abs(dZ) < self.tol:
                break

        c_mix = self._c_mix(comp)
        if c_mix != 0.0:
            Z = Z - c_mix * P / (R_UNIVERSAL * T_eff)

        return float(max(Z, 1e-4))

    def density(self, P: float, T: float, composition: dict[str, float] | None = None) -> float:
        """Mixture gas density [kg/m³]: ρ = P·M_mix / (Z·R·T)."""
        P_val = max(float(P), 1.0)
        T_val = max(float(T), T_FLOOR)
        comp = self._get_comp(composition)
        Z = self.z_factor(P_val, T_val, comp)
        M_mix = self.molar_mass(comp)
        return P_val * M_mix / (Z * R_UNIVERSAL * T_val)

    def molar_volume(
        self, P: float, T: float, composition: dict[str, float] | None = None
    ) -> float:
        """Mixture molar volume Vm [m³/mol]: Vm = Z·R·T / P."""
        P_val = max(float(P), 1.0)
        T_val = max(float(T), T_FLOOR)
        Z = self.z_factor(P_val, T_val, composition)
        return Z * R_UNIVERSAL * T_val / P_val

    def drho_dP(self, P: float, T: float, composition: dict[str, float] | None = None) -> float:
        """Isothermal derivative (∂ρ/∂P)_T,y [kg/(m³·Pa)]."""
        P_val = max(float(P), 1.0)
        T_val = max(float(T), T_FLOOR)
        comp = self._get_comp(composition)
        dP = max(P_val * 1.0e-6, 1.0)
        P_lo = max(P_val - dP, 1.0)
        P_hi = P_val + dP
        return (self.density(P_hi, T_val, comp) - self.density(P_lo, T_val, comp)) / (P_hi - P_lo)

    def isothermal_compressibility(
        self, P: float, T: float, composition: dict[str, float] | None = None
    ) -> float:
        """Isothermal gas compressibility cg = (1/ρ)·(∂ρ/∂P)_T [Pa⁻¹]."""
        P_val = max(float(P), 1.0)
        T_val = max(float(T), T_FLOOR)
        comp = self._get_comp(composition)
        rho = self.density(P_val, T_val, comp)
        return self.drho_dP(P_val, T_val, comp) / max(rho, 1e-12)

    def component_fugacity_coefficients(
        self, P: float, T: float, composition: dict[str, float] | None = None
    ) -> dict[str, float]:
        """Component fugacity coefficients φ̂_i = f̂_i / (y_i·P) in the mixture."""
        P_val = max(float(P), 1.0)
        T_val = max(float(T), T_FLOOR)
        comp = self._get_comp(composition)
        Z = self.z_factor(P_val, T_val, comp)
        c_mix = self._c_mix(comp)
        Z_pr = Z + c_mix * P_val / (R_UNIVERSAL * T_val)

        a_mix = self._a_mix(T_val, comp)
        b_mix = self._b_mix(comp)
        A = a_mix * P_val / (R_UNIVERSAL**2 * T_val**2)
        B = b_mix * P_val / (R_UNIVERSAL * T_val)
        s2 = math.sqrt(2.0)

        species = list(comp.keys())
        a_pure = {sp: GAS_COMPONENTS[sp].a_T(T_val) for sp in species}
        log_term = math.log(max((Z_pr + (1.0 + s2) * B) / max(Z_pr + (1.0 - s2) * B, 1e-12), 1e-12))

        phi_hat: dict[str, float] = {}
        for sp_i in species:
            spec_i = GAS_COMPONENTS[sp_i]
            b_ratio = spec_i.b / max(b_mix, 1e-20)
            sum_aij = 0.0
            for sp_j, y_j in comp.items():
                k_ij = 0.0 if sp_i == sp_j else self._k_ij.get((sp_i, sp_j), 0.0)
                sum_aij += y_j * math.sqrt(a_pure[sp_i] * a_pure[sp_j]) * (1.0 - k_ij)
            ln_phi = (
                b_ratio * (Z_pr - 1.0)
                - math.log(max(Z_pr - B, 1e-12))
                - (A / (2.0 * s2 * max(B, 1e-12)))
                * (2.0 * sum_aij / max(a_mix, 1e-30) - b_ratio)
                * log_term
            )
            c_i = (
                float(self._explicit_vt)
                if self._explicit_vt is not None
                else (spec_i.c_peneloux if self.peneloux else 0.0)
            )
            if c_i != 0.0:
                ln_phi -= c_i * P_val / (R_UNIVERSAL * T_val)
            phi_hat[sp_i] = math.exp(max(-50.0, min(50.0, ln_phi)))
        return phi_hat

    def fugacity_coefficient(self, P: float, T: float) -> float:
        """Effective mixture fugacity coefficient (or H2 fugacity coefficient)."""
        phis = self.component_fugacity_coefficients(P, T)
        return phis.get("H2", next(iter(phis.values())))

    def fugacity(self, P: float, T: float) -> float:
        """Hydrogen partial fugacity f̂_H2 = y_H2 · P · φ̂_H2 [Pa] for adsorption."""
        P_val = max(float(P), 1.0)
        phis = self.component_fugacity_coefficients(P_val, T)
        y_h2 = self.composition.get("H2", 0.0)
        return y_h2 * P_val * phis.get("H2", 1.0)

    def viscosity(
        self,
        P: float,
        T: float,
        rho: float | None = None,
        composition: dict[str, float] | None = None,
    ) -> float:
        """
        Mixture dynamic viscosity [Pa·s] at (P, T) via Wilke (1950) dilute-gas
        mixing rule + mass-weighted dense-fluid residual correction.
        """
        T_val = max(float(T), T_FLOOR)
        P_val = max(float(P), 1.0)
        comp = self._get_comp(composition)

        species = list(comp.items())
        mu0_pure: dict[str, float] = {}
        for sp, _ in species:
            spec = GAS_COMPONENTS[sp]
            M_g = spec.molar_mass * 1000.0
            T_star = T_val / spec.lj_eps_kb
            omega_v = (
                1.16145 / T_star**0.14874
                + 0.52487 / math.exp(0.7732 * T_star)
                + 2.16178 / math.exp(2.43787 * T_star)
            )
            mu0_pure[sp] = (
                2.6693e-6 * math.sqrt(M_g * T_val) / (spec.lj_sigma**2 * max(omega_v, 1e-6))
            )

        # Wilke mixing rule for dilute gas mixture
        mu0_mix = 0.0
        for sp_i, y_i in species:
            M_i = GAS_COMPONENTS[sp_i].molar_mass
            mu_i = mu0_pure[sp_i]
            denom = 0.0
            for sp_j, y_j in species:
                if sp_i == sp_j:
                    phi_ij = 1.0
                else:
                    M_j = GAS_COMPONENTS[sp_j].molar_mass
                    mu_j = mu0_pure[sp_j]
                    num_ij = (1.0 + math.sqrt(mu_i / mu_j) * (M_j / M_i) ** 0.25) ** 2
                    phi_ij = num_ij / math.sqrt(8.0 * (1.0 + M_i / M_j))
                denom += y_j * phi_ij
            mu0_mix += y_i * mu_i / max(denom, 1e-12)

        rho_val = max(float(rho), 0.0) if rho is not None else self.density(P_val, T_val, comp)
        w_map = self.mass_fractions(comp)
        c1_mix = sum(w * GAS_COMPONENTS[sp].mu_res_c1 for sp, w in w_map.items())
        c2_mix = sum(w * GAS_COMPONENTS[sp].mu_res_c2 for sp, w in w_map.items())
        mu = mu0_mix + c1_mix * rho_val + c2_mix * (rho_val**2)
        return float(min(max(mu, 6.0e-6), 2.0e-4))

    def heating_values_volumetric(
        self, composition: dict[str, float] | None = None
    ) -> tuple[float, float]:
        """
        Return volumetric (HHV, LHV) in [MJ/sm³] at metric standard conditions
        (1.01325 bar, 15 °C).
        """
        comp = self._get_comp(composition)
        hhv_mol = sum(y * GAS_COMPONENTS[sp].hhv_molar for sp, y in comp.items())
        lhv_mol = sum(y * GAS_COMPONENTS[sp].lhv_molar for sp, y in comp.items())
        hhv_vol = (hhv_mol / V_MOLAR_STD) * 1.0e-6
        lhv_vol = (lhv_mol / V_MOLAR_STD) * 1.0e-6
        return float(hhv_vol), float(lhv_vol)

    def specific_gravity(self, composition: dict[str, float] | None = None) -> float:
        """Gas specific gravity relative to dry air G_s = M_mix / M_air [-]."""
        return self.molar_mass(composition) / M_AIR

    def wobbe_index(self, composition: dict[str, float] | None = None) -> float:
        """Upper Wobbe Index WI = HHV_vol / sqrt(G_s) [MJ/sm³]."""
        comp = self._get_comp(composition)
        hhv_vol, _ = self.heating_values_volumetric(comp)
        gs = max(self.specific_gravity(comp), 1e-12)
        return float(hhv_vol / math.sqrt(gs))


# ── Two-Phase Aquifer Rock-Fluid Physics (Direction 3) ───────────────────────


@dataclass
class RelativePermeabilityModel:
    """
    Two-phase gas-brine relative permeability, capillary pressure, and
    Buckley-Leverett fractional flow model for aquifer/depleted-reservoir UHS.

    Supports:
    - ``"corey"``: Power-law Corey (1954) relative permeability curves
        k_rg(S_g) = k_rg,max * S_ge^(n_g),   k_rw(S_w) = k_rw,max * S_we^(n_w)
    - ``"brooks_corey"``: Brooks-Corey (1964) pore-size distribution model with
        drainage capillary pressure P_c(S_g) = P_e * S_we^(-1/λ).

    When ``krg_max`` is ``None`` (default), the endpoint gas relative
    permeability at connate water saturation ``S_g = 1 - swr`` defaults to the
    water-wet reduction ``(1 - swr) ** n_g`` (which equals ``1.0`` for ``swr = 0``).
    """

    model: str = "corey"  # 'corey' | 'brooks_corey'
    swr: float = 0.0  # [-] connate / immobile residual water saturation S_wr in [0, 0.95)
    sgc: float = 0.0  # [-] critical / trapped gas saturation S_gc in [0, 1 - S_wr)
    krg_max: float | None = None  # [-] endpoint k_rg at S_g = 1 - S_wr
    krw_max: float = 1.0  # [-] endpoint k_rw at S_w = 1 - S_gc
    n_g: float = 2.0  # [-] Corey gas exponent (> 0)
    n_w: float = 3.0  # [-] Corey water exponent (> 0)
    brooks_corey_lambda: float = 2.0  # [-] Brooks-Corey pore-size distribution index λ
    entry_pressure_Pa: float = 5.0e3  # [Pa] capillary entry pressure P_e
    active_gas_saturation: float | None = None  # [-] operating gas saturation S_g (None -> 1 - swr)

    def __post_init__(self) -> None:
        self.model = self.model.strip().lower()
        if self.model not in ("corey", "brooks_corey"):
            raise ValueError(
                f"Unknown relative permeability model {self.model!r}. "
                "Choose 'corey' or 'brooks_corey'."
            )
        if not (0.0 <= self.swr < 0.95):
            raise ValueError(f"swr must be in [0.0, 0.95), got {self.swr}")
        if not (0.0 <= self.sgc < 1.0 - self.swr):
            raise ValueError(f"sgc must be in [0.0, 1 - swr), got {self.sgc}")
        if self.n_g <= 0.0 or self.n_w <= 0.0:
            raise ValueError("Corey exponents n_g and n_w must be positive")
        if self.brooks_corey_lambda <= 0.0:
            raise ValueError("brooks_corey_lambda must be positive")

    @property
    def endpoint_krg(self) -> float:
        """Effective endpoint k_rg,max at S_g = 1 - S_wr."""
        if self.krg_max is not None:
            return min(max(float(self.krg_max), 0.0), 1.0)
        return float((1.0 - self.swr) ** self.n_g)

    def effective_gas_saturation(self, S_g: float) -> float:
        """Normalized mobile gas saturation S_ge in [0, 1]."""
        denom = max(1.0 - self.swr - self.sgc, 1e-12)
        return min(max((float(S_g) - self.sgc) / denom, 0.0), 1.0)

    def effective_water_saturation(self, S_w: float) -> float:
        """Normalized mobile water saturation S_we in [0, 1]."""
        denom = max(1.0 - self.swr - self.sgc, 1e-12)
        return min(max((float(S_w) - self.swr) / denom, 0.0), 1.0)

    def krg(self, S_g: float | None = None) -> float:
        """Gas relative permeability k_rg(S_g) in [0, 1]."""
        sg_val = (
            float(S_g)
            if S_g is not None
            else (
                float(self.active_gas_saturation)
                if self.active_gas_saturation is not None
                else (1.0 - self.swr)
            )
        )
        s_ge = self.effective_gas_saturation(sg_val)
        if s_ge <= 0.0:
            return 0.0
        k_max = self.endpoint_krg
        if self.model == "brooks_corey":
            s_we = 1.0 - s_ge
            lam = self.brooks_corey_lambda
            val = k_max * ((1.0 - s_we) ** 2) * (1.0 - s_we ** ((2.0 + lam) / lam))
            return float(min(max(val, 0.0), 1.0))
        return float(min(max(k_max * (s_ge**self.n_g), 0.0), 1.0))

    def krw(self, S_w: float | None = None) -> float:
        """Water/brine relative permeability k_rw(S_w) in [0, 1]."""
        sw_val = float(S_w) if S_w is not None else self.swr
        s_we = self.effective_water_saturation(sw_val)
        if s_we <= 0.0:
            return 0.0
        if self.model == "brooks_corey":
            lam = self.brooks_corey_lambda
            val = self.krw_max * (s_we ** ((2.0 + 3.0 * lam) / lam))
            return float(min(max(val, 0.0), 1.0))
        return float(min(max(self.krw_max * (s_we**self.n_w), 0.0), 1.0))

    def capillary_pressure_Pa(self, S_g: float | None = None) -> float:
        """Brooks-Corey drainage capillary pressure P_c(S_g) = P_g - P_w [Pa]."""
        sg_val = float(S_g) if S_g is not None else (1.0 - self.swr)
        s_we = max(1.0 - self.effective_gas_saturation(sg_val), 0.01)
        return float(self.entry_pressure_Pa * (s_we ** (-1.0 / self.brooks_corey_lambda)))

    def gas_saturation_from_capillary_pressure(self, P_c: float) -> float:
        """Inverse Brooks-Corey equilibrium gas saturation S_g(P_c) in [S_gc, 1 - S_wr] (0 below entry)."""
        if P_c <= self.entry_pressure_Pa or self.entry_pressure_Pa <= 0.0:
            return 0.0
        s_we = (self.entry_pressure_Pa / float(P_c)) ** self.brooks_corey_lambda
        s_we = min(max(s_we, 0.0), 1.0)
        return float(self.sgc + (1.0 - self.swr - self.sgc) * (1.0 - s_we))

    def fractional_flow_gas(
        self,
        S_g: float,
        mu_g: float = 9.0e-6,
        mu_w: float = 8.0e-4,
    ) -> float:
        """Buckley-Leverett gas fractional flow f_g(S_g) [-]."""
        lam_g = self.krg(S_g) / max(float(mu_g), 1e-12)
        lam_w = self.krw(1.0 - float(S_g)) / max(float(mu_w), 1e-12)
        denom = lam_g + lam_w
        if denom <= 0.0:
            return 0.0
        return float(lam_g / denom)

    def dfg_dSg(
        self,
        S_g: float,
        mu_g: float = 9.0e-6,
        mu_w: float = 8.0e-4,
    ) -> float:
        """Derivative df_g/dS_g [-] for Buckley-Leverett frontal advance speed."""
        sg_min = self.sgc
        sg_max = 1.0 - self.swr
        if S_g <= sg_min or S_g >= sg_max:
            return 0.0
        ds = max((sg_max - sg_min) * 1.0e-5, 1.0e-6)
        s_lo = max(S_g - ds, sg_min)
        s_hi = min(S_g + ds, sg_max)
        if s_hi <= s_lo:
            return 0.0
        return (
            self.fractional_flow_gas(s_hi, mu_g, mu_w) - self.fractional_flow_gas(s_lo, mu_g, mu_w)
        ) / (s_hi - s_lo)

    def buckley_leverett_front(
        self,
        mu_g: float = 9.0e-6,
        mu_w: float = 8.0e-4,
        n_samples: int = 200,
    ) -> tuple[float, float]:
        """
        Welge tangent construction for hydrogen drainage into a brine aquifer
        (initial S_g = 0). Returns ``(S_g_front, dfg_dSg_front)``.
        """
        sg_max = 1.0 - self.swr
        if sg_max <= self.sgc + 1e-9:
            return self.sgc, 0.0
        best_sg = sg_max
        best_secant = 0.0
        n_pts = max(int(n_samples), 20)
        for idx in range(1, n_pts + 1):
            sg = self.sgc + (sg_max - self.sgc) * (idx / n_pts)
            fg = self.fractional_flow_gas(sg, mu_g, mu_w)
            secant = fg / max(sg, 1e-12)
            if secant >= best_secant:
                best_secant = secant
                best_sg = sg
        return float(best_sg), float(best_secant)

    def effective_porosity_factor(self) -> float:
        """Fraction of pore volume occupied by the gas phase (1 - S_wr)."""
        if self.active_gas_saturation is not None:
            return min(max(float(self.active_gas_saturation), 0.0), 1.0 - self.swr)
        return max(1.0 - self.swr, 0.0)

    def effective_permeability_factor(self, S_g: float | None = None) -> float:
        """Effective gas relative permeability multiplier k_rg(S_g) in [0, 1]."""
        return self.krg(S_g)


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
