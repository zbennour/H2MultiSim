"""
H2Store-Unified  —  Hydrogen Subsurface Storage Simulator
          Notation aligned with Code B where equivalent variables exist
==========================================================
A fully self-contained hydrogen storage simulator.
No installation required beyond Python itself.

HOW TO RUN:
    1. Make sure Python 3.8 or newer is installed
       (check by typing:  python --version  in your terminal)
    2. Save this file anywhere on your computer
    3. Open a terminal / command prompt in that folder
    4. Type:  python h2store_v3_fast.py
    5. Results appear on screen and are saved as CSV files
       in a folder called  h2store_output/

WHAT IT MODELS:
    - Real-gas hydrogen behaviour (Peng-Robinson equation of state)
    - Pressure-dependent pore volume (rigid / elastic / breathing-MOF)
    - Hydrogen adsorption on sorbent materials (Langmuir / DSL isotherm)
    - Well injection and withdrawal performance (Darcy flow)
    - Complete operation schedules (cushion fill, cycles, shut-in)
    - Geological uncertainty (Monte Carlo over porosity, permeability etc.)

SCENARIOS INCLUDED:
    1. Teesside-style salt cavern         (UK, fast cycles)
    2. Okoroafor depleted reservoir       (porous media, seasonal)
    3. Rough-style depleted gas field     (large North Sea field)
    4. MOF cryogenic pressure vessel      (sorbent, 77 K)

DEPENDENCIES:
    None — uses only Python standard library (math, csv, os, random)
"""

# ── standard library only ────────────────────────────────────────────────────
import csv
import math
import os
import random
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Tuple


# ═══════════════════════════════════════════════════════════════════════════
# 1.  PHYSICAL CONSTANTS
# ═══════════════════════════════════════════════════════════════════════════

R_u    = 8.314462618     # J mol⁻¹ K⁻¹   universal gas constant
M      = 2.016e-3        # kg mol⁻¹       molar mass of hydrogen
R_spec = R_u / M   # J kg⁻¹ K⁻¹    specific gas constant for H2

# Hydrogen critical properties (NIST)
T_c      = 33.19         # K      critical temperature
P_c      = 1.2964e6      # Pa     critical pressure
omega    = -0.219        # —      acentric factor

# Minimum temperature the EOS is valid above
T_FLOOR  = T_c + 0.5  # 33.69 K


# ═══════════════════════════════════════════════════════════════════════════
# 2.  UNIT CONVERSION HELPERS
#     All functions take a value in SI and return the converted value.
#     Named clearly so non-coders can read them.
# ═══════════════════════════════════════════════════════════════════════════

def pa_to_bar(P):        return P / 1e5
def pa_to_MPa(P):        return P / 1e6
def pa_to_psi(P):        return P / 6894.757293
def K_to_C(T):return T - 273.15
def K_to_F(T): return K_to_C(T) * 9.0 / 5.0 + 32.0
def kg_to_tonnes(m):     return m / 1000.0
def kg_to_short_tons(m): return m / 907.18474
def kg_to_lb(m):         return m * 2.20462262185
def m3_to_cubic_feet(V): return V * 35.3147
def millidarcy_to_m2(k): return k * 9.869233e-16


# ═══════════════════════════════════════════════════════════════════════════
# 3.  PURE-PYTHON MATH UTILITIES
#     Replacements for numpy functions — no external libraries needed.
# ═══════════════════════════════════════════════════════════════════════════

def _cbrt(x: float) -> float:
    """Cube root that works for negative numbers (math.pow can't)."""
    if x >= 0.0:
        return x ** (1.0 / 3.0)
    return -((-x) ** (1.0 / 3.0))


def _cubic_real_roots(p: float, q: float, r: float) -> List[float]:
    """
    Solve monic cubic:  x³ + p·x² + q·x + r = 0

    Uses the analytical trigonometric / Cardano method.
    Returns a sorted list of all real roots.
    This approach never fails, never requires a matrix solver,
    and works with pure Python math functions only.

    Three cases:
      D > 0  →  one real root   (Cardano)
      D < 0  →  three real roots (trigonometric / Viète)
      D ≈ 0  →  repeated root
    """
    # Transform to depressed cubic  t³ + A·t + B = 0  via  x = t - p/3
    p3 = p / 3.0
    A  = q - p3 * p
    B  = r + p3 * (2.0 * p3 ** 2 - q)
    D  = (B / 2.0) ** 2 + (A / 3.0) ** 3    # discriminant

    if D > 1e-14:
        # One real root
        sq = math.sqrt(D)
        u  = _cbrt(-B / 2.0 + sq)
        v  = _cbrt(-B / 2.0 - sq)
        return [u + v - p3]

    elif D < -1e-14:
        # Three distinct real roots — trigonometric form
        m     = 2.0 * math.sqrt(max(-A / 3.0, 0.0))
        denom = A * m
        arg   = (3.0 * B / denom) if abs(denom) > 1e-30 else 0.0
        arg   = max(-1.0, min(1.0, arg))   # clamp for acos
        theta = math.acos(arg) / 3.0
        return sorted([
            m * math.cos(theta - 2.0 * math.pi * k / 3.0) - p3
            for k in range(3)
        ])

    else:
        # Repeated root
        u = _cbrt(-B / 2.0)
        return sorted([2.0 * u - p3, -u - p3])


def _bisect(func: Callable, lo: float, hi: float,
            xtol: float = 1.0, max_iter: int = 80) -> float:
    """
    Find root of func in [lo, hi] by bisection.
    Pure Python replacement for scipy.optimize.brentq.
    Reliable but slightly slower — fine for our use.
    """
    f_lo = func(lo)
    f_hi = func(hi)
    if f_lo * f_hi > 0.0:
        # Root not bracketed — return whichever end is closer to zero
        return lo if abs(f_lo) < abs(f_hi) else hi
    for _ in range(max_iter):
        mid   = (lo + hi) / 2.0
        f_mid = func(mid)
        if abs(f_mid) < 1e-10 or (hi - lo) < xtol:
            return mid
        if f_lo * f_mid < 0.0:
            hi   = mid
            f_hi = f_mid
        else:
            lo   = mid
            f_lo = f_mid
    return (lo + hi) / 2.0


def _box_muller(rng: random.Random) -> float:
    """
    Generate one standard normal sample using the Box-Muller transform.
    Pure Python replacement for numpy.random.normal.
    """
    u1 = max(rng.random(), 1e-15)
    u2 = rng.random()
    return math.sqrt(-2.0 * math.log(u1)) * math.cos(2.0 * math.pi * u2)


def _lognormal_sample(rng: random.Random,
                       mean: float, cv: float) -> float:
    """
    Sample from a lognormal distribution with given mean and
    coefficient of variation (cv = std / mean).
    Lognormal is used for geological parameters because they are
    always positive and naturally span orders of magnitude.
    """
    sigma = math.sqrt(math.log(1.0 + cv ** 2))
    mu    = math.log(mean) - 0.5 * sigma ** 2
    z     = _box_muller(rng)
    return math.exp(mu + sigma * z)


def _runge_kutta_4(func: Callable, t: float, y: List[float],
                    dt: float) -> List[float]:
    """
    One step of the classic 4th-order Runge-Kutta method.
    Pure Python replacement for scipy.integrate.solve_ivp.

    More accurate than Euler because it evaluates the derivative
    at four points within the step and takes a weighted average.

    func(t, y) → dy/dt  (list of derivatives)
    Returns new y after time step dt.
    """
    k1 = func(t,            y)
    k2 = func(t + dt / 2.0, [y[i] + dt / 2.0 * k1[i] for i in range(len(y))])
    k3 = func(t + dt / 2.0, [y[i] + dt / 2.0 * k2[i] for i in range(len(y))])
    k4 = func(t + dt,       [y[i] + dt        * k3[i] for i in range(len(y))])
    return [y[i] + dt / 6.0 * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i])
            for i in range(len(y))]


def _percentile(sorted_values: List[float], p: float) -> float:
    """
    Compute percentile p (0–100) from a pre-sorted list.
    Pure Python replacement for numpy.percentile.
    """
    if not sorted_values:
        return float('nan')
    p   = max(0.0, min(100.0, p))
    idx = (len(sorted_values) - 1) * p / 100.0
    lo  = int(idx)
    hi  = min(lo + 1, len(sorted_values) - 1)
    return sorted_values[lo] + (sorted_values[hi] - sorted_values[lo]) * (idx - lo)


# ═══════════════════════════════════════════════════════════════════════════
# 4.  PENG-ROBINSON EOS FOR PURE HYDROGEN
#
#     The Peng-Robinson equation of state corrects the ideal gas law
#     for real-gas behaviour at high pressures. For hydrogen at 700 bar,
#     the ideal gas law is wrong by ~40%; PR EOS is wrong by <2%.
# ═══════════════════════════════════════════════════════════════════════════

class PengRobinsonH2:
    """
    Peng-Robinson EOS for pure hydrogen.

    Provides:
      z_factor(P, T)  — compressibility factor Z  [-]
      density(P, T)   — real-gas density           [kg/m³]
      fugacity(P, T)  — thermodynamic fugacity      [Pa]
      viscosity(P, T) — dynamic viscosity           [Pa·s]

    All inputs: P in Pa, T in K.
    Uses analytical cubic solver — never crashes or diverges.
    """

    def __init__(self):
        kappa   = 0.37464 + 1.54226 * omega - 0.26992 * omega ** 2
        self.kappa = kappa
        self.a_c    = 0.45724 * R_u ** 2 * T_c ** 2 / P_c
        self.b      = 0.07780 * R_u * T_c / P_c

    # ── internal helpers ─────────────────────────────────────────────────

    def _safe_T(self, T: float) -> float:
        """Clamp temperature above H2 critical point to keep EOS valid."""
        return max(float(T), T_FLOOR)

    def _alpha(self, T: float) -> float:
        """Temperature correction factor for molecular attraction."""
        Tr = self._safe_T(T) / T_c
        return (1.0 + self.kappa * (1.0 - math.sqrt(Tr))) ** 2

    def _AB(self, P: float, T: float) -> Tuple[float, float]:
        """Dimensionless PR parameters A and B."""
        T = self._safe_T(T)
        aT = self.a_c * self._alpha(T)
        A = aT * P / (R_u * T) ** 2
        B = self.b * P / (R_u * T)
        return A, B

    # ── public interface ─────────────────────────────────────────────────

    def z_factor(self, P: float, T: float) -> float:
        """
        Compressibility factor Z.
        Z = 1 for ideal gas; Z ≠ 1 for real gas.
        At 700 bar, 300 K: Z ≈ 1.45 for hydrogen.
        """
        if P <= 0.0 or T <= 0.0:
            return 1.0
        A, B = self._AB(P, T)
        if not (math.isfinite(A) and math.isfinite(B)):
            return 1.0

        # PR cubic: Z³ + c2·Z² + c1·Z + c0 = 0
        c2 = -(1.0 - B)
        c1 =  A - 3.0 * B ** 2 - 2.0 * B
        c0 = -(A * B - B ** 2 - B ** 3)

        roots     = _cubic_real_roots(c2, c1, c0)
        gas_roots = [z for z in roots if z > B + 1e-10 and math.isfinite(z)]

        if not gas_roots:
            return max(1.0, B + 1e-6)    # safe fallback
        return max(gas_roots)             # gas phase = largest real root

    def density(self, P: float, T: float) -> float:
        """Real-gas density [kg/m³]."""
        P = max(float(P), 1.0)
        T = self._safe_T(T)
        Z = self.z_factor(P, T)
        return P / (Z * R_spec * T)

    def fugacity(self, P: float, T: float) -> float:
        """
        Fugacity [Pa] — the thermodynamically correct driving force
        for adsorption at high pressure. At low P, fugacity ≈ pressure.
        At 100 bar, fugacity may differ from pressure by 10–20%.
        """
        P  = max(float(P), 1.0)
        T  = self._safe_T(T)
        Z  = self.z_factor(P, T)
        A, B = self._AB(P, T)
        sq2 = math.sqrt(2.0)
        d1  = Z + (1.0 + sq2) * B
        d2  = Z + (1.0 - sq2) * B
        if d2 <= 0.0 or (Z - B) <= 0.0 or B <= 0.0:
            return P                      # ideal gas fallback
        ln_phi = ((Z - 1.0)
                  - math.log(Z - B)
                  - A / (2.0 * sq2 * B) * math.log(d1 / d2))
        ln_phi = max(-50.0, min(50.0, ln_phi))
        return P * math.exp(ln_phi)

    def viscosity(self, P: float, T: float) -> float:
        """
        Dynamic viscosity [Pa·s] computed from P and T.
        Uses a Chung-type correlation fitted for H2.
        Typical values: ~8×10⁻⁶ Pa·s at 300 K, 1 bar
                        ~2×10⁻⁵ Pa·s at 300 K, 500 bar
        """
        T   = self._safe_T(T)
        rho = self.density(P, T)
        mu0 = 8.76e-6 * (T / 300.0) ** 0.67    # dilute-gas viscosity
        b_e = 2.6e-5                             # Enskog volume parameter
        return mu0 * (1.0 + 0.8 * b_e * rho + 0.35 * (b_e * rho) ** 2)


# ═══════════════════════════════════════════════════════════════════════════
# 5.  PORE VOLUME MODEL
#
#     Pore volume is the space available for hydrogen inside the reservoir.
#     It changes slightly with pressure because the rock / vessel compresses.
# ═══════════════════════════════════════════════════════════════════════════

@dataclass
class PoreVolumeModel:
    """
    How pore volume changes with pressure.

    model = 'rigid'
        Volume is constant. Use for steel pressure vessels.
        V = V0

    model = 'elastic'
        Volume grows exponentially with pressure (rock compressibility).
        Use for salt caverns and depleted reservoirs.
        V = V0 × exp(ct × (P − P0))

    model = 'breathing'
        Volume jumps at a threshold pressure (gate-opening effect).
        Use for flexible Metal-Organic Framework (MOF) materials.
        V follows a smooth S-curve from V0 to V0 × V_ratio.

    Parameters
    ----------
    V0      : reference pore volume [m³]
    P0      : reference pressure [Pa]
    ct      : total compressibility [Pa⁻¹]  (elastic model)
    P_gate  : gate-opening pressure [Pa]     (breathing model)
    k_gate  : gate sharpness [Pa⁻¹]         (breathing model)
    V_ratio : open/closed volume ratio       (breathing model)
    """
    model:   str   = 'elastic'
    V0:      float = 1.0
    P0:      float = 1e5
    total_compressibility: float = 1e-9
    P_gate:  float = 5e6
    k_gate:  float = 1e-6
    V_ratio: float = 2.0

    def volume(self, P: float) -> float:
        """Return pore volume [m³] at pressure P [Pa]."""
        P = max(float(P), 1.0)
        if self.model == 'rigid':
            return self.V0
        if self.model == 'elastic':
            return max(self.V0 * math.exp(self.total_compressibility * (P - self.P0)), 1e-12)
        if self.model == 'breathing':
            # Sigmoid: smoothly transitions from V0 to V0*V_ratio
            f = 1.0 / (1.0 + math.exp(-self.k_gate * (P - self.P_gate)))
            return self.V0 + (self.V0 * self.V_ratio - self.V0) * f
        raise ValueError(f"Unknown pore model: {self.model!r}. "
                         f"Choose 'rigid', 'elastic', or 'breathing'.")


# ═══════════════════════════════════════════════════════════════════════════
# 6.  ADSORPTION ISOTHERM
#
#     Adsorption is when hydrogen molecules stick to the surface of a
#     porous material (like MOF or activated carbon). This extra storage
#     mechanism can significantly increase capacity at moderate pressures.
# ═══════════════════════════════════════════════════════════════════════════

@dataclass
class IsothermModel:
    """
    Adsorption isotherm model.

    model = 'none'
        No adsorption. Use this for caverns and rock reservoirs
        where there is no sorbent material.

    model = 'langmuir'
        Single-site Langmuir isotherm.
        Good for activated carbons and simple porous materials.
        n = n_max × K×f / (1 + K×f)

    model = 'dsl'
        Dual-Site Langmuir — two independent adsorption sites.
        More accurate for MOFs that have two types of binding sites.
        n = n1 + n2  (sum of two Langmuir terms)

    All isotherms use fugacity (not raw pressure) as the driving force,
    which is thermodynamically correct at high pressures.

    Parameters
    ----------
    n_max1, n_max2 : saturation capacity  [mol / kg_sorbent]
    K0_1,   K0_2   : affinity constant at T_ref  [Pa⁻¹]
    T_ref          : reference temperature  [K]
    dH1,    dH2    : adsorption enthalpy  [J/mol]  (negative = exothermic)
    """
    model:  str   = 'none'
    n_max1: float = 20.0
    K0_1:   float = 1e-7
    n_max2: float = 5.0
    K0_2:   float = 5e-8
    T_ref:  float = 298.0
    dH1:    float = -6000.0
    dH2:    float = -4000.0

    def _K1(self, T: float) -> float:
        """Temperature-dependent affinity constant, site 1 (van't Hoff)."""
        return self.K0_1 * math.exp(-self.dH1 / R_u * (1.0/T - 1.0/self.T_ref))

    def _K2(self, T: float) -> float:
        """Temperature-dependent affinity constant, site 2 (van't Hoff)."""
        return self.K0_2 * math.exp(-self.dH2 / R_u * (1.0/T - 1.0/self.T_ref))

    def loading(self, fugacity_Pa: float, T_K: float) -> float:
        """
        Adsorbed amount [mol / kg_sorbent] at given fugacity and temperature.
        Returns 0 if model is 'none'.
        """
        if self.model == 'none':
            return 0.0
        f = max(float(fugacity_Pa), 0.0)
        T = max(float(T_K), T_FLOOR)
        k1 = self._K1(T)
        n1 = self.n_max1 * k1 * f / (1.0 + k1 * f)
        if self.model == 'langmuir':
            return n1
        k2 = self._K2(T)
        n2 = self.n_max2 * k2 * f / (1.0 + k2 * f)
        return n1 + n2


# ═══════════════════════════════════════════════════════════════════════════
# 7.  RESERVOIR PROPERTIES  (user inputs go here)
# ═══════════════════════════════════════════════════════════════════════════

@dataclass
class ReservoirProperties:
    """
    All geological and engineering parameters for one storage site.

    In a real application, the user would fill these in from:
      - Seismic surveys      → area, thickness, depth
      - Core analysis        → porosity, S_wi
      - Well tests           → permeability (via WellParams)
      - Pressure gauges      → reference_pressure
      - Temperature logs     → T_surface, geothermal gradient
      - Regulatory limits    → P_min, P_max

    Parameters
    ----------
    area               : areal extent of the reservoir [m²]
    thickness          : net pay or cavern height [m]
    porosity           : fraction of rock that is void space [—]
    S_wi               : fraction of pore space occupied by water [—]
    depth              : depth to reservoir mid-point [m]
    reference_pressure : initial / starting pressure [Pa]
    P_min, P_max       : operating pressure limits [Pa]
    T_surface          : ground surface temperature [K]
    geothermal         : temperature increase per metre of depth [K/m]
    c_p_rock           : rock heat capacity [J/kg/K]
    rho_rock           : rock density [kg/m³]
    k_therm            : rock thermal conductivity [W/m/K]
    m_sorbent          : mass of adsorbing material in reservoir [kg]
                         (set to 0 for caverns and plain reservoirs)
    dH_ads             : heat released per mole adsorbed [J/mol]
    tau_T_days         : thermal relaxation time [days]
                         (how long it takes temperature to return
                          to geothermal equilibrium after a pressure event)
    """
    # Geometry
    area:               float
    thickness:          float
    porosity:           float
    S_wi:               float = 0.0

    # Location
    depth:              float = 0.0
    T_surface:          float = 288.0   # K  (= 15°C)
    geothermal:         float = 0.03    # K/m

    # Pressure limits
    reference_pressure: float = 1e5    # Pa
    P_min:              float = 1e5    # Pa
    P_max:              float = 2e7    # Pa

    # Rock thermal properties
    c_p_rock:           float = 900.0
    rho_rock:           float = 2500.0
    k_therm:            float = 2.5

    # Sorbent (set m_sorbent = 0 if no adsorbing material)
    m_sorbent:          float = 0.0
    dH_ads:             float = -6000.0

    # Thermal relaxation
    tau_T_days:         float = 30.0

    # ── computed properties ──────────────────────────────────────────────

    @property
    def V_bulk(self) -> float:
        """Total bulk volume of the reservoir [m³]."""
        return self.area * self.thickness

    @property
    def V_pore_0(self) -> float:
        """
        Available pore volume at reference conditions [m³].
        = bulk volume × porosity × (1 − water saturation)
        """
        return self.V_bulk * self.porosity * (1.0 - self.S_wi)

    @property
    def T_reservoir(self) -> float:
        """Geothermal temperature at reservoir depth [K]."""
        return self.T_surface + self.geothermal * self.depth

    @property
    def tau_T(self) -> float:
        """Thermal relaxation time in seconds."""
        return self.tau_T_days * 86400.0


# ═══════════════════════════════════════════════════════════════════════════
# 8.  WELL PARAMETERS
# ═══════════════════════════════════════════════════════════════════════════

@dataclass
class WellParams:
    """
    Well configuration for injection and withdrawal.

    Uses simple Darcy radial flow to compute the mass rate:
        q_mass = WI × (P_wellbore − P_reservoir) / viscosity × density

    Parameters
    ----------
    permeability    : rock permeability [millidarcy, mD]
                      Higher = easier for gas to flow
                      Typical: cavern ~5000 mD, sandstone 10–500 mD
    thickness       : perforated interval length [m]
    drainage_radius : distance gas must travel from reservoir to well [m]
    well_radius     : borehole radius [m]  (typically 0.1–0.3 m)
    skin            : damage factor [—]  (0 = undamaged, >0 = damaged)
    n_wells         : number of wells
    max_mass_rate   : facility limit on injection/withdrawal [kg/s]
                      (represents compressor / pipeline capacity)
    """
    permeability:    float = 50.0
    thickness:       float = 20.0
    drainage_radius: float = 500.0
    well_radius:     float = 0.15
    skin:            float = 0.0
    n_wells:         int   = 1
    max_mass_rate:   float = 1e9    # kg/s  (effectively unlimited by default)

    def well_index(self) -> float:
        """
        Well Index [m³ Pa⁻¹ s⁻¹] — how easily fluid flows into the well.
        Combines permeability, geometry, and number of wells.
        """
        k   = millidarcy_to_m2(self.permeability)
        re  = max(self.drainage_radius, self.well_radius * 10.0)
        den = math.log(re / self.well_radius) + self.skin
        den = max(den, 1e-10)
        return self.n_wells * k * self.thickness / den

    def mass_rate(self, P_reservoir: float, bhp: float,
                  T: float, eos: PengRobinsonH2) -> float:
        """
        Mass flow rate [kg/s].
        Positive = injection into reservoir.
        Negative = withdrawal / production from reservoir.
        Capped by max_mass_rate (facility limit).
        bhp = bottomhole pressure [Pa]  (Code B notation)
        """
        mu     = eos.viscosity(P_reservoir, T)
        rho    = eos.density(P_reservoir, T)
        PI_vol = self.well_index()           # productivity index [m³/s/Pa]
        delta_p = bhp - P_reservoir           # pressure drawdown [Pa]
        q_vol   = PI_vol * delta_p / max(mu, 1e-15)
        q_mass  = q_vol * rho
        cap     = self.max_mass_rate
        return max(-cap, min(cap, q_mass))


# ═══════════════════════════════════════════════════════════════════════════
# 9.  SCHEDULE SYSTEM
#
#     A schedule is a list of operation periods.
#     Each period specifies: what to do, how long, and at what pressure.
# ═══════════════════════════════════════════════════════════════════════════

@dataclass
class SchedulePeriod:
    """
    One operation period.

    name                : label (appears in output CSV)
    mode                : 'inject'   — push hydrogen in
                          'withdraw' — pull hydrogen out
                          'shut_in'  — no flow, let pressure equilibrate
    duration_days       : how long this period lasts [days]
    bottomhole_pressure : target pressure at the wellbore [Pa]
                          For inject: set higher than reservoir pressure
                          For withdraw: set lower than reservoir pressure
    max_mass_rate       : maximum flow rate for this period [kg/s]
    """
    name:                str
    mode:                str
    duration_days:       float
    bottomhole_pressure: float
    max_mass_rate:       float


# ═══════════════════════════════════════════════════════════════════════════
# 10.  SIMULATION RESULT CONTAINER
# ═══════════════════════════════════════════════════════════════════════════

@dataclass
class SimulationResult:
    """Stores all time-series output and summary metrics."""
    times:        List[float] = field(default_factory=list)   # seconds
    pressures:    List[float] = field(default_factory=list)   # Pa
    temperatures: List[float] = field(default_factory=list)   # K
    masses:       List[float] = field(default_factory=list)   # kg
    mass_rates:   List[float] = field(default_factory=list)   # kg/s
    phases:       List[str]   = field(default_factory=list)   # period name
    metrics:      Dict[str, float] = field(default_factory=dict)


# ═══════════════════════════════════════════════════════════════════════════
# 11.  RESERVOIR MODEL  (the main simulator)
# ═══════════════════════════════════════════════════════════════════════════

class ReservoirModel:
    """
    Lumped (0-D) hydrogen storage reservoir simulator.

    Integrates three coupled processes:
      1. Material balance  — mass of H2 = free gas + adsorbed gas
      2. Pressure update   — from mass change via EOS inversion
      3. Temperature update— adiabatic compression + geothermal relaxation
                             + adsorption/desorption heat

    Time integration uses 4th-order Runge-Kutta (pure Python).
    """

    def __init__(self,
                 props: ReservoirProperties,
                 pore:  PoreVolumeModel,
                 iso:   IsothermModel,
                 well:  WellParams,
                 eos:   Optional[PengRobinsonH2] = None,
                 gamma: float = 1.4):
        self.props = props
        self.pore  = pore
        self.iso   = iso
        self.well  = well
        self.eos   = eos or PengRobinsonH2()
        self.gamma = gamma   # adiabatic index for H2

    # ── mass helpers ─────────────────────────────────────────────────────

    def _adsorbed_moles(self, P: float, T: float) -> float:
        """Total adsorbed H2 [mol] in the reservoir."""
        if self.props.m_sorbent <= 0.0:
            return 0.0
        f = self.eos.fugacity(P, T)
        return self.iso.loading(f, T) * self.props.m_sorbent

    def total_mass(self, P: float, T: float) -> float:
        """
        Total H2 mass [kg] = free gas + adsorbed gas.
        This is what the material balance tracks.
        """
        Vp    = self.pore.volume(P)
        m_gas = self.eos.density(P, T) * Vp
        m_ads = self._adsorbed_moles(P, T) * M
        return m_gas + m_ads



    # ── temperature model ────────────────────────────────────────────────

    def _update_temperature(self, T: float, P_old: float, P_new: float,
                             dt: float, dn_ads_dt: float) -> float:
        """
        Three-part temperature update:

        1. Adiabatic step
           When pressure rises (injection), temperature rises like a
           compressed gas. When pressure falls (withdrawal), it cools.
           T_adiabatic = T × (P_new / P_old) ^ ((γ-1)/γ)

        2. Geothermal relaxation
           The reservoir rock slowly brings temperature back to the
           natural geothermal temperature over timescale tau_T.
           T_new = T_geo + (T_adiabatic - T_geo) × exp(-dt/tau_T)

        3. Adsorption heat
           Adsorption releases heat; desorption absorbs heat.
           The temperature change = Q_ads / total_heat_capacity
        """
        # 1. Adiabatic
        if P_old > 0.0 and P_new > 0.0 and math.isfinite(P_new / P_old):
            exponent    = (self.gamma - 1.0) / self.gamma
            T_adiabatic = T * (P_new / P_old) ** exponent
        else:
            T_adiabatic = T

        # 2. Geothermal relaxation
        T_geo = self.props.T_reservoir
        tau   = self.props.tau_T
        if tau > 0.0:
            relax = math.exp(-dt / tau)
            T_new = T_geo + (T_adiabatic - T_geo) * relax
        else:
            T_new = T_adiabatic

        # 3. Adsorption heat
        if self.props.m_sorbent > 0.0 and abs(dn_ads_dt) > 1e-15:
            m_now  = max(self.total_mass(P_new, T_new), 1.0)
            cp_H2  = 14300.0    # J kg⁻¹ K⁻¹
            m_rock = (self.props.V_bulk
                      * (1.0 - self.props.porosity)
                      * self.props.rho_rock)
            C_tot  = m_now * cp_H2 + m_rock * self.props.c_p_rock
            Q_ads  = dn_ads_dt * (-self.props.dH_ads)   # W
            T_new += Q_ads * dt / max(C_tot, 1.0)

        return max(T_new, T_FLOOR)

    # ── fast derivatives (P,T as state — no bisection needed) ────────────

    def _dm_dP(self, P: float, T: float, dP: float = 1e4) -> float:
        """
        Numerical derivative dm/dP at constant T.
        Uses two EOS calls — very fast, no iteration.
        This is how Code B achieves its speed: instead of inverting
        mass→pressure (slow bisection), we differentiate mass w.r.t.
        pressure (two direct EOS evaluations).
        """
        P1 = max(P - dP, 1.0)
        P2 = P + dP
        return (self.total_mass(P2, T) - self.total_mass(P1, T)) / (P2 - P1)

    def _derivatives(self, P: float, T: float,
                     P_wf: float) -> Tuple[float, float, float]:
        """
        Compute rates of change:  dP/dt, dT/dt, q_mass  [SI units/s]

        Strategy (fast — no bisection):
          1. q_mass  = well rate from Darcy IPR          (2 EOS calls)
          2. dP/dt   = q_mass / (dm/dP)   linearised MB  (2 EOS calls)
          3. dT/dt   = temperature model applied over dt_probe  (few ops)

        Total cost per call: ~6 EOS evaluations.
        Compare to old approach: ~160 EOS evaluations per step (bisection).
        """
        # Clamp state
        P = max(P, 1e3)
        T = max(T, T_FLOOR)

        # 1. Mass rate from well [kg/s]
        q_mass = self.well.mass_rate(P, P_wf, T, self.eos)  # P_wf = bhp

        # 2. Pressure rate via linearised material balance [Pa/s]
        #    dP/dt = (dm/dt) / (dm/dP)
        dm_dP = self._dm_dP(P, T)
        if abs(dm_dP) < 1e-15:
            # Fallback: ideal-gas estimate
            dm_dP = self.total_mass(P, T) / max(P, 1e3)
        dP_dt = q_mass / dm_dP

        # 3. Temperature rate [K/s]
        #    Project one small probe step to estimate dT/dt
        dt_probe = 600.0          # 10 minutes — small enough to be accurate
        P2 = max(P + dP_dt * dt_probe, 1e3)

        n_old  = self._adsorbed_moles(P,  T)
        n_new  = self._adsorbed_moles(P2, T)
        dn_dt  = (n_new - n_old) / dt_probe

        T2     = self._update_temperature(T, P, P2, dt_probe, dn_dt)
        dT_dt  = (T2 - T) / dt_probe

        return dP_dt, dT_dt, q_mass

    # ── run one schedule period  (fast Euler + RK4 hybrid) ───────────────

    def _run_period(self, P0: float, T0: float,
                    period: SchedulePeriod,
                    dt_s: float) -> Tuple[List, List, List, List, float, float]:
        """
        Integrate one SchedulePeriod.

        Time integration strategy:
          - State vector is (P, T) — pressure and temperature directly.
          - Derivatives dP/dt and dT/dt are computed analytically
            (no bisection, no matrix solve).
          - Uses RK4 for accuracy, but each RK4 step costs only
            ~24 EOS calls instead of ~640 in the old approach.
          - For shut-in periods: simple exponential T relaxation,
            no flow calculation needed at all.

        Speed vs Code B:
          - Code B Euler:   ~3  operations / step
          - This RK4:       ~24 EOS calls  / step  (≈ Code B × 8)
          - Old RK4+bisect: ~640 EOS calls / step  (≈ Code B × 200)
          So this is ~25× faster than the old unified code while
          still being more accurate than Code B's Euler method.
        """
        t_total = period.duration_days * 86400.0
        n_steps = max(int(t_total / dt_s), 1)
        dt_use  = t_total / n_steps

        # Determine wellbore pressure and flow cap
        if period.mode == 'shut_in':
            P_wf      = P0
            flow_cap  = 0.0
        else:
            P_wf      = period.bottomhole_pressure
            flow_cap  = min(period.max_mass_rate, self.well.max_mass_rate)

        original_cap          = self.well.max_mass_rate
        self.well.max_mass_rate = flow_cap

        P, T  = P0, T0
        t_now = 0.0

        ts = [0.0];  Ps = [P0];  Ts = [T0];  qs = [0.0]

        if period.mode == 'shut_in':
            # ── shut-in: pure thermal relaxation, no flow ─────────────────
            # Temperature decays exponentially — skip all EOS calls.
            T_geo = self.props.T_reservoir
            tau   = self.props.tau_T
            for _ in range(n_steps):
                t_now += dt_use
                if tau > 0.0:
                    T = T_geo + (T - T_geo) * math.exp(-dt_use / tau)
                T = max(T, T_FLOOR)
                ts.append(t_now); Ps.append(P)
                Ts.append(T);     qs.append(0.0)

        else:
            # ── injection / withdrawal: RK4 on (P, T) state ───────────────
            for _ in range(n_steps):

                # RK4 — four derivative evaluations per step
                dP1, dT1, q1 = self._derivatives(P,                T,                P_wf)
                dP2, dT2, _  = self._derivatives(P + dP1*dt_use/2, T + dT1*dt_use/2, P_wf)
                dP3, dT3, _  = self._derivatives(P + dP2*dt_use/2, T + dT2*dt_use/2, P_wf)
                dP4, dT4, _  = self._derivatives(P + dP3*dt_use,   T + dT3*dt_use,   P_wf)

                P = P + dt_use / 6.0 * (dP1 + 2*dP2 + 2*dP3 + dP4)
                T = T + dt_use / 6.0 * (dT1 + 2*dT2 + 2*dT3 + dT4)

                # Safety clamps
                P = max(min(P, self.props.P_max * 3.0), 1e3)
                T = max(T, T_FLOOR)

                t_now += dt_use
                ts.append(t_now); Ps.append(P)
                Ts.append(T);     qs.append(q1)

        self.well.max_mass_rate = original_cap
        return ts, Ps, Ts, qs, Ps[-1], Ts[-1]

    # ── main public method ────────────────────────────────────────────────

    def run_schedule(self,
                     schedule:       List[SchedulePeriod],
                     dt_days:        float = 1.0,
                     collect_series: bool  = True) -> SimulationResult:
        """
        Run the complete operation schedule.

        Parameters
        ----------
        schedule       : list of SchedulePeriod objects
        dt_days        : time step size in days (smaller = more accurate)
        collect_series : if True, store full time series; if False, only metrics

        Returns
        -------
        SimulationResult with time series and summary metrics
        """
        result   = SimulationResult()
        P        = self.props.reference_pressure
        T        = self.props.T_reservoir
        dt_s     = dt_days * 86400.0
        t_offset = 0.0

        max_P = P;    min_P = P
        max_m = 0.0;  min_m = 1e30
        total_inj = 0.0
        total_wdr = 0.0

        for period in schedule:
            ts, Ps, Ts, qs, P, T = self._run_period(P, T, period, dt_s)

            for i, (t, Pi, Ti, qi) in enumerate(zip(ts, Ps, Ts, qs)):
                mi    = self.total_mass(Pi, Ti)
                t_abs = t_offset + t
                min_P = min(min_P, Pi)
                max_m = max(max_m, mi)
                min_m = min(min_m, mi)

                dt_i = (ts[i] - ts[i-1]) if i > 0 else 0.0
                if qi > 0.0:
                    total_inj += qi * dt_i
                else:
                    total_wdr += -qi * dt_i

                if collect_series:
                    result.times.append(t_abs)
                    result.pressures.append(Pi)
                    result.temperatures.append(Ti)
                    result.masses.append(mi)
                    result.mass_rates.append(qi)
                    result.phases.append(period.name)

            t_offset += period.duration_days * 86400.0

        result.metrics = {
            'max_pressure_Pa':    max_P,
            'min_pressure_Pa':    min_P,
            'max_mass_kg':        max_m,
            'min_mass_kg':        min_m,
            'working_gas_kg':     max_m - min_m,
            'total_injected_kg':  total_inj,
            'total_withdrawn_kg': total_wdr,
            'cycle_balance_kg':   total_inj - total_wdr,
        }
        return result


# ═══════════════════════════════════════════════════════════════════════════
# 12.  SCHEDULE BUILDERS
# ═══════════════════════════════════════════════════════════════════════════

def schedule_salt_cavern(ref_pressure: float) -> List[SchedulePeriod]:
    """
    Salt cavern operation schedule.
    Fast, high-rate cycles typical of salt cavern storage.

    Phase 1: 60-day cushion gas fill
    Phase 2: 10 short cycles
              4 days inject → 3 days withdraw → 3 days shut-in
    """
    P_inj  = min(ref_pressure * 2.0,  ref_pressure + 50e5)
    P_prod = max(ref_pressure * 0.30, ref_pressure - 25e5)
    rate   = 200.0   # kg/s

    sched = [SchedulePeriod('cushion_fill', 'inject', 60.0, P_inj, rate)]
    for i in range(10):
        c = i + 1
        sched += [
            SchedulePeriod(f'c{c}_inject',   'inject',   4.0, P_inj,  rate),
            SchedulePeriod(f'c{c}_withdraw', 'withdraw', 3.0, P_prod, rate),
            SchedulePeriod(f'c{c}_shut_in',  'shut_in',  3.0, 0.0,   0.0),
        ]
    return sched


def schedule_porous_seasonal(ref_pressure: float) -> List[SchedulePeriod]:
    """
    Porous reservoir seasonal schedule.
    Inject in summer, withdraw in winter — 3 annual cycles.

    Phase 1: 365-day cushion gas fill
    Phase 2: 3 annual cycles
              180 days inject → 150 days withdraw → 35 days shut-in
    """
    P_inj  = min(ref_pressure * 1.5,  ref_pressure + 80e5)
    P_prod = max(ref_pressure * 0.50, ref_pressure - 60e5)
    rate   = 50.0

    sched = [SchedulePeriod('cushion_fill', 'inject', 365.0, P_inj, rate)]
    for i in range(3):
        c = i + 1
        sched += [
            SchedulePeriod(f'c{c}_inject',   'inject',   180.0, P_inj,  rate),
            SchedulePeriod(f'c{c}_withdraw', 'withdraw', 150.0, P_prod, rate),
            SchedulePeriod(f'c{c}_shut_in',  'shut_in',   35.0, 0.0,   0.0),
        ]
    return sched


def schedule_rough_field(ref_pressure: float) -> List[SchedulePeriod]:
    """
    Large depleted gas field schedule.
    Longer cycles suited to high-volume, lower-deliverability fields.

    Phase 1: 270-day cushion gas fill
    Phase 2: 3 cycles
              210 days inject → 120 days withdraw → 35 days shut-in
    """
    P_inj  = min(ref_pressure * 1.5, 250e5)
    P_prod = max(ref_pressure * 0.5,  80e5)
    rate   = 75.0

    sched = [SchedulePeriod('cushion_fill', 'inject', 270.0, P_inj, rate)]
    for i in range(3):
        c = i + 1
        sched += [
            SchedulePeriod(f'c{c}_inject',   'inject',   210.0, P_inj,  rate),
            SchedulePeriod(f'c{c}_withdraw', 'withdraw', 120.0, P_prod, rate),
            SchedulePeriod(f'c{c}_shut_in',  'shut_in',   35.0, 0.0,   0.0),
        ]
    return sched


def schedule_mof_vessel(ref_pressure: float) -> List[SchedulePeriod]:
    """
    MOF cryogenic vessel schedule.
    Short fill-hold-discharge cycles at 77 K.
    """
    return [
        SchedulePeriod('fill',      'inject',   5.0, 100e5, 5.0),
        SchedulePeriod('hold',      'shut_in',  1.0,   0.0, 0.0),
        SchedulePeriod('discharge', 'withdraw', 5.0,   1e5, 5.0),
        SchedulePeriod('hold2',     'shut_in',  1.0,   0.0, 0.0),
    ]


# ═══════════════════════════════════════════════════════════════════════════
# 13.  REAL-WORLD SCENARIOS
# ═══════════════════════════════════════════════════════════════════════════

def scenario_teesside_cavern() -> Dict:
    """
    Teesside-style salt cavern (UK).
    One solution-mined cavern ~70,000 m³, depth ~380 m, 20–100 bar.
    Reference: Teesside hydrogen storage project, UK.
    """
    props = ReservoirProperties(
        area=700.0, thickness=100.0, porosity=0.99, S_wi=0.0,
        depth=380.0, T_surface=283.0, geothermal=0.030,
        reference_pressure=40e5, P_min=20e5, P_max=100e5,
        m_sorbent=0.0, tau_T_days=15.0,
    )
    pore = PoreVolumeModel(model='elastic', V0=props.V_pore_0,
                            P0=40e5, total_compressibility=1e-10)
    iso  = IsothermModel(model='none')
    well = WellParams(permeability=5000.0, thickness=100.0,
                       drainage_radius=50.0, well_radius=0.20,
                       skin=0.0, n_wells=2, max_mass_rate=200.0)
    return dict(name='Teesside Salt Cavern',
                props=props, pore=pore, iso=iso, well=well,
                schedule_fn=schedule_salt_cavern)


def scenario_okoroafor_reservoir() -> Dict:
    """
    Depleted gas reservoir (Okoroafor et al. study parameters).
    3 km × 3 km × 100 m, porosity 0.20, k ~490 mD, 80 bar, ~43°C.
    Reference: Okoroafor et al. (2022), International Journal of
               Hydrogen Energy.
    """
    props = ReservoirProperties(
        area=9e6, thickness=100.0, porosity=0.20, S_wi=0.30,
        depth=1000.0, T_surface=288.0, geothermal=0.028,
        reference_pressure=80e5, P_min=40e5, P_max=160e5,
        m_sorbent=0.0, tau_T_days=45.0,
    )
    pore = PoreVolumeModel(model='elastic', V0=props.V_pore_0,
                            P0=80e5, total_compressibility=5e-10)
    iso  = IsothermModel(model='none')
    well = WellParams(permeability=490.0, thickness=100.0,
                       drainage_radius=1500.0, well_radius=0.10,
                       skin=0.0, n_wells=3, max_mass_rate=50.0)
    return dict(name='Okoroafor Depleted Reservoir',
                props=props, pore=pore, iso=iso, well=well,
                schedule_fn=schedule_porous_seasonal)


def scenario_rough_field() -> Dict:
    """
    Rough-style depleted gas field (UK North Sea).
    10 km × 3 km × 30 m, depth ~2700 m, k ~100 mD, ~150 bar.
    Reference: Rough gas field, offshore UK (now decommissioned UGS).
    """
    props = ReservoirProperties(
        area=30e6, thickness=30.0, porosity=0.15, S_wi=0.35,
        depth=2700.0, T_surface=288.0, geothermal=0.030,
        reference_pressure=150e5, P_min=80e5, P_max=250e5,
        m_sorbent=0.0, tau_T_days=60.0,
    )
    pore = PoreVolumeModel(model='elastic', V0=props.V_pore_0,
                            P0=150e5, total_compressibility=4e-10)
    iso  = IsothermModel(model='none')
    well = WellParams(permeability=100.0, thickness=30.0,
                       drainage_radius=2500.0, well_radius=0.15,
                       skin=1.0, n_wells=4, max_mass_rate=75.0)
    return dict(name='Rough-style Depleted Gas Field',
                props=props, pore=pore, iso=iso, well=well,
                schedule_fn=schedule_rough_field)


def scenario_mof_vessel() -> Dict:
    """
    MOF-packed cryogenic pressure vessel at 77 K.
    MOF-5 analogue, Langmuir isotherm, breathing pore model.
    Reference: Rowsell & Yaghi (2005), DOE hydrogen storage targets.
    """
    props = ReservoirProperties(
        area=0.5, thickness=1.0, porosity=0.40, S_wi=0.0,
        depth=0.0, T_surface=77.0, geothermal=0.0,
        reference_pressure=1e5, P_min=1e5, P_max=100e5,
        m_sorbent=200.0, dH_ads=-6500.0, tau_T_days=1.0,
    )
    pore = PoreVolumeModel(model='breathing', V0=props.V_pore_0,
                            P0=1e5, P_gate=20e5, k_gate=5e-5, V_ratio=2.2)
    iso  = IsothermModel(model='langmuir', n_max1=40.0,
                          K0_1=2e-7, T_ref=77.0, dH1=-6500.0)
    well = WellParams(permeability=5000.0, thickness=1.0,
                       drainage_radius=0.35, well_radius=0.01,
                       n_wells=1, max_mass_rate=5.0)
    return dict(name='MOF Cryogenic Vessel 77K',
                props=props, pore=pore, iso=iso, well=well,
                schedule_fn=schedule_mof_vessel)


# ═══════════════════════════════════════════════════════════════════════════
# 14.  MONTE CARLO
# ═══════════════════════════════════════════════════════════════════════════

@dataclass
class MonteCarloConfig:
    """
    Settings for geological uncertainty analysis.

    Each parameter is sampled many times from a probability distribution.
    The simulator runs once per sample, and the spread in results shows
    how sensitive storage capacity is to geological uncertainty.

    cv = coefficient of variation = standard deviation / mean
       cv = 0.10 means ±10% variability (low uncertainty)
       cv = 0.40 means ±40% variability (high uncertainty)

    Permeability is sampled in log10 space because it can vary over
    many orders of magnitude in nature (0.001 mD to 10,000 mD).
    """
    n_realizations: int   = 50

    porosity_mean:  float = 0.20;  porosity_cv:  float = 0.15
    thickness_mean: float = 20.0;  thickness_cv: float = 0.20
    log10_perm_mean: float = 1.0;  log10_perm_std: float = 0.35  # log10(mD)
    area_mean:      float = 1e6;   area_cv:      float = 0.10
    geothermal_T_mean: float = 330.0; geothermal_T_std: float = 5.0   # K


def run_monte_carlo(base_scenario: Dict,
                    mc_cfg:        MonteCarloConfig,
                    dt_days:       float = 1.0,
                    seed:          int   = 42) -> Dict[str, Dict[str, float]]:
    """
    Run Monte Carlo ensemble for geological uncertainty.

    For each realization:
      1. Sample porosity, thickness, permeability, area, temperature
      2. Build a perturbed reservoir model
      3. Run the full schedule
      4. Record summary metrics

    Returns statistics (mean, p10, p50, p90 etc.) for each metric.
    """
    rng    = random.Random(seed)
    values: Dict[str, List[float]] = {}
    base   = base_scenario['props']

    for i in range(mc_cfg.n_realizations):
        # Sample each uncertain parameter
        phi  = min(max(_lognormal_sample(rng, mc_cfg.porosity_mean,
                                          mc_cfg.porosity_cv), 0.01), 0.50)
        h    = max(_lognormal_sample(rng, mc_cfg.thickness_mean,
                                      mc_cfg.thickness_cv), 1.0)
        logk = rng.gauss(mc_cfg.log10_perm_mean, mc_cfg.log10_perm_std)
        k_mD = max(10.0 ** logk, 0.001)
        area = max(_lognormal_sample(rng, mc_cfg.area_mean,
                                      mc_cfg.area_cv), 1e4)
        T_geo = max(rng.gauss(mc_cfg.geothermal_T_mean, mc_cfg.geothermal_T_std), T_FLOOR)

        # Build perturbed reservoir
        props_i = ReservoirProperties(
            area=area, thickness=h, porosity=phi,
            S_wi=base.S_wi, depth=base.depth,
            T_surface=T_geo - base.geothermal * base.depth,
            geothermal=base.geothermal,
            reference_pressure=base.reference_pressure,
            P_min=base.P_min, P_max=base.P_max,
            m_sorbent=base.m_sorbent, dH_ads=base.dH_ads,
            tau_T_days=base.tau_T_days,
        )
        pore_i = PoreVolumeModel(
            model=base_scenario['pore'].model,
            V0=props_i.V_pore_0,
            P0=base_scenario['pore'].P0,
            total_compressibility=base_scenario['pore'].total_compressibility,
        )
        well_i = WellParams(
            permeability=k_mD,
            thickness=h,
            drainage_radius=base_scenario['well'].drainage_radius,
            well_radius=base_scenario['well'].well_radius,
            skin=base_scenario['well'].skin,
            n_wells=base_scenario['well'].n_wells,
            max_mass_rate=base_scenario['well'].max_mass_rate,
        )
        model_i = ReservoirModel(props_i, pore_i,
                                  base_scenario['iso'], well_i)
        sched_i = base_scenario['schedule_fn'](props_i.reference_pressure)

        try:
            res = model_i.run_schedule(sched_i, dt_days=dt_days,
                                        collect_series=False)
            for key, val in res.metrics.items():
                values.setdefault(key, []).append(val)
        except Exception:
            pass   # skip failed realizations silently

        # Progress indicator every 10 samples
        if (i + 1) % 10 == 0:
            print(f'    Monte Carlo: {i+1}/{mc_cfg.n_realizations} done')

    # Compute statistics
    stats: Dict[str, Dict[str, float]] = {}
    for name, vals in values.items():
        sv = sorted(vals)
        stats[name] = {
            'mean': sum(sv) / len(sv),
            'min':  sv[0],
            'p10':  _percentile(sv, 10.0),
            'p50':  _percentile(sv, 50.0),
            'p90':  _percentile(sv, 90.0),
            'max':  sv[-1],
        }
    return stats


# ═══════════════════════════════════════════════════════════════════════════
# 15.  CSV EXPORT
# ═══════════════════════════════════════════════════════════════════════════

def export_timeseries_csv(result: SimulationResult,
                           filepath: str,
                           label: str = '') -> None:
    """
    Export simulation time series to CSV.

    Columns include both SI units (Pa, K, kg) and imperial units
    (psi, °F, lb) so the file is useful to engineers from any background.
    The 'phase' column shows which operation period each row belongs to.
    """
    if not result.times:
        return

    folder = os.path.dirname(os.path.abspath(filepath))
    if folder:
        os.makedirs(folder, exist_ok=True)

    with open(filepath, 'w', newline='') as fh:
        fh.write(f'# H2Store-Unified output | {label}\n')
        fh.write('# SI units: Pa, K, kg, kg/s\n')
        fh.write('# Imperial: psi, degF, tonnes, short_tons, lb, lb/hr\n#\n')
        w = csv.writer(fh)
        w.writerow([
            'time_s', 'time_days', 'phase',
            'pressure_Pa', 'pressure_bar', 'pressure_MPa', 'pressure_psi',
            'temperature_K', 'temperature_C', 'temperature_F',
            'mass_kg', 'mass_tonnes', 'mass_short_tons', 'mass_lb',
            'rate_kg_s', 'rate_kg_day', 'rate_tonnes_day', 'rate_lb_hr',
        ])
        phases = result.phases or [''] * len(result.times)
        for t, P, T, m, q, ph in zip(
                result.times, result.pressures, result.temperatures,
                result.masses, result.mass_rates, phases):
            w.writerow([
                f'{t:.1f}',
                f'{t/86400:.5f}',
                ph,
                f'{P:.3e}',
                f'{pa_to_bar(P):.4f}',
                f'{pa_to_MPa(P):.5f}',
                f'{pa_to_psi(P):.4f}',
                f'{T:.3f}',
                f'{K_to_C(T):.3f}',
                f'{K_to_F(T):.3f}',
                f'{m:.6e}',
                f'{kg_to_tonnes(m):.6e}',
                f'{kg_to_short_tons(m):.6e}',
                f'{kg_to_lb(m):.6e}',
                f'{q:.6e}',
                f'{q*86400:.6e}',
                f'{kg_to_tonnes(q*86400):.6e}',
                f'{kg_to_lb(q)*3600:.6e}',
            ])

    print(f'    Saved: {os.path.abspath(filepath)}  ({len(result.times)} rows)')


def export_montecarlo_csv(stats: Dict,
                           filepath: str,
                           label: str = '') -> None:
    """Export Monte Carlo statistics table to CSV."""
    folder = os.path.dirname(os.path.abspath(filepath))
    if folder:
        os.makedirs(folder, exist_ok=True)
    with open(filepath, 'w', newline='') as fh:
        fh.write(f'# H2Store-Unified Monte Carlo | {label}\n#\n')
        w = csv.writer(fh)
        w.writerow(['metric', 'mean', 'min', 'p10', 'p50', 'p90', 'max'])
        for name, st in stats.items():
            w.writerow([name] + [
                f'{st[k]:.4e}'
                for k in ('mean', 'min', 'p10', 'p50', 'p90', 'max')
            ])
    print(f'    Saved: {os.path.abspath(filepath)}')


# ═══════════════════════════════════════════════════════════════════════════
# 16.  PRINT HELPERS
# ═══════════════════════════════════════════════════════════════════════════

def print_metrics(metrics: Dict[str, float], name: str) -> None:
    """Print simulation summary metrics in both SI and imperial."""
    print(f'\n  Summary — {name}')
    print(f'  {"─"*55}')
    for k, v in metrics.items():
        if 'pressure' in k:
            print(f'  {k:<30}  {v:>12.3e} Pa'
                  f'  ({pa_to_bar(v):.1f} bar,'
                  f'  {pa_to_psi(v):.1f} psi)')
        elif any(x in k for x in ('mass', 'gas', 'injected',
                                   'withdrawn', 'balance')):
            print(f'  {k:<30}  {v:>12.3e} kg'
                  f'  ({kg_to_tonnes(v):.2f} t,'
                  f'  {kg_to_short_tons(v):.2f} short ton)')
        else:
            print(f'  {k:<30}  {v:>12.3e}')


def print_mc_table(stats: Dict, label: str) -> None:
    """Print Monte Carlo statistics as a readable table."""
    print(f'\n  Monte Carlo Results — {label}')
    print(f'  {"─"*70}')
    print(f'  {"Metric":<28} {"mean":>11} {"p10":>11} '
          f'{"p50":>11} {"p90":>11}')
    print(f'  {"─"*70}')
    for name, st in stats.items():
        print(f'  {name:<28} {st["mean"]:>11.3e} {st["p10"]:>11.3e}'
              f' {st["p50"]:>11.3e} {st["p90"]:>11.3e}')


# ═══════════════════════════════════════════════════════════════════════════
# 17.  MAIN — runs when you type:  python h2store_v3_fast.py
# ═══════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':

    OUT = 'h2store_output'
    os.makedirs(OUT, exist_ok=True)

    print('=' * 65)
    print('  H2Store-Unified — Hydrogen Storage Simulator')
    print('  Pure Python — no external libraries required')
    print('=' * 65)
    print(f'  Output folder: ./{OUT}/')

    # List of all scenarios to run
    scenarios = [
        scenario_teesside_cavern(),
        scenario_okoroafor_reservoir(),
        scenario_rough_field(),
        scenario_mof_vessel(),
    ]

    for sc in scenarios:
        name = sc['name']
        slug = (name.lower()
                .replace(' ', '_')
                .replace('(', '').replace(')', '')
                .replace('/', '_'))

        print(f'\n{"═"*65}')
        print(f'  Scenario: {name}')
        print(f'{"═"*65}')

        # ── build and run ─────────────────────────────────────────────────
        model  = ReservoirModel(sc['props'], sc['pore'],
                                sc['iso'],   sc['well'])
        sched  = sc['schedule_fn'](sc['props'].reference_pressure)
        dt     = 0.5 if ('Cavern' in name or 'MOF' in name) else 1.0
        result = model.run_schedule(sched, dt_days=dt)

        # ── print summary ─────────────────────────────────────────────────
        print_metrics(result.metrics, name)

        if result.times:
            print(f'\n  Final state after'
                  f' {result.times[-1]/86400:.0f} days:')
            print(f'    Pressure    : '
                  f'{pa_to_bar(result.pressures[-1]):.2f} bar  '
                  f'({pa_to_psi(result.pressures[-1]):.1f} psi)')
            print(f'    Temperature : '
                  f'{K_to_C(result.temperatures[-1]):.2f} °C  '
                  f'({K_to_F(result.temperatures[-1]):.2f} °F)')
            print(f'    H2 mass     : '
                  f'{kg_to_tonnes(result.masses[-1]):.3f} t  '
                  f'({kg_to_short_tons(result.masses[-1]):.3f} short ton)')

        # ── export time series CSV ────────────────────────────────────────
        export_timeseries_csv(
            result,
            f'{OUT}/{slug}_timeseries.csv',
            label=name
        )

        # ── Monte Carlo (skip MOF — different physics, slow) ──────────────
        if 'MOF' not in name:
            print(f'\n  Running Monte Carlo'
                  f' ({50} realizations)...')
            mc_cfg = MonteCarloConfig(
                n_realizations = 50,
                porosity_mean  = sc['props'].porosity,
                porosity_cv    = 0.15,
                thickness_mean = sc['props'].thickness,
                thickness_cv   = 0.20,
                log10_perm_mean = math.log10(sc['well'].permeability),
                log10_perm_std  = 0.35,
                area_mean      = sc['props'].area,
                area_cv        = 0.10,
                geothermal_T_mean = sc['props'].T_reservoir,
                geothermal_T_std  = 5.0,
            )
            stats = run_monte_carlo(sc, mc_cfg, dt_days=dt)
            print_mc_table(stats, name)
            export_montecarlo_csv(
                stats,
                f'{OUT}/{slug}_montecarlo.csv',
                label=name
            )

    print(f'\n{"="*65}')
    print(f'  All done.  Check the ./{OUT}/ folder for CSV files.')
    print(f'{"="*65}')
