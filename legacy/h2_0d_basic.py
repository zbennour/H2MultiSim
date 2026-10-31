"""
Hydrogen subsurface storage toy simulator (dependency-free, single file).

Features:
- Lumped (0D) reservoir model with pressure/temperature coupling
- Real-gas effects using Peng-Robinson EOS for pure H2
- Simple well inflow performance to model injection / withdrawal dynamics
- Geological Monte Carlo (porosity, permeability, thickness variability)
- Complete cycle analysis (cushion gas fill, injection/withdrawal cycles)

Everything is implemented with Python standard library only.
This is not a design tool - it's a didactic, easily hackable prototype.
"""

import math
import random
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional


# ---------------------------------------------------------------------------
# Thermodynamics: Peng-Robinson EOS for pure hydrogen
# ---------------------------------------------------------------------------


class H2PengRobinsonEOS:
    """
    Simple Peng-Robinson EOS implementation for pure hydrogen.

    Uses:
      P in Pa, T in K
      Returns compressibility factor Z (dimensionless) and density in kg/m3.
    """

    # Universal constants / critical properties for H2
    R_u = 8.314462618  # J/mol/K
    M = 2.016e-3       # kg/mol (molar mass of H2)
    R_spec = R_u / M   # J/kg/K

    # Hydrogen critical properties (approx)
    T_c = 33.2         # K
    P_c = 1.293e6      # Pa
    omega = -0.216     # Acentric factor

    def __init__(self, max_iter: int = 50, tol: float = 1e-6):
        self.max_iter = max_iter
        self.tol = tol

        # Peng-Robinson parameters (molar basis)
        self.a_c = 0.45724 * self.R_u ** 2 * self.T_c ** 2 / self.P_c
        self.b = 0.07780 * self.R_u * self.T_c / self.P_c

        self.kappa = (
            0.37464
            + 1.54226 * self.omega
            - 0.26992 * self.omega ** 2
        )

    def _alpha(self, T: float) -> float:
        Tr = T / self.T_c
        return (1.0 + self.kappa * (1.0 - math.sqrt(Tr))) ** 2

    def z_factor(self, P: float, T: float) -> float:
        """
        Solve cubic equation of state for Z with Newton-Raphson
        starting from ideal gas guess Z=1.0.

        For gas phase we take the largest real root, which Newton with
        Z0=1 usually finds directly at relevant conditions.
        """
        if P <= 0.0:
            return 1.0

        alpha = self._alpha(T)
        aT = self.a_c * alpha

        A = aT * P / (self.R_u ** 2 * T ** 2)
        B = self.b * P / (self.R_u * T)

        # Cubic in Z: Z^3 - (1-B)Z^2 + (A - 3B^2 - 2B)Z - (AB - B^2 - B^3) = 0
        def f(Z: float) -> float:
            return (
                Z ** 3
                - (1.0 - B) * Z ** 2
                + (A - 3.0 * B ** 2 - 2.0 * B) * Z
                - (A * B - B ** 2 - B ** 3)
            )

        def df(Z: float) -> float:
            return (
                3.0 * Z ** 2
                - 2.0 * (1.0 - B) * Z
                + (A - 3.0 * B ** 2 - 2.0 * B)
            )

        Z = 1.0  # initial guess (ideal gas)
        for _ in range(self.max_iter):
            F = f(Z)
            dF = df(Z)
            if abs(dF) < 1e-12:
                break
            dZ = -F / dF
            Z += dZ
            if abs(dZ) < self.tol:
                break

        # Ensure reasonable bounds
        if Z <= 0.0 or Z > 10.0:
            Z = 1.0
        return Z

    def density(self, P: float, T: float) -> float:
        """Return gas density [kg/m3] at P,T."""
        if T <= 0.0:
            raise ValueError("Temperature must be positive")
        Z = self.z_factor(P, T)
        rho = P / (Z * self.R_spec * T)
        return rho


# ---------------------------------------------------------------------------
# Reservoir properties / state
# ---------------------------------------------------------------------------


@dataclass
class ReservoirProperties:
    area: float               # m2
    thickness: float          # m
    porosity: float           # fraction
    total_compressibility: float  # 1/Pa (rock + brine, very approximate)
    permeability: float       # m2
    reference_pressure: float # Pa
    reference_temperature: float  # K
    well_radius: float        # m
    drainage_radius: float    # m
    gas_viscosity: float      # Pa.s (approx, used in PI)


@dataclass
class ReservoirState:
    pressure: float           # Pa
    temperature: float        # K
    mass: float               # kg (gas in pore volume)


@dataclass
class SchedulePeriod:
    name: str
    mode: str                 # "inject", "withdraw", "shut_in"
    duration_days: float
    bottomhole_pressure: float  # Pa (used for inject/withdraw)
    max_mass_rate: float        # kg/s (facility limit)


@dataclass
class SimulationResult:
    times: List[float] = field(default_factory=list)        # seconds
    pressures: List[float] = field(default_factory=list)    # Pa
    temperatures: List[float] = field(default_factory=list) # K
    masses: List[float] = field(default_factory=list)       # kg
    mass_rates: List[float] = field(default_factory=list)   # kg/s

    # Aggregated metrics (filled after simulation)
    metrics: Dict[str, float] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Reservoir model (0D lumped, single tank) with simple well model
# ---------------------------------------------------------------------------


class ReservoirModel:
    def __init__(
        self,
        props: ReservoirProperties,
        eos: Optional[H2PengRobinsonEOS] = None,
        gamma_gas: float = 1.4,
        thermal_relaxation_time_days: float = 30.0,
    ):
        self.props = props
        self.eos = eos or H2PengRobinsonEOS()
        self.gamma = gamma_gas

        self.volume_bulk = self.props.area * self.props.thickness  # m3
        self.volume_pore_ref = self.volume_bulk * self.props.porosity
        self.tau_T = thermal_relaxation_time_days * 86400.0

    # ---- Pore volume and mass helpers ----

    def pore_volume(self, pressure: float) -> float:
        """Pore volume as function of pressure (simple exponential compressibility)."""
        dp = pressure - self.props.reference_pressure
        Vp = self.volume_pore_ref * math.exp(self.props.total_compressibility * dp)
        return Vp

    def mass_from_state(self, pressure: float, temperature: float) -> float:
        """Mass of gas in the reservoir at given P,T (assuming full gas-filled pores)."""
        Vp = self.pore_volume(pressure)
        rho = self.eos.density(pressure, temperature)
        return rho * Vp

    def dm_dP_numeric(self, pressure: float, temperature: float, dP: float = 1e4) -> float:
        """Numeric derivative dm/dP at constant T."""
        P1 = max(pressure - dP, 1e3)
        P2 = pressure + dP
        m1 = self.mass_from_state(P1, temperature)
        m2 = self.mass_from_state(P2, temperature)
        return (m2 - m1) / (P2 - P1)

    # ---- Well / rate model ----

    def productivity_index(self) -> float:
        """
        Very simple radial-flow productivity index (volumetric).

        PI_vol [m3/s/Pa] ≈ 2π k h / (μ ln(re/rw))

        Returns:
            PI_mass [kg/s/Pa] = PI_vol * density (but here we return volumetric PI;
            we convert to mass later with current density).
        """
        k = self.props.permeability
        h = self.props.thickness
        mu = self.props.gas_viscosity
        re = self.props.drainage_radius
        rw = self.props.well_radius
        if re <= rw:
            re = rw * 10.0
        denom = math.log(re / rw)
        if denom <= 0.0:
            denom = 1.0
        PI_vol = 2.0 * math.pi * k * h / (mu * denom)
        return PI_vol

    def compute_mass_rate(
        self,
        state: ReservoirState,
        period: SchedulePeriod,
    ) -> float:
        """
        Compute mass rate [kg/s] for the current schedule period.

        Sign convention:
          positive -> injection into reservoir
          negative -> withdrawal / production
        """
        if period.mode == "shut_in":
            return 0.0

        PI_vol = self.productivity_index()
        rho = self.eos.density(state.pressure, state.temperature)

        bhp = period.bottomhole_pressure

        if period.mode == "inject":
            delta_p = max(bhp - state.pressure, 0.0)
            q_vol = PI_vol * delta_p
            q_mass = q_vol * rho
            q_mass = min(q_mass, max(0.0, period.max_mass_rate))
            return q_mass
        elif period.mode == "withdraw":
            delta_p = max(state.pressure - bhp, 0.0)
            q_vol = PI_vol * delta_p
            q_mass = q_vol * rho
            q_mass = -min(q_mass, max(0.0, period.max_mass_rate))
            return q_mass
        else:
            return 0.0

    # ---- Time stepping ----

    def step(
        self,
        state: ReservoirState,
        period: SchedulePeriod,
        dt: float,
    ) -> Tuple[ReservoirState, float]:
        """
        Advance reservoir state by dt seconds.

        Returns:
            new_state, q_mass (kg/s)
        """
        P = state.pressure
        T = state.temperature

        # (1) Compute current mass and rate
        m = self.mass_from_state(P, T)
        q_mass = self.compute_mass_rate(state, period)

        # (2) Material balance -> new pressure using linearization
        if abs(q_mass) < 1e-18:
            # No flow: relax T only
            new_P = P
        else:
            dm = q_mass * dt
            dm_dP = self.dm_dP_numeric(P, T)
            if abs(dm_dP) < 1e-12:
                # Fallback: assume gas-dominated compressibility ~ m / P
                dm_dP = m / max(P, 1e3)
            dP = dm / dm_dP
            new_P = max(P + dP, 1e3)

        # (3) Temperature update: adiabatic compression/expansion + thermal relaxation
        if new_P <= 0.0 or P <= 0.0:
            new_P = max(new_P, 1e3)
            P = max(P, 1e3)

        exponent = (self.gamma - 1.0) / self.gamma
        T_adiabatic = T * (new_P / P) ** exponent

        T_geo = self.props.reference_temperature
        if self.tau_T > 0.0:
            relax = math.exp(-dt / self.tau_T)
            new_T = T_geo + (T_adiabatic - T_geo) * relax
        else:
            new_T = T_adiabatic

        # (4) Update mass from new P,T (keeps internal consistency)
        new_m = self.mass_from_state(new_P, new_T)

        new_state = ReservoirState(pressure=new_P, temperature=new_T, mass=new_m)
        return new_state, q_mass

    def run_schedule(
        self,
        initial_state: ReservoirState,
        schedule: List[SchedulePeriod],
        dt_seconds: float,
        collect_time_series: bool = True,
    ) -> SimulationResult:
        """
        Run full schedule and return SimulationResult.
        """
        result = SimulationResult()
        state = initial_state
        t = 0.0

        max_pressure = state.pressure
        min_pressure = state.pressure
        max_mass = state.mass
        min_mass = state.mass

        total_injected = 0.0   # kg
        total_withdrawn = 0.0  # kg

        if collect_time_series:
            result.times.append(t)
            result.pressures.append(state.pressure)
            result.temperatures.append(state.temperature)
            result.masses.append(state.mass)
            result.mass_rates.append(0.0)

        for period in schedule:
            n_steps = int(period.duration_days * 86400.0 / dt_seconds)
            for _ in range(max(n_steps, 1)):
                new_state, q_mass = self.step(state, period, dt_seconds)
                t += dt_seconds

                if q_mass >= 0.0:
                    total_injected += q_mass * dt_seconds
                else:
                    total_withdrawn += -q_mass * dt_seconds

                max_pressure = max(max_pressure, new_state.pressure)
                min_pressure = min(min_pressure, new_state.pressure)
                max_mass = max(max_mass, new_state.mass)
                min_mass = min(min_mass, new_state.mass)

                state = new_state

                if collect_time_series:
                    result.times.append(t)
                    result.pressures.append(state.pressure)
                    result.temperatures.append(state.temperature)
                    result.masses.append(state.mass)
                    result.mass_rates.append(q_mass)

        # Working gas capacity = max - min mass
        working_gas = max_mass - min_mass

        result.metrics = {
            "max_pressure_Pa": max_pressure,
            "min_pressure_Pa": min_pressure,
            "max_mass_kg": max_mass,
            "min_mass_kg": min_mass,
            "working_gas_kg": working_gas,
            "total_injected_kg": total_injected,
            "total_withdrawn_kg": total_withdrawn,
            "cycle_gas_balance_kg": total_injected - total_withdrawn,
        }
        return result


# ---------------------------------------------------------------------------
# Monte Carlo driver for geological uncertainty
# ---------------------------------------------------------------------------


@dataclass
class MonteCarloConfig:
    n_realizations: int
    # Distributions are simple (mean, std) or (mean, log10_std) etc.
    porosity_mean: float
    porosity_std: float

    thickness_mean: float
    thickness_std: float

    log10_perm_mean: float   # permeability [m2] in log10 space
    log10_perm_std: float

    geothermal_T_mean: float  # K
    geothermal_T_std: float

    area_mean: float
    area_std: float


def sample_geology(base: ReservoirProperties, cfg: MonteCarloConfig, rng: random.Random) -> ReservoirProperties:
    """Sample a new ReservoirProperties from Monte Carlo distributions around a base case."""
    phi = max(0.01, min(0.35, rng.gauss(cfg.porosity_mean, cfg.porosity_std)))
    h = max(5.0, rng.gauss(cfg.thickness_mean, cfg.thickness_std))
    logk = rng.gauss(cfg.log10_perm_mean, cfg.log10_perm_std)
    k = 10.0 ** logk
    T_geo = max(250.0, rng.gauss(cfg.geothermal_T_mean, cfg.geothermal_T_std))
    area = max(1e5, rng.gauss(cfg.area_mean, cfg.area_std))

    return ReservoirProperties(
        area=area,
        thickness=h,
        porosity=phi,
        total_compressibility=base.total_compressibility,
        permeability=k,
        reference_pressure=base.reference_pressure,
        reference_temperature=T_geo,
        well_radius=base.well_radius,
        drainage_radius=base.drainage_radius,
        gas_viscosity=base.gas_viscosity,
    )


def percentile(sorted_values: List[float], p: float) -> float:
    """Simple percentile for sorted list (0-100)."""
    if not sorted_values:
        return float("nan")
    if p <= 0.0:
        return sorted_values[0]
    if p >= 100.0:
        return sorted_values[-1]
    k = (len(sorted_values) - 1) * p / 100.0
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return sorted_values[int(k)]
    d0 = sorted_values[f] * (c - k)
    d1 = sorted_values[c] * (k - f)
    return d0 + d1


def run_monte_carlo(
    base_props: ReservoirProperties,
    initial_pressure: float,
    initial_temperature: float,
    schedule: List[SchedulePeriod],
    dt_seconds: float,
    mc_cfg: MonteCarloConfig,
    seed: int = 1234,
) -> Dict[str, Dict[str, float]]:
    """
    Run Monte Carlo ensemble and return aggregated statistics for key metrics.

    Returns:
        stats[metric_name] = {
            "mean": ...,
            "min": ...,
            "max": ...,
            "p10": ...,
            "p50": ...,
            "p90": ...,
        }
    """
    rng = random.Random(seed)
    eos = H2PengRobinsonEOS()

    # Accumulators
    metric_values: Dict[str, List[float]] = {}

    for _ in range(mc_cfg.n_realizations):
        rp = sample_geology(base_props, mc_cfg, rng)
        model = ReservoirModel(rp, eos=eos)
        # Initialize mass consistent with initial P,T and this geology
        init_mass = model.mass_from_state(initial_pressure, rp.reference_temperature)
        state = ReservoirState(
            pressure=initial_pressure,
            temperature=rp.reference_temperature,
            mass=init_mass,
        )
        result = model.run_schedule(state, schedule, dt_seconds, collect_time_series=False)

        for k, v in result.metrics.items():
            metric_values.setdefault(k, []).append(v)

    # Compute statistics
    stats: Dict[str, Dict[str, float]] = {}
    for name, values in metric_values.items():
        values_sorted = sorted(values)
        mean_v = sum(values) / len(values_sorted)
        stats[name] = {
            "mean": mean_v,
            "min": values_sorted[0],
            "max": values_sorted[-1],
            "p10": percentile(values_sorted, 10.0),
            "p50": percentile(values_sorted, 50.0),
            "p90": percentile(values_sorted, 90.0),
        }

    return stats


# ---------------------------------------------------------------------------
# Example usage / basic CLI
# ---------------------------------------------------------------------------


def build_default_base_reservoir() -> ReservoirProperties:
    """
    Define a simple base-case storage aquifer.

    Tunable to your site; numbers here are purely illustrative.
    """
    return ReservoirProperties(
        area=5e6,              # m2 (~2.2 km radius)
        thickness=40.0,        # m
        porosity=0.18,         # fraction
        total_compressibility=5e-10,  # 1/Pa
        permeability=100e-15,  # 100 mD in SI
        reference_pressure=150e5,  # Pa (150 bar)
        reference_temperature=320.0,  # K (~47 C)
        well_radius=0.1,       # m
        drainage_radius=2000.0,  # m
        gas_viscosity=1e-5,    # Pa.s
    )


def build_example_schedule(ref_pressure: float) -> List[SchedulePeriod]:
    """
    Build a simple storage cycle schedule:

    - Cushion gas fill: 180 days injection
    - Three yearly cycles of injection (120 d), withdrawal (120 d), shut-in (125 d)
    """
    bar_to_Pa = 1e5

    # Safety limit: 1.8 * reference pressure as rough cap
    inject_bhp = min(ref_pressure * 1.8, ref_pressure + 100.0 * bar_to_Pa)
    withdraw_bhp = max(ref_pressure * 0.6, ref_pressure - 80.0 * bar_to_Pa)

    max_mass_rate = 50.0  # kg/s (facility constraint)

    schedule: List[SchedulePeriod] = []

    # Cushion gas fill
    schedule.append(
        SchedulePeriod(
            name="cushion_fill",
            mode="inject",
            duration_days=180.0,
            bottomhole_pressure=inject_bhp,
            max_mass_rate=max_mass_rate,
        )
    )

    # Operation cycles (e.g., 3 years)
    for i in range(3):
        schedule.append(
            SchedulePeriod(
                name=f"cycle{i+1}_inject",
                mode="inject",
                duration_days=120.0,
                bottomhole_pressure=inject_bhp,
                max_mass_rate=max_mass_rate,
            )
        )
        schedule.append(
            SchedulePeriod(
                name=f"cycle{i+1}_withdraw",
                mode="withdraw",
                duration_days=120.0,
                bottomhole_pressure=withdraw_bhp,
                max_mass_rate=max_mass_rate,
            )
        )
        schedule.append(
            SchedulePeriod(
                name=f"cycle{i+1}_shut_in",
                mode="shut_in",
                duration_days=125.0,
                bottomhole_pressure=0.0,
                max_mass_rate=0.0,
            )
        )

    return schedule


def run_single_realization_demo():
    base = build_default_base_reservoir()
    eos = H2PengRobinsonEOS()
    model = ReservoirModel(base, eos=eos)

    # Initial state at hydrostatic conditions
    P0 = base.reference_pressure
    T0 = base.reference_temperature
    m0 = model.mass_from_state(P0, T0)

    initial_state = ReservoirState(pressure=P0, temperature=T0, mass=m0)
    schedule = build_example_schedule(P0)

    # 1-day time step
    dt = 86400.0

    result = model.run_schedule(initial_state, schedule, dt_seconds=dt, collect_time_series=True)

    print("Single realization results:")
    for k, v in result.metrics.items():
        print(f"  {k}: {v: .3e}")

    # Print last state
    print(f"Final time [days]: {result.times[-1] / 86400.0: .2f}")
    print(f"Final pressure [bar]: {result.pressures[-1] / 1e5: .2f}")
    print(f"Final temperature [degC]: {result.temperatures[-1] - 273.15: .2f}")
    print(f"Final mass [kt]: {result.masses[-1] / 1e6: .3f}")


def run_monte_carlo_demo():
    base = build_default_base_reservoir()
    schedule = build_example_schedule(base.reference_pressure)

    mc_cfg = MonteCarloConfig(
        n_realizations=20,
        porosity_mean=base.porosity,
        porosity_std=0.03,
        thickness_mean=base.thickness,
        thickness_std=5.0,
        log10_perm_mean=math.log10(base.permeability),
        log10_perm_std=0.3,
        geothermal_T_mean=base.reference_temperature,
        geothermal_T_std=5.0,
        area_mean=base.area,
        area_std=1e6,
    )

    dt = 86400.0  # 1 day
    stats = run_monte_carlo(
        base_props=base,
        initial_pressure=base.reference_pressure,
        initial_temperature=base.reference_temperature,
        schedule=schedule,
        dt_seconds=dt,
        mc_cfg=mc_cfg,
        seed=42,
    )

    print("\nMonte Carlo ensemble statistics (20 realizations):")
    for name, st in stats.items():
        print(f"\nMetric: {name}")
        print(f"  mean: {st['mean']: .3e}")
        print(f"  min : {st['min']:  .3e}")
        print(f"  p10 : {st['p10']:  .3e}")
        print(f"  p50 : {st['p50']:  .3e}")
        print(f"  p90 : {st['p90']:  .3e}")
        print(f"  max : {st['max']:  .3e}")


if __name__ == "__main__":
    # Run a single deterministic realization
    run_single_realization_demo()

    # Run a small Monte Carlo ensemble