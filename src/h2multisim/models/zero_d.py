"""
h2multisim.models.zero_d
========================
Zero-dimensional (0D) lumped-parameter reservoir model for underground
hydrogen storage (UHS) in salt caverns, depleted gas fields, and aquifers.
"""

from __future__ import annotations

import csv
import math
import os
import random
from dataclasses import dataclass, field
from typing import Callable

from h2multisim.thermo import (
    IsothermModel,
    PengRobinsonH2,
    PoreVolumeModel,
    coupled_temperature_update,
)
from h2multisim.units import (
    K_to_C,
    K_to_F,
    bar_to_pa,
    kg_to_lb,
    kg_to_short_tons,
    kg_to_tonnes,
    mD_to_m2,
    pa_to_bar,
    pa_to_MPa,
    pa_to_psi,
    percentile,
)


@dataclass
class ReservoirProperties:
    """
    Geometric and rock properties for the 0D reservoir model.

    For salt caverns:
        porosity = 1.0, total_compressibility = 0.0, is_cavern = True
    """

    area: float  # m² areal extent
    thickness: float  # m net pay or cavern height
    porosity: float  # fraction (1.0 for caverns)
    total_compressibility: float  # Pa⁻¹ total rock+fluid compressibility (0.0 for caverns)
    permeability: float  # m² intrinsic permeability
    reference_pressure: float  # Pa reference pressure
    reference_temperature: float  # K reference geothermal temperature
    well_radius: float  # m wellbore radius
    drainage_radius: float  # m effective drainage radius
    is_cavern: bool = False
    n_wells: int = 1
    pore_model: PoreVolumeModel | None = None
    isotherm: IsothermModel | None = None
    m_sorbent: float = 0.0  # kg sorbent mass for MOF/adsorbent storage


# Alias for naming consistency
ReservoirProperties0D = ReservoirProperties


@dataclass
class ReservoirState:
    """Instantaneous thermodynamic state of the reservoir control volume."""

    pressure: float  # Pa
    temperature: float  # K
    mass: float  # kg


@dataclass
class SchedulePeriod:
    """One operational period in an injection/withdrawal schedule."""

    name: str
    mode: str  # 'inject' | 'withdraw' | 'shut_in'
    duration_days: float  # days
    bhp_target: float  # Pa (P_max for inject, P_min for withdraw)
    max_mass_rate: float  # kg/s maximum surface rate


# Alias
SchedulePeriod0D = SchedulePeriod


@dataclass
class SimulationResult:
    """Simulation time series and aggregated metrics."""

    times: list[float] = field(default_factory=list)  # s
    pressures: list[float] = field(default_factory=list)  # Pa
    temperatures: list[float] = field(default_factory=list)  # K
    masses: list[float] = field(default_factory=list)  # kg
    mass_rates: list[float] = field(default_factory=list)  # kg/s
    metrics: dict[str, float] = field(default_factory=dict)


# Alias
SimulationResult0D = SimulationResult


class UHSReservoir0D:
    """
    Zero-dimensional lumped-parameter UHS reservoir model.
    Enforces strict BHP clamping to eliminate linearisation divergence.
    """

    def __init__(
        self,
        props: ReservoirProperties,
        eos: PengRobinsonH2 | None = None,
        tau_days: float | None = None,
    ) -> None:
        self.props = props
        self.eos = eos or PengRobinsonH2()

        if tau_days is None:
            tau_days = 3.0 if props.is_cavern else 30.0
        self.tau_s: float = tau_days * 86400.0

        self.V_bulk: float = props.area * props.thickness
        self.V_pore0: float = self.V_bulk * props.porosity

        if props.pore_model is None:
            model_type = "rigid" if props.is_cavern else "elastic"
            self.pore_model = PoreVolumeModel(
                model=model_type,
                V0=self.V_pore0,
                P0=props.reference_pressure,
                total_compressibility=props.total_compressibility,
            )
        else:
            self.pore_model = props.pore_model

        self.isotherm = props.isotherm or IsothermModel(model="none")

    def pore_volume(self, P: float) -> float:
        """Pore volume [m³] at pressure P."""
        return self.pore_model.volume(P)

    def mass_from_state(self, P: float, T: float) -> float:
        """Total gas mass [kg] including free gas and adsorbed gas."""
        Vp = self.pore_volume(P)
        rho = self.eos.density(P, T)
        m_free = rho * Vp

        m_ads = 0.0
        if self.props.m_sorbent > 0.0 and self.isotherm.model != "none":
            moles_per_kg = self.isotherm.adsorbed_moles_per_kg(P, T, self.eos)
            m_ads = moles_per_kg * self.props.m_sorbent * 2.01594e-3

        return m_free + m_ads

    def _dm_dP(self, P: float, T: float, dP: float = 1.0e4) -> float:
        P1 = max(P - dP, 1e3)
        P2 = P + dP
        return (self.mass_from_state(P2, T) - self.mass_from_state(P1, T)) / (P2 - P1)

    def _productivity_index(self) -> float:
        if self.props.is_cavern:
            return 1.0e6
        k = self.props.permeability
        h = self.props.thickness
        re = max(self.props.drainage_radius, self.props.well_radius * 2.0)
        rw = self.props.well_radius
        mu = 9.0e-6  # Pa·s
        ln_re_rw = math.log(re / rw)
        PI_single = 2.0 * math.pi * k * h / (mu * max(ln_re_rw, 1e-3))
        return PI_single * self.props.n_wells

    def _mass_rate(self, state: ReservoirState, period: SchedulePeriod) -> float:
        if period.mode == "shut_in":
            return 0.0

        PI = self._productivity_index()
        rho = self.eos.density(state.pressure, state.temperature)
        bhp = period.bhp_target
        qmax = period.max_mass_rate

        if period.mode == "inject":
            dp = max(bhp - state.pressure, 0.0)
            q_vol = PI * dp
            q_kg = min(q_vol * rho, qmax)
            return max(q_kg, 0.0)

        if period.mode == "withdraw":
            dp = max(state.pressure - bhp, 0.0)
            q_vol = PI * dp
            q_kg = min(q_vol * rho, qmax)
            return -max(q_kg, 0.0)

        return 0.0

    def step(
        self, state: ReservoirState, period: SchedulePeriod, dt: float
    ) -> tuple[ReservoirState, float]:
        """Advance reservoir state by time step dt [s]."""
        P = state.pressure
        T = state.temperature

        q_mass = self._mass_rate(state, period)

        if abs(q_mass) < 1e-18:
            new_P = P
        else:
            dm = q_mass * dt
            dm_dP = self._dm_dP(P, T)
            if abs(dm_dP) < 1e-12:
                m_cur = self.mass_from_state(P, T)
                dm_dP = m_cur / max(P, 1e3)
            dP = dm / dm_dP
            new_P = P + dP

        # BHP Clamping
        if period.mode == "inject":
            new_P = min(new_P, period.bhp_target)
        elif period.mode == "withdraw":
            new_P = max(new_P, period.bhp_target)

        new_P = max(new_P, 1.0e3)

        new_T = coupled_temperature_update(
            T_old=T,
            P_old=max(P, 1.0),
            P_new=new_P,
            T_geothermal=self.props.reference_temperature,
            dt_s=dt,
            tau_s=self.tau_s,
            eos=self.eos,
        )

        new_m = self.mass_from_state(new_P, new_T)
        return ReservoirState(pressure=new_P, temperature=new_T, mass=new_m), q_mass

    def run_schedule(
        self,
        initial_state: ReservoirState,
        schedule: list[SchedulePeriod],
        dt_s: float,
        collect_time_series: bool = True,
    ) -> SimulationResult:
        result = SimulationResult()
        state = initial_state
        t = 0.0

        P_max = state.pressure
        P_min = state.pressure
        m_max = state.mass
        m_min = state.mass
        inj_kg = 0.0
        wdr_kg = 0.0

        if collect_time_series:
            result.times.append(t)
            result.pressures.append(state.pressure)
            result.temperatures.append(state.temperature)
            result.masses.append(state.mass)
            result.mass_rates.append(0.0)

        for period in schedule:
            n_steps = max(int(period.duration_days * 86400.0 / dt_s), 1)
            for _ in range(n_steps):
                new_state, q = self.step(state, period, dt_s)
                t += dt_s

                if q >= 0.0:
                    inj_kg += q * dt_s
                else:
                    wdr_kg -= q * dt_s

                P_max = max(P_max, new_state.pressure)
                P_min = min(P_min, new_state.pressure)
                m_max = max(m_max, new_state.mass)
                m_min = min(m_min, new_state.mass)
                state = new_state

                if collect_time_series:
                    result.times.append(t)
                    result.pressures.append(state.pressure)
                    result.temperatures.append(state.temperature)
                    result.masses.append(state.mass)
                    result.mass_rates.append(q)

        wg = m_max - m_min
        eff = wg / m_max * 100.0 if m_max > 0.0 else 0.0

        result.metrics = {
            "max_pressure_Pa": P_max,
            "min_pressure_Pa": P_min,
            "max_mass_kg": m_max,
            "min_mass_kg": m_min,
            "cushion_gas_kg": m_min,
            "working_gas_kg": wg,
            "storage_efficiency_pct": eff,
            "total_injected_kg": inj_kg,
            "total_withdrawn_kg": wdr_kg,
            "cycle_balance_kg": inj_kg - wdr_kg,
            "final_pressure_Pa": state.pressure,
            "final_temperature_K": state.temperature,
            "final_mass_kg": state.mass,
        }
        return result


H2Reservoir0DModel = UHSReservoir0D


# ── Canonical Scenarios ──────────────────────────────────────────────────────
def build_teesside_cavern() -> ReservoirProperties:
    """Teesside salt cavern benchmark (70,000 m³, 40–120 bar)."""
    cavern_volume = 70_000.0
    height = 100.0
    area = cavern_volume / height
    return ReservoirProperties(
        area=area,
        thickness=height,
        porosity=1.0,
        total_compressibility=0.0,
        permeability=1.0e-9,
        reference_pressure=bar_to_pa(40.0),
        reference_temperature=293.15,
        well_radius=0.20,
        drainage_radius=50.0,
        is_cavern=True,
        n_wells=1,
    )


def build_okoroafor_reservoir() -> ReservoirProperties:
    """Okoroafor depleted gas reservoir benchmark (80 bar, 316 K, 100 mD)."""
    re = 500.0
    return ReservoirProperties(
        area=math.pi * re**2,
        thickness=100.0,
        porosity=0.20,
        total_compressibility=5.0e-10,
        permeability=mD_to_m2(100.0),
        reference_pressure=bar_to_pa(80.0),
        reference_temperature=316.0,
        well_radius=0.10,
        drainage_radius=re,
        is_cavern=False,
        n_wells=1,
    )


def build_rough_reservoir() -> ReservoirProperties:
    """Rough-style depleted gas field benchmark (150 bar, 369 K, 100 mD)."""
    re = 500.0
    return ReservoirProperties(
        area=math.pi * re**2,
        thickness=30.0,
        porosity=0.15,
        total_compressibility=4.0e-10,
        permeability=mD_to_m2(100.0),
        reference_pressure=bar_to_pa(150.0),
        reference_temperature=369.0,
        well_radius=0.15,
        drainage_radius=re,
        is_cavern=False,
        n_wells=1,
    )


def schedule_cavern(P_ref: float) -> list[SchedulePeriod]:
    bar = 1.0e5
    P_max = min(2.0 * P_ref, P_ref + 80.0 * bar)
    P_min = max(0.35 * P_ref, P_ref - 25.0 * bar)
    qmax = 200.0
    sched = [SchedulePeriod("cushion_fill", "inject", 90.0, P_max, qmax)]
    for i in range(180):
        sched.extend(
            [
                SchedulePeriod(f"c{i + 1}_inj", "inject", 0.5, P_max, qmax),
                SchedulePeriod(f"c{i + 1}_wdr", "withdraw", 0.4, P_min, qmax),
                SchedulePeriod(f"c{i + 1}_sht", "shut_in", 0.1, 0.0, 0.0),
            ]
        )
    return sched


def schedule_porous_reservoir(P_ref: float) -> list[SchedulePeriod]:
    bar = 1.0e5
    P_max = min(1.35 * P_ref, P_ref + 50.0 * bar)
    P_min = max(0.60 * P_ref, P_ref - 40.0 * bar)
    qmax = 35.0
    sched = [SchedulePeriod("cushion_fill", "inject", 730.0, P_max, qmax)]
    for i in range(3):
        sched.extend(
            [
                SchedulePeriod(f"c{i + 1}_inj", "inject", 210.0, P_max, qmax),
                SchedulePeriod(f"c{i + 1}_wdr", "withdraw", 120.0, P_min, qmax),
                SchedulePeriod(f"c{i + 1}_sht", "shut_in", 35.0, 0.0, 0.0),
            ]
        )
    return sched


# ── Monte Carlo Geological Uncertainty ───────────────────────────────────────
@dataclass
class MonteCarloConfig:
    n_realizations: int
    porosity_mean: float
    porosity_std: float
    thickness_mean: float
    thickness_std: float
    log10_perm_mean: float
    log10_perm_std: float
    geothermal_T_mean: float
    geothermal_T_std: float
    area_mean: float
    area_std: float


def run_monte_carlo(
    base_props: ReservoirProperties,
    schedule_fn: Callable[[float], list[SchedulePeriod]],
    mc_cfg: MonteCarloConfig,
    dt_days: float = 1.0,
    tau_days: float | None = None,
    seed: int = 42,
) -> dict[str, dict[str, float]]:
    rng = random.Random(seed)
    eos = PengRobinsonH2()
    dt = dt_days * 86400.0
    metric_values: dict[str, list[float]] = {}

    for _ in range(mc_cfg.n_realizations):
        phi = max(0.01, min(0.40, rng.gauss(mc_cfg.porosity_mean, mc_cfg.porosity_std)))
        h = max(5.0, rng.gauss(mc_cfg.thickness_mean, mc_cfg.thickness_std))
        logk = rng.gauss(mc_cfg.log10_perm_mean, mc_cfg.log10_perm_std)
        k = 10.0**logk
        T_g = max(250.0, rng.gauss(mc_cfg.geothermal_T_mean, mc_cfg.geothermal_T_std))
        A = max(1.0e4, rng.gauss(mc_cfg.area_mean, mc_cfg.area_std))

        rp = ReservoirProperties(
            area=A,
            thickness=h,
            porosity=phi,
            total_compressibility=base_props.total_compressibility,
            permeability=k,
            reference_pressure=base_props.reference_pressure,
            reference_temperature=T_g,
            well_radius=base_props.well_radius,
            drainage_radius=base_props.drainage_radius,
            is_cavern=base_props.is_cavern,
            n_wells=base_props.n_wells,
        )

        model = UHSReservoir0D(rp, eos=eos, tau_days=tau_days)
        P0 = rp.reference_pressure
        T0 = rp.reference_temperature
        m0 = model.mass_from_state(P0, T0)
        state = ReservoirState(pressure=P0, temperature=T0, mass=m0)
        sched = schedule_fn(P0)
        result = model.run_schedule(state, sched, dt, collect_time_series=False)

        for k_met, v_met in result.metrics.items():
            metric_values.setdefault(k_met, []).append(v_met)

    stats: dict[str, dict[str, float]] = {}
    for name, vals in metric_values.items():
        sv = sorted(vals)
        mean_v = sum(sv) / len(sv)
        var_v = sum((x - mean_v) ** 2 for x in sv) / len(sv)

        stats[name] = {
            "mean": mean_v,
            "std": math.sqrt(var_v),
            "min": sv[0],
            "p10": percentile(sv, 10.0),
            "p50": percentile(sv, 50.0),
            "p90": percentile(sv, 90.0),
            "max": sv[-1],
        }
    return stats


def export_timeseries_csv(result: SimulationResult, filepath: str) -> None:
    """Export time series to CSV with dual units."""
    if not result.times:
        return
    os.makedirs(os.path.dirname(os.path.abspath(filepath)), exist_ok=True)
    with open(filepath, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(
            [
                "time_s",
                "time_days",
                "pressure_Pa",
                "pressure_bar",
                "pressure_MPa",
                "pressure_psi",
                "T_K",
                "T_C",
                "T_F",
                "mass_kg",
                "mass_tonnes",
                "mass_short_tons",
                "mass_lb",
                "mass_rate_kg_s",
                "mass_rate_tonnes_day",
            ]
        )
        for t, P, T, m, q in zip(
            result.times,
            result.pressures,
            result.temperatures,
            result.masses,
            result.mass_rates,
            strict=False,
        ):
            w.writerow(
                [
                    f"{t:.1f}",
                    f"{t / 86400:.5f}",
                    f"{P:.3e}",
                    f"{pa_to_bar(P):.4f}",
                    f"{pa_to_MPa(P):.4f}",
                    f"{pa_to_psi(P):.4f}",
                    f"{T:.3f}",
                    f"{K_to_C(T):.3f}",
                    f"{K_to_F(T):.3f}",
                    f"{m:.6e}",
                    f"{kg_to_tonnes(m):.6e}",
                    f"{kg_to_short_tons(m):.6e}",
                    f"{kg_to_lb(m):.6e}",
                    f"{q:.6e}",
                    f"{kg_to_tonnes(q * 86400):.6e}",
                ]
            )
