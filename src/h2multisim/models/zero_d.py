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
    volume_translation: float = 0.0  # m³/mol Peneloux volume translation (e.g. -4.0e-6)
    peneloux: bool = False  # If True and volume_translation == 0, sets volume_translation = -4.0e-6

    def __post_init__(self) -> None:
        if self.peneloux and self.volume_translation == 0.0:
            self.volume_translation = -4.0e-6


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
        if eos is None:
            vt = getattr(props, "volume_translation", 0.0)
            self.eos = PengRobinsonH2(volume_translation=vt)
        else:
            self.eos = eos

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

    def _mass_rate_at(self, P: float, T: float, period: SchedulePeriod) -> float:
        """Signed well mass rate [kg/s] at pressure P and temperature T."""
        if period.mode == "shut_in":
            return 0.0

        PI = self._productivity_index()
        rho = self.eos.density(P, T)
        bhp = period.bhp_target
        qmax = period.max_mass_rate

        if period.mode == "inject":
            dp = max(bhp - P, 0.0)
            q_kg = min(PI * dp * rho, qmax)
            return max(q_kg, 0.0)

        if period.mode == "withdraw":
            dp = max(P - bhp, 0.0)
            q_kg = min(PI * dp * rho, qmax)
            return -max(q_kg, 0.0)

        return 0.0

    def _mass_rate(self, state: ReservoirState, period: SchedulePeriod) -> float:
        """Signed well mass rate [kg/s] at a reservoir state (diagnostic helper)."""
        return self._mass_rate_at(state.pressure, state.temperature, period)

    def _step_residual(
        self,
        P_new: float,
        P_old: float,
        T_old: float,
        m_old: float,
        period: SchedulePeriod,
        dt: float,
    ) -> tuple[float, float, float]:
        """Mass-balance residual and associated temperature/rate at P_new."""
        T_new = coupled_temperature_update(
            T_old=T_old,
            P_old=max(P_old, 1.0),
            P_new=max(P_new, 1.0),
            T_geothermal=self.props.reference_temperature,
            dt_s=dt,
            tau_s=self.tau_s,
            eos=self.eos,
        )
        q = self._mass_rate_at(P_new, T_new, period)
        m_new = self.mass_from_state(P_new, T_new)
        return m_new - m_old - q * dt, q, T_new

    def step(
        self, state: ReservoirState, period: SchedulePeriod, dt: float
    ) -> tuple[ReservoirState, float]:
        """
        Advance the reservoir state by dt [s] with an exact mass balance.

        Solving ``m(P_new, T_new) - m_old = q(P_new, T_new)·dt`` removes the
        explicit linearisation and the old post-hoc BHP clamp: the well rate is
        the rate that actually changes the inventory, so the reported throughput
        and the stored mass cannot diverge.
        """
        P_old = state.pressure
        T_old = state.temperature
        m_old = self.mass_from_state(P_old, T_old)

        if period.mode == "shut_in":
            lo = max(1.0e3, P_old * 0.5)
            hi = max(P_old * 1.5, lo * 1.01)
        else:
            lo = max(1.0e3, min(P_old, period.bhp_target) * 0.98)
            hi = max(P_old, period.bhp_target) * 1.02 + 1.0e3

        f_lo = self._step_residual(lo, P_old, T_old, m_old, period, dt)[0]
        f_hi = self._step_residual(hi, P_old, T_old, m_old, period, dt)[0]
        for _ in range(80):
            if f_lo <= 0.0:
                break
            lo = max(1.0e3, lo * 0.5)
            f_lo = self._step_residual(lo, P_old, T_old, m_old, period, dt)[0]
        for _ in range(80):
            if f_hi >= 0.0:
                break
            hi *= 2.0
            f_hi = self._step_residual(hi, P_old, T_old, m_old, period, dt)[0]

        P_new = 0.5 * (lo + hi)
        q = 0.0
        T_new = T_old
        for _ in range(80):
            f, q, T_new = self._step_residual(P_new, P_old, T_old, m_old, period, dt)
            if abs(f) <= 1.0e-9 * max(abs(m_old), 1.0):
                break
            if f > 0.0:
                hi = P_new
            else:
                lo = P_new
            if hi - lo <= 1.0e-7 * max(P_new, 1.0):
                break

            dP = max(1.0, 1.0e-6 * P_new)
            f2 = self._step_residual(P_new + dP, P_old, T_old, m_old, period, dt)[0]
            slope = (f2 - f) / dP
            if slope > 0.0 and math.isfinite(slope):
                P_try = P_new - f / slope
            else:
                P_try = 0.5 * (lo + hi)
            if not (lo < P_try < hi):
                P_try = 0.5 * (lo + hi)
            P_new = max(P_try, 1.0e3)

        m_new = self.mass_from_state(P_new, T_new)
        q_report = (m_new - m_old) / dt if dt > 0.0 else 0.0
        return ReservoirState(pressure=P_new, temperature=T_new, mass=m_new), q_report

    def run_schedule(
        self,
        initial_state: ReservoirState,
        schedule: list[SchedulePeriod],
        dt_s: float,
        collect_time_series: bool = True,
    ) -> SimulationResult:
        result = SimulationResult()
        m0 = self.mass_from_state(initial_state.pressure, initial_state.temperature)
        state = ReservoirState(
            pressure=initial_state.pressure,
            temperature=initial_state.temperature,
            mass=m0,
        )
        t = 0.0

        P_max = state.pressure
        P_min = state.pressure
        m_max = m0
        m_min = m0
        inj_kg = 0.0
        wdr_kg = 0.0

        if collect_time_series:
            result.times.append(t)
            result.pressures.append(state.pressure)
            result.temperatures.append(state.temperature)
            result.masses.append(state.mass)
            result.mass_rates.append(0.0)

        for period in schedule:
            n_steps = max(round(period.duration_days * 86400.0 / dt_s), 1)
            dt_act = period.duration_days * 86400.0 / n_steps
            for _ in range(n_steps):
                new_state, q = self.step(state, period, dt_act)
                t += dt_act

                if q >= 0.0:
                    inj_kg += q * dt_act
                else:
                    wdr_kg -= q * dt_act

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
        net_kg = inj_kg - wdr_kg
        mb_err = abs((state.mass - m0) - net_kg) / max(abs(net_kg), 1.0)

        result.metrics = {
            "max_pressure_Pa": P_max,
            "min_pressure_Pa": P_min,
            "max_mass_kg": m_max,
            "min_mass_kg": m_min,
            "initial_mass_kg": m0,
            "cushion_gas_kg": m_min,
            "working_gas_kg": wg,
            "storage_efficiency_pct": eff,
            "total_injected_kg": inj_kg,
            "total_withdrawn_kg": wdr_kg,
            "net_stored_kg": net_kg,
            "cycled_gas_kg": inj_kg + wdr_kg,
            "cycle_balance_kg": net_kg,
            "mass_balance_error": mb_err,
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
    """
    Canonical salt-cavern schedule.

    ``P_ref`` is the minimum operating pressure.  The cavern is filled by
    80 bar above that (Teesside: 40 → 120 bar), then cycled 180 times.  The
    15 kg/s rate limit is representative of a single Teesside well; the exact
    step solver lets the pressure targets, not the rate cap, terminate each
    half-cycle.
    """
    bar = 1.0e5
    P_min = P_ref
    P_max = P_ref + 80.0 * bar
    qmax = 15.0
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


def schedule_rough(P_ref: float) -> list[SchedulePeriod]:
    """
    Rough-style large depleted gas field schedule.

    Operating window:
        P_min = max(0.65 · P_ref, 100 bar)
        P_max = min(1.45 · P_ref, 220 bar)

    Phase 1: 270-day cushion gas fill (q_max = 90 kg/s)
    Phase 2: 3 annual cycles (210d injection, 120d withdrawal, 35d shut-in)
    """
    bar = 1.0e5
    P_max = min(1.45 * P_ref, 220.0 * bar)
    P_min = max(0.65 * P_ref, 100.0 * bar)
    qmax = 90.0
    sched = [SchedulePeriod("cushion_fill", "inject", 270.0, P_max, qmax)]
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
