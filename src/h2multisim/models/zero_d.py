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

from h2multisim.models._common import check_inputs, validate_mode
from h2multisim.thermo import (
    IsothermModel,
    PengRobinsonH2,
    PengRobinsonMixture,
    PoreVolumeModel,
    RelativePermeabilityModel,
    coupled_temperature_update,
    normalize_composition,
    viscosity_h2,
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
    cushion_gas: str = "H2"  # Native cushion gas species ('H2', 'CH4', 'CO2', 'N2')
    initial_composition: dict[str, float] | None = None  # Optional initial mole fractions
    injected_composition: dict[str, float] | None = None  # Optional injected stream mole fractions
    plume_segregation: float = (
        0.0  # [0, 1] buoyant plume segregation during withdrawal (0 = well-mixed)
    )
    swr: float = 0.0  # [-] connate / residual immobile water saturation S_wr
    sgc: float = 0.0  # [-] critical / trapped gas saturation S_gc
    rel_perm: RelativePermeabilityModel | None = (
        None  # Optional two-phase relative permeability model
    )

    def __post_init__(self) -> None:
        check_inputs(
            "ReservoirProperties",
            positive={
                "area": self.area,
                "thickness": self.thickness,
                "permeability": self.permeability,
                "reference_pressure": self.reference_pressure,
                "reference_temperature": self.reference_temperature,
                "well_radius": self.well_radius,
                "drainage_radius": self.drainage_radius,
                "n_wells": float(self.n_wells),
            },
            non_negative={
                "total_compressibility": self.total_compressibility,
                "m_sorbent": self.m_sorbent,
            },
            fraction={"porosity": self.porosity},
        )
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
    composition: dict[str, float] | None = None  # Active mole-fraction composition


@dataclass
class SchedulePeriod:
    """One operational period in an injection/withdrawal schedule."""

    name: str
    mode: str  # 'inject' | 'withdraw' | 'shut_in'
    duration_days: float  # days
    bhp_target: float  # Pa (P_max for inject, P_min for withdraw)
    max_mass_rate: float  # kg/s maximum surface rate

    def __post_init__(self) -> None:
        validate_mode(self.mode, f"schedule period {self.name!r}")
        check_inputs(
            f"schedule period {self.name!r}",
            positive={"duration_days": self.duration_days},
            non_negative={"bhp_target": self.bhp_target, "max_mass_rate": self.max_mass_rate},
        )


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
    h2_mole_fractions: list[float] = field(default_factory=list)  # [-] H2 mole fraction y_H2(t)
    wobbe_indices: list[float] = field(default_factory=list)  # MJ/sm³ Wobbe Index
    hhv_values: list[float] = field(default_factory=list)  # MJ/sm³ volumetric HHV
    metrics: dict[str, float] = field(default_factory=dict)


# Alias
SimulationResult0D = SimulationResult


class UHSReservoir0D:
    """
    Zero-dimensional lumped-parameter UHS reservoir model.
    Solves the coupled nonlinear mass-balance root at each time step,
    including multicomponent cushion-gas dilution (H2 / CH4 / CO2 / N2).
    """

    def __init__(
        self,
        props: ReservoirProperties,
        eos: PengRobinsonH2 | PengRobinsonMixture | None = None,
        tau_days: float | None = None,
    ) -> None:
        self.props = props
        cg = getattr(props, "cushion_gas", "H2") or "H2"
        init_comp = normalize_composition(getattr(props, "initial_composition", None) or {cg: 1.0})
        inj_comp = normalize_composition(
            getattr(props, "injected_composition", None) or {"H2": 1.0}
        )
        self._init_comp: dict[str, float] = init_comp
        self._inj_comp: dict[str, float] = inj_comp
        self._is_multicomponent: bool = not (
            len(init_comp) == 1 and "H2" in init_comp and len(inj_comp) == 1 and "H2" in inj_comp
        )

        if eos is None:
            vt = getattr(props, "volume_translation", 0.0)
            if self._is_multicomponent:
                self.eos: PengRobinsonH2 | PengRobinsonMixture = PengRobinsonMixture(
                    composition=init_comp,
                    peneloux=getattr(props, "peneloux", False),
                    volume_translation=vt
                    if vt != 0.0 and not getattr(props, "peneloux", False)
                    else None,
                )
            else:
                self.eos = PengRobinsonH2(volume_translation=vt)
        else:
            self.eos = eos
            if isinstance(eos, PengRobinsonMixture) and not eos.is_pure_h2():
                self._is_multicomponent = True
                self._init_comp = dict(eos.composition)

        self._mix_helper = (
            self.eos
            if isinstance(self.eos, PengRobinsonMixture)
            else PengRobinsonMixture(
                composition=self._init_comp,
                peneloux=getattr(props, "peneloux", False),
            )
        )

        if tau_days is None:
            tau_days = 3.0 if props.is_cavern else 30.0
        self.tau_s: float = tau_days * 86400.0

        self.rel_perm: RelativePermeabilityModel = getattr(
            props, "rel_perm", None
        ) or RelativePermeabilityModel(
            swr=float(getattr(props, "swr", 0.0)),
            sgc=float(getattr(props, "sgc", 0.0)),
        )
        self.V_bulk: float = props.area * props.thickness
        self.V_pore0: float = (
            self.V_bulk * props.porosity * self.rel_perm.effective_porosity_factor()
        )

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
        self._last_step_converged: bool = True

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

    def _productivity_index(
        self,
        P: float | None = None,
        T: float | None = None,
        rho: float | None = None,
    ) -> float:
        if self.props.is_cavern:
            return 1.0e6
        k = self.props.permeability * self.rel_perm.effective_permeability_factor()
        h = self.props.thickness
        re = max(self.props.drainage_radius, self.props.well_radius * 2.0)
        rw = self.props.well_radius
        P_val = P if P is not None else self.props.reference_pressure
        T_val = T if T is not None else self.props.reference_temperature
        mu = viscosity_h2(T_val, P_val, rho=rho, eos=self.eos)
        ln_re_rw = math.log(re / rw)
        PI_single = 2.0 * math.pi * k * h / (mu * max(ln_re_rw, 1e-3))
        return PI_single * self.props.n_wells

    def _mass_rate_at(self, P: float, T: float, period: SchedulePeriod) -> float:
        """Signed well mass rate [kg/s] at pressure P and temperature T."""
        if period.mode == "shut_in":
            return 0.0

        rho = self.eos.density(P, T)
        PI = self._productivity_index(P, T, rho)
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

    def _trial_composition(
        self,
        comp_m_old: dict[str, float],
        w_inj: dict[str, float],
        w_prod: dict[str, float],
        q_val: float,
        dt: float,
    ) -> dict[str, float]:
        dm = q_val * dt
        all_sp = set(comp_m_old) | set(w_inj) | set(w_prod)
        if dm >= 0.0:
            m_map = {sp: comp_m_old.get(sp, 0.0) + w_inj.get(sp, 0.0) * dm for sp in all_sp}
        else:
            m_map = {
                sp: max(comp_m_old.get(sp, 0.0) + w_prod.get(sp, 0.0) * dm, 0.0) for sp in all_sp
            }
        return PengRobinsonMixture.mole_fractions_from_mass(m_map)

    def _step_residual(
        self,
        P_new: float,
        P_old: float,
        T_old: float,
        m_old: float,
        period: SchedulePeriod,
        dt: float,
        comp_info: tuple[dict[str, float], dict[str, float], dict[str, float], dict[str, float]]
        | None = None,
    ) -> tuple[float, float, float]:
        """Mass-balance residual and associated temperature/rate at P_new."""
        if comp_info is not None and isinstance(self.eos, PengRobinsonMixture):
            y_old, comp_m_old, w_inj, w_prod = comp_info
            self.eos.set_composition(y_old)
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
            y_trial = self._trial_composition(comp_m_old, w_inj, w_prod, q, dt)
            self.eos.set_composition(y_trial)
            q = self._mass_rate_at(P_new, T_new, period)
            y_trial = self._trial_composition(comp_m_old, w_inj, w_prod, q, dt)
            self.eos.set_composition(y_trial)
            m_new = self.mass_from_state(P_new, T_new)
            return m_new - m_old - q * dt, q, T_new

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
        explicit linearisation and the old post-hoc BHP clamp.  The returned
        rate is the step-average rate (m_new − m_old)/dt: when the well hits its
        BHP target mid-step the instantaneous well rate is effectively
        discontinuous, and the average rate is the physically meaningful one.
        Because that makes the run-level mass balance an identity, solver
        health is tracked separately: ``_last_step_converged`` is False when
        the root was neither bracketed nor located within tolerance.
        """
        P_old = state.pressure
        T_old = state.temperature

        comp_info = None
        y_old = normalize_composition(state.composition or self._init_comp)
        if self._is_multicomponent and isinstance(self.eos, PengRobinsonMixture):
            self.eos.set_composition(y_old)
            m_old = self.mass_from_state(P_old, T_old)
            w_old = self.eos.mass_fractions(y_old)
            comp_m_old = {sp: w * m_old for sp, w in w_old.items()}
            w_inj = self.eos.mass_fractions(self._inj_comp)
            seg = min(max(float(getattr(self.props, "plume_segregation", 0.0)), 0.0), 1.0)
            if seg > 0.0 and comp_m_old.get("H2", 0.0) > 1e-6:
                all_sp = set(w_old) | set(w_inj)
                w_prod = {
                    sp: (1.0 - seg) * w_old.get(sp, 0.0) + seg * w_inj.get(sp, 0.0) for sp in all_sp
                }
            else:
                w_prod = w_old
            comp_info = (y_old, comp_m_old, w_inj, w_prod)
        else:
            m_old = self.mass_from_state(P_old, T_old)

        if period.mode == "shut_in":
            lo = max(1.0e3, P_old * 0.5)
            hi = max(P_old * 1.5, lo * 1.01)
        else:
            lo = max(1.0e3, min(P_old, period.bhp_target) * 0.98)
            hi = max(P_old, period.bhp_target) * 1.02 + 1.0e3

        f_lo = self._step_residual(lo, P_old, T_old, m_old, period, dt, comp_info)[0]
        f_hi = self._step_residual(hi, P_old, T_old, m_old, period, dt, comp_info)[0]
        for _ in range(80):
            if f_lo <= 0.0:
                break
            lo = max(1.0e3, lo * 0.5)
            f_lo = self._step_residual(lo, P_old, T_old, m_old, period, dt, comp_info)[0]
        for _ in range(80):
            if f_hi >= 0.0:
                break
            hi *= 2.0
            f_hi = self._step_residual(hi, P_old, T_old, m_old, period, dt, comp_info)[0]
        converged = f_lo <= 0.0 <= f_hi
        located = False

        P_new = 0.5 * (lo + hi)
        q = 0.0
        T_new = T_old
        for _ in range(80):
            f, q, T_new = self._step_residual(P_new, P_old, T_old, m_old, period, dt, comp_info)
            if abs(f) <= 1.0e-9 * max(abs(m_old), 1.0):
                located = True
                break
            if f > 0.0:
                hi = P_new
            else:
                lo = P_new
            if hi - lo <= 1.0e-7 * max(P_new, 1.0):
                located = True
                break

            dP = max(1.0, 1.0e-6 * P_new)
            f2 = self._step_residual(P_new + dP, P_old, T_old, m_old, period, dt, comp_info)[0]
            slope = (f2 - f) / dP
            if slope > 0.0 and math.isfinite(slope):
                P_try = P_new - f / slope
            else:
                P_try = 0.5 * (lo + hi)
            if not (lo < P_try < hi):
                P_try = 0.5 * (lo + hi)
            P_new = max(P_try, 1.0e3)

        _, q, T_new = self._step_residual(P_new, P_old, T_old, m_old, period, dt, comp_info)
        y_new: dict[str, float] | None = state.composition
        if comp_info is not None and isinstance(self.eos, PengRobinsonMixture):
            y_new = dict(self.eos.composition)
        m_new = self.mass_from_state(P_new, T_new)
        q_report = (m_new - m_old) / dt if dt > 0.0 else 0.0
        self._last_step_converged = converged and located

        return (
            ReservoirState(pressure=P_new, temperature=T_new, mass=m_new, composition=y_new),
            q_report,
        )

    def run_schedule(
        self,
        initial_state: ReservoirState,
        schedule: list[SchedulePeriod],
        dt_s: float,
        collect_time_series: bool = True,
    ) -> SimulationResult:
        result = SimulationResult()
        y0 = normalize_composition(initial_state.composition or self._init_comp)
        if isinstance(self.eos, PengRobinsonMixture):
            self.eos.set_composition(y0)
        m0 = self.mass_from_state(initial_state.pressure, initial_state.temperature)
        state = ReservoirState(
            pressure=initial_state.pressure,
            temperature=initial_state.temperature,
            mass=m0,
            composition=y0 if self._is_multicomponent else initial_state.composition,
        )
        t = 0.0

        P_max = state.pressure
        P_min = state.pressure
        m_max = m0
        m_min = m0
        inj_kg = 0.0
        wdr_kg = 0.0
        unconverged_steps = 0

        y_cur = y0
        y_h2_init = y_cur.get("H2", 1.0)
        min_wdr_h2 = y_h2_init

        if collect_time_series:
            result.times.append(t)
            result.pressures.append(state.pressure)
            result.temperatures.append(state.temperature)
            result.masses.append(state.mass)
            result.mass_rates.append(0.0)
            result.h2_mole_fractions.append(y_h2_init)
            result.wobbe_indices.append(self._mix_helper.wobbe_index(y_cur))
            result.hhv_values.append(self._mix_helper.heating_values_volumetric(y_cur)[0])

        for period in schedule:
            n_steps = max(round(period.duration_days * 86400.0 / dt_s), 1)
            dt_act = period.duration_days * 86400.0 / n_steps
            for _ in range(n_steps):
                new_state, q = self.step(state, period, dt_act)
                if not self._last_step_converged:
                    unconverged_steps += 1
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

                y_cur = normalize_composition(state.composition or self._init_comp)
                y_h2_cur = y_cur.get("H2", 1.0)
                if q < 0.0:
                    min_wdr_h2 = min(min_wdr_h2, y_h2_cur)

                if collect_time_series:
                    result.times.append(t)
                    result.pressures.append(state.pressure)
                    result.temperatures.append(state.temperature)
                    result.masses.append(state.mass)
                    result.mass_rates.append(q)
                    result.h2_mole_fractions.append(y_h2_cur)
                    result.wobbe_indices.append(self._mix_helper.wobbe_index(y_cur))
                    result.hhv_values.append(self._mix_helper.heating_values_volumetric(y_cur)[0])

        wg = m_max - m_min
        eff = wg / m_max * 100.0 if m_max > 0.0 else 0.0
        net_kg = inj_kg - wdr_kg
        # Bookkeeping check only: step rates are inventory-derived, so this is ~0 by
        # construction.  Solver health is reported by ``converged``.
        mb_err = abs((state.mass - m0) - net_kg) / max(abs(net_kg), 1.0)
        y_final = normalize_composition(state.composition or self._init_comp)
        if wdr_kg <= 0.0:
            min_wdr_h2 = y_final.get("H2", 1.0)

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
            "unconverged_steps": float(unconverged_steps),
            "converged": 1.0 if unconverged_steps == 0 else 0.0,
            "final_pressure_Pa": state.pressure,
            "final_temperature_K": state.temperature,
            "final_mass_kg": state.mass,
            "initial_h2_purity_mol_pct": y_h2_init * 100.0,
            "final_h2_purity_mol_pct": y_final.get("H2", 1.0) * 100.0,
            "min_withdrawn_h2_purity_mol_pct": min_wdr_h2 * 100.0,
            "final_wobbe_index_MJ_sm3": self._mix_helper.wobbe_index(y_final),
            "final_hhv_MJ_sm3": self._mix_helper.heating_values_volumetric(y_final)[0],
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
    schedule_fn: Callable[[float], list[SchedulePeriod]] | list[SchedulePeriod],
    mc_cfg: MonteCarloConfig,
    dt_days: float = 1.0,
    tau_days: float | None = None,
    seed: int = 42,
) -> dict[str, dict[str, float]]:
    rng = random.Random(seed)
    vt = getattr(base_props, "volume_translation", 0.0)
    eos = PengRobinsonH2(volume_translation=vt)
    dt = dt_days * 86400.0
    metric_values: dict[str, list[float]] = {}

    for _ in range(mc_cfg.n_realizations):
        phi = max(0.01, min(0.40, rng.gauss(mc_cfg.porosity_mean, mc_cfg.porosity_std)))
        h = max(5.0, rng.gauss(mc_cfg.thickness_mean, mc_cfg.thickness_std))
        logk = rng.gauss(mc_cfg.log10_perm_mean, mc_cfg.log10_perm_std)
        k = 10.0**logk
        T_g = max(250.0, rng.gauss(mc_cfg.geothermal_T_mean, mc_cfg.geothermal_T_std))
        A = max(1.0e4, rng.gauss(mc_cfg.area_mean, mc_cfg.area_std))
        re = base_props.drainage_radius if base_props.is_cavern else math.sqrt(A / math.pi)
        sat_factor = (
            base_props.rel_perm or RelativePermeabilityModel(swr=base_props.swr, sgc=base_props.sgc)
        ).effective_porosity_factor()
        pm = (
            PoreVolumeModel(
                model=base_props.pore_model.model,
                V0=A * h * phi * sat_factor,
                P0=base_props.pore_model.P0,
                total_compressibility=base_props.pore_model.total_compressibility,
                P_gate=base_props.pore_model.P_gate,
                k_gate=base_props.pore_model.k_gate,
                V_ratio=base_props.pore_model.V_ratio,
            )
            if base_props.pore_model is not None
            else None
        )

        rp = ReservoirProperties(
            area=A,
            thickness=h,
            porosity=phi,
            total_compressibility=base_props.total_compressibility,
            permeability=k,
            reference_pressure=base_props.reference_pressure,
            reference_temperature=T_g,
            well_radius=base_props.well_radius,
            drainage_radius=re,
            is_cavern=base_props.is_cavern,
            n_wells=base_props.n_wells,
            pore_model=pm,
            isotherm=base_props.isotherm,
            m_sorbent=base_props.m_sorbent,
            volume_translation=base_props.volume_translation,
            peneloux=base_props.peneloux,
            cushion_gas=base_props.cushion_gas,
            initial_composition=base_props.initial_composition,
            injected_composition=base_props.injected_composition,
            plume_segregation=base_props.plume_segregation,
            swr=base_props.swr,
            sgc=base_props.sgc,
            rel_perm=base_props.rel_perm,
        )

        model = UHSReservoir0D(
            rp,
            eos=None if rp.cushion_gas.upper() != "H2" or rp.initial_composition else eos,
            tau_days=tau_days,
        )
        P0 = rp.reference_pressure
        T0 = rp.reference_temperature
        m0 = model.mass_from_state(P0, T0)
        state = ReservoirState(pressure=P0, temperature=T0, mass=m0)
        sched = schedule_fn(P0) if callable(schedule_fn) else schedule_fn
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
