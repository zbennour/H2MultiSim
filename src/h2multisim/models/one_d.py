"""
h2multisim.models.one_d
=======================
One-dimensional (1D) radial finite-volume reservoir model for near-wellbore
dynamics and transient pressure analysis in underground hydrogen storage (UHS).
"""

from __future__ import annotations

import csv
import math
import os
from dataclasses import dataclass, field

from h2multisim.thermo import (
    GAS_COMPONENTS,
    M_H2,
    PengRobinsonH2,
    PengRobinsonMixture,
    RelativePermeabilityModel,
    coupled_temperature_update,
    viscosity_h2,
)
from h2multisim.units import (
    K_to_C,
    bar_to_pa,
    kg_to_tonnes,
    mD_to_m2,
    pa_to_bar,
    pa_to_psi,
)


@dataclass
class RadialReservoirProps:
    """Reservoir properties for the 1D radial model."""

    r_well: float  # m wellbore radius
    r_ext: float  # m outer boundary radius
    h: float  # m net pay thickness
    phi: float  # fraction porosity
    k: float  # m² permeability
    Ct: float  # Pa⁻¹ total compressibility
    P_init: float  # Pa initial reservoir pressure
    T_init: float  # K initial / geothermal temperature
    skin: float = 0.0  # mechanical skin factor
    outer_bc: str = "no_flow"  # 'no_flow' | 'const_pressure'
    tau_days: float = 30.0  # days thermal relaxation constant
    volume_translation: float = 0.0  # m³/mol Peneloux volume translation (e.g. -4.0e-6)
    peneloux: bool = False  # If True and volume_translation == 0, sets volume_translation = -4.0e-6
    cushion_gas: str = "H2"  # Native cushion gas species ('H2', 'CH4', 'CO2', 'N2')
    initial_h2_mole_fraction: float | None = (
        None  # Initial reservoir H2 mole fraction (None -> 1.0 if H2 else 0.0)
    )
    dispersivity_m: float = 0.5  # m radial longitudinal mechanical dispersivity α_L
    molecular_diffusivity_m2_s: float = 1.0e-5  # m²/s effective pore molecular diffusivity D_mol
    swr: float = 0.0  # [-] connate / residual immobile water saturation S_wr
    sgc: float = 0.0  # [-] critical / trapped gas saturation S_gc
    rel_perm: RelativePermeabilityModel | None = (
        None  # Optional two-phase relative permeability model
    )

    def __post_init__(self) -> None:
        if self.peneloux and self.volume_translation == 0.0:
            self.volume_translation = -4.0e-6


ReservoirProperties1D = RadialReservoirProps


@dataclass
class WellControl:
    """Operational period control for 1D radial model."""

    name: str
    mode: str  # 'inject' | 'withdraw' | 'shut_in'
    duration_d: float  # days
    bhp_target: float  # Pa
    q_max: float  # kg/s


SchedulePeriod1D = WellControl


@dataclass
class RadialSnapshot:
    """Radial pressure, temperature, and composition profile at one instant."""

    time_s: float
    r_m: list[float]
    P_Pa: list[float]
    T_K: list[float]
    h2_mole_fraction: list[float] = field(default_factory=list)


@dataclass
class RadialResult:
    """Output from 1D radial simulation."""

    times_s: list[float] = field(default_factory=list)
    BHP_Pa: list[float] = field(default_factory=list)
    q_kg_s: list[float] = field(default_factory=list)
    cum_inj_kg: list[float] = field(default_factory=list)
    cum_wdr_kg: list[float] = field(default_factory=list)
    P_avg_Pa: list[float] = field(default_factory=list)
    T_avg_K: list[float] = field(default_factory=list)
    mass_kg: list[float] = field(default_factory=list)
    well_h2_mole_fraction: list[float] = field(default_factory=list)
    well_wobbe_index: list[float] = field(default_factory=list)
    snapshots: list[RadialSnapshot] = field(default_factory=list)
    metrics: dict[str, float] = field(default_factory=dict)


SimulationResult1D = RadialResult


def build_log_grid(r_w: float, r_e: float, N: int) -> tuple[list[float], list[float]]:
    """Logarithmic radial grid: N cell-centres and N+1 interfaces between r_w and r_e."""
    if N < 3:
        raise ValueError("N must be >= 3")
    ratio = r_e / r_w
    interfaces = [r_w * ratio ** (i / N) for i in range(N + 1)]
    centres = [math.sqrt(interfaces[i] * interfaces[i + 1]) for i in range(N)]
    return centres, interfaces


def thomas_solve(
    a: list[float],
    b: list[float],
    c: list[float],
    d: list[float],
) -> list[float]:
    """
    Solve tridiagonal system a[i]*x[i-1] + b[i]*x[i] + c[i]*x[i+1] = d[i]
    using the Thomas algorithm (O(N) operations).
    """
    N = len(d)
    c_ = [0.0] * N
    d_ = [0.0] * N
    c_[0] = c[0] / b[0]
    d_[0] = d[0] / b[0]
    for i in range(1, N):
        denom = b[i] - a[i] * c_[i - 1]
        if abs(denom) < 1e-30:
            denom = 1e-30
        c_[i] = c[i] / denom
        d_[i] = (d[i] - a[i] * d_[i - 1]) / denom

    x = [0.0] * N
    x[N - 1] = d_[N - 1]
    for i in range(N - 2, -1, -1):
        x[i] = d_[i] - c_[i] * x[i + 1]
    return x


class UHSRadial1D:
    """
    Implicit 1D radial finite-volume UHS simulator.
    Unconditionally stable backward-Euler solver with Thomas algorithm.
    """

    def __init__(
        self,
        props: RadialReservoirProps,
        N: int = 40,
        eos: PengRobinsonH2 | PengRobinsonMixture | None = None,
        max_picard: int = 20,
        picard_tol_Pa: float = 0.1,
    ) -> None:
        self.props = props
        self.N = N
        if eos is None:
            vt = getattr(props, "volume_translation", 0.0)
            self.eos: PengRobinsonH2 | PengRobinsonMixture = PengRobinsonH2(volume_translation=vt)
        else:
            self.eos = eos
        self.tau_s = props.tau_days * 86400.0
        self.max_picard = max_picard
        self.picard_tol_Pa = picard_tol_Pa
        self._last_picard_residual_Pa: float = 0.0

        self.r_c, self.r_i = build_log_grid(props.r_well, props.r_ext, N)
        self.V_cell = [
            math.pi * (self.r_i[i + 1] ** 2 - self.r_i[i] ** 2) * props.h for i in range(N)
        ]
        self.rel_perm: RelativePermeabilityModel = getattr(
            props, "rel_perm", None
        ) or RelativePermeabilityModel(
            swr=float(getattr(props, "swr", 0.0)),
            sgc=float(getattr(props, "sgc", 0.0)),
        )
        self._phi_eff: float = props.phi * self.rel_perm.effective_porosity_factor()
        self._k_eff: float = props.k * self.rel_perm.effective_permeability_factor()
        self.Vp = [V * self._phi_eff for V in self.V_cell]

        cg = (getattr(props, "cushion_gas", "H2") or "H2").strip().upper()
        if cg not in GAS_COMPONENTS:
            raise ValueError(f"Unknown cushion_gas {cg!r}. Choose from {list(GAS_COMPONENTS)}")
        self._cushion_gas: str = cg
        y_init_raw = getattr(props, "initial_h2_mole_fraction", None)
        self._y_init: float = (
            min(max(float(y_init_raw), 0.0), 1.0)
            if y_init_raw is not None
            else (1.0 if cg == "H2" else 0.0)
        )
        self._is_multicomponent: bool = not (cg == "H2" and self._y_init >= 1.0 - 1e-12)
        vt = getattr(props, "volume_translation", 0.0)
        self._mix_eos = PengRobinsonMixture(
            composition={"H2": self._y_init, cg: 1.0 - self._y_init} if cg != "H2" else {"H2": 1.0},
            peneloux=getattr(props, "peneloux", False),
            volume_translation=vt if vt != 0.0 and not getattr(props, "peneloux", False) else None,
        )
        self.y_h2: list[float] = [self._y_init] * N
        self._last_y_new: list[float] = list(self.y_h2)

    def _comp_at(self, y_val: float) -> dict[str, float]:
        if self._cushion_gas == "H2" or y_val >= 1.0 - 1e-12:
            return {"H2": 1.0}
        if y_val <= 1e-12:
            return {self._cushion_gas: 1.0}
        return {"H2": y_val, self._cushion_gas: 1.0 - y_val}

    def _cell_density(self, P_val: float, T_val: float, y_val: float = 1.0) -> float:
        if not self._is_multicomponent:
            return self.eos.density(P_val, T_val)
        return self._mix_eos.density(P_val, T_val, self._comp_at(y_val))

    def _cell_drho_dP(self, P_val: float, T_val: float, y_val: float = 1.0) -> float:
        if not self._is_multicomponent:
            return self.eos.drho_dP(P_val, T_val)
        return self._mix_eos.drho_dP(P_val, T_val, self._comp_at(y_val))

    def _cell_viscosity(
        self, T_val: float, P_val: float, rho_val: float, y_val: float = 1.0
    ) -> float:
        if not self._is_multicomponent:
            return viscosity_h2(T_val, P_val, rho=rho_val, eos=self.eos)
        return self._mix_eos.viscosity(P_val, T_val, rho=rho_val, composition=self._comp_at(y_val))

    def _trans(
        self,
        P: list[float],
        T: list[float],
        i: int,
        y_prof: list[float] | None = None,
    ) -> float:
        # i = N-1 is the link to the ghost cell used by the const-pressure
        # outer boundary; P and T must then contain N+1 entries.
        if i < 0 or i > self.N - 1:
            return 0.0
        y_arr = y_prof if y_prof is not None else self.y_h2
        y_L = y_arr[i] if i < len(y_arr) else self._y_init
        y_R = y_arr[i + 1] if i + 1 < len(y_arr) else self._y_init
        rho_L = self._cell_density(P[i], T[i], y_L)
        rho_R = self._cell_density(P[i + 1], T[i + 1], y_R)
        mu_L = self._cell_viscosity(T[i], P[i], rho_L, y_L)
        mu_R = self._cell_viscosity(T[i + 1], P[i + 1], rho_R, y_R)
        if (rho_L + rho_R) < 1e-20 or (mu_L + mu_R) < 1e-30:
            return 0.0
        rho_h = 2.0 * rho_L * rho_R / (rho_L + rho_R)
        mu_h = 2.0 * mu_L * mu_R / (mu_L + mu_R)
        j = min(i + 1, self.N - 1)
        ln_r = math.log(self.r_i[j + 1] / self.r_i[j])
        if ln_r < 1e-12:
            return 0.0
        return 2.0 * math.pi * self._k_eff * self.props.h * rho_h / (mu_h * ln_r)

    def _well_index(self, T_w: float, P_w: float, y_w: float | None = None) -> float:
        # Peaceman radial equivalent radius: the innermost cell centre is the
        # well block radius, which removes the arbitrary 0.2·Δr approximation.
        r_eq = self.r_c[0]
        yw = self.y_h2[0] if y_w is None else y_w
        rho_w = self._cell_density(P_w, T_w, yw)
        mu_w = self._cell_viscosity(T_w, P_w, rho_w, yw)
        ln_r = math.log(max(r_eq, self.props.r_well * 1.001) / self.props.r_well)
        ln_r = max(ln_r + self.props.skin, 0.1)
        return 2.0 * math.pi * self._k_eff * self.props.h * rho_w / (mu_w * ln_r)

    def _solve_transport_step(
        self,
        rho_old: list[float],
        Pk: list[float],
        Tk: list[float],
        T_half: list[float],
        q_well: float,
        dt: float,
        y_old: list[float],
    ) -> list[float]:
        """
        Implicit backward-Euler upwind advection-dispersion solver for radial
        hydrogen mass fraction w_H2, returning updated mole fractions y_H2.
        """
        N = self.N
        M_cg = GAS_COMPONENTS[self._cushion_gas].molar_mass
        w_old = [(y * M_H2) / max(y * M_H2 + (1.0 - y) * M_cg, 1e-30) for y in y_old]
        m_cell_old = [max(self.Vp[i] * rho_old[i], 1e-12) for i in range(N)]

        # Inter-cell mass flow rate F_i [kg/s] from cell i to cell i+1
        F = [T_half[i] * (Pk[i] - Pk[i + 1]) for i in range(N - 1)]

        # Radial dispersive mass conductance G_i [kg/s] across interface i+1/2
        alpha_L = max(float(getattr(self.props, "dispersivity_m", 0.5)), 0.0)
        D_mol = max(float(getattr(self.props, "molecular_diffusivity_m2_s", 1.0e-5)), 0.0)
        G = [0.0] * (N - 1)
        for i in range(N - 1):
            r_int = self.r_i[i + 1]
            area_int = 2.0 * math.pi * r_int * self.props.h * self._phi_eff
            rho_avg = 0.5 * (
                self._cell_density(Pk[i], Tk[i], y_old[i])
                + self._cell_density(Pk[i + 1], Tk[i + 1], y_old[i + 1])
            )
            v_int = abs(F[i]) / max(area_int * rho_avg, 1e-12)
            D_eff = D_mol + alpha_L * v_int
            ln_rc = math.log(max(self.r_c[i + 1] / self.r_c[i], 1.0001))
            G[i] = 2.0 * math.pi * self.props.h * self._phi_eff * rho_avg * D_eff / ln_rc

        a = [0.0] * N
        b = [0.0] * N
        c = [0.0] * N
        d = [0.0] * N

        q_inj = max(q_well, 0.0)
        for i in range(N):
            a_i = -(max(F[i - 1], 0.0) + G[i - 1]) if i > 0 else 0.0
            c_i = -(max(-F[i], 0.0) + G[i]) if i < N - 1 else 0.0
            well_in = q_inj if i == 0 else 0.0
            base_m = m_cell_old[i] / dt
            a[i] = a_i
            c[i] = c_i
            b[i] = base_m + well_in - a_i - c_i
            d[i] = base_m * w_old[i] + well_in * 1.0

        w_new = thomas_solve(a, b, c, d)
        y_new: list[float] = []
        for w in w_new:
            wc = min(max(w, 0.0), 1.0)
            n_h2 = wc / M_H2
            n_cg = (1.0 - wc) / M_cg
            y_new.append(n_h2 / max(n_h2 + n_cg, 1e-30))
        return y_new

    def _solve_mass_step(
        self,
        P_old: list[float],
        T_old: list[float],
        ctrl: WellControl,
        dt: float,
        fixed_q: float | None,
    ) -> tuple[list[float], list[float], float]:
        """
        Solve one implicit step with a mass-conservative accumulation term.

        The accumulation is linearised around the current pressure iterate:
            Vp·(ρ(P) - ρ_old)/dt ≈ A·(P - P_iter) + Vp·(ρ_iter - ρ_old)/dt
        with A = Vp·(dρ/dP).  Picard iterations remove the linearisation error,
        so the discrete mass balance closes to the pressure tolerance.  The
        temperature update is part of the iterate because density depends on T;
        this keeps the well source and the stored mass consistent.

        ``fixed_q`` is a signed mass rate [kg/s] for rate-controlled periods;
        ``None`` means the well is controlled by its BHP target.
        """
        N = self.N
        y_old = list(self.y_h2)
        T_half = [self._trans(P_old, T_old, i, y_old) if i < N - 1 else 0.0 for i in range(N)]
        rho_old = [self._cell_density(P_old[i], T_old[i], y_old[i]) for i in range(N)]
        Pk = list(P_old)
        Tk = list(T_old)
        y_k = list(y_old)
        T_geo = self.props.T_init

        for _picard in range(self.max_picard):
            Tk = [
                coupled_temperature_update(
                    T_old=T_old[i],
                    P_old=max(P_old[i], 1.0),
                    P_new=max(Pk[i], 1.0),
                    T_geothermal=T_geo,
                    dt_s=dt,
                    tau_s=self.tau_s,
                    eos=self.eos,
                )
                for i in range(N)
            ]
            rho_k = [self._cell_density(Pk[i], Tk[i], y_k[i]) for i in range(N)]
            drho = [self._cell_drho_dP(Pk[i], Tk[i], y_k[i]) for i in range(N)]

            a = [0.0] * N
            b = [0.0] * N
            c = [0.0] * N
            d = [0.0] * N

            WI = 0.0
            if ctrl.mode != "shut_in" and fixed_q is None:
                WI = self._well_index(Tk[0], Pk[0], y_k[0])

            for i in range(N):
                T_left = T_half[i - 1] if i > 0 else 0.0
                T_right = T_half[i] if i < N - 1 else 0.0
                T_right_outer = 0.0
                if i == N - 1 and self.props.outer_bc == "const_pressure":
                    T_right_outer = self._trans(
                        Pk + [self.props.P_init],
                        Tk + [self.props.T_init],
                        N - 1,
                        y_k + [self._y_init],
                    )

                A = self.Vp[i] * drho[i] / dt
                # Exact accumulation linearised around Pk.
                accum = self.Vp[i] * (rho_k[i] - rho_old[i]) / dt

                a[i] = -T_left
                c[i] = -T_right
                b[i] = A + T_left + T_right + T_right_outer
                d[i] = A * Pk[i] - accum

                if i == N - 1 and T_right_outer > 0.0:
                    d[N - 1] += T_right_outer * self.props.P_init

            if ctrl.mode != "shut_in":
                if fixed_q is None:
                    b[0] += WI
                    d[0] += WI * ctrl.bhp_target
                else:
                    d[0] += fixed_q

            P_new = thomas_solve(a, b, c, d)
            P_new = [max(p, 1.0e3) for p in P_new]
            change = max(abs(P_new[i] - Pk[i]) for i in range(N))
            self._last_picard_residual_Pa = change
            Pk = P_new

            if self._is_multicomponent and _picard < 2:
                if ctrl.mode == "shut_in":
                    q_est = 0.0
                elif fixed_q is not None:
                    q_est = fixed_q
                else:
                    q_est = WI * (ctrl.bhp_target - Pk[0])
                y_next = self._solve_transport_step(rho_old, Pk, Tk, T_half, q_est, dt, y_old)
                y_diff = max(abs(y_next[i] - y_k[i]) for i in range(N))
                y_k = y_next
                if y_diff > 1e-9:
                    continue

            if change < self.picard_tol_Pa:
                break

        # Synchronize final temperature with the converged pressure Pk.
        Tk = [
            coupled_temperature_update(
                T_old=T_old[i],
                P_old=max(P_old[i], 1.0),
                P_new=max(Pk[i], 1.0),
                T_geothermal=T_geo,
                dt_s=dt,
                tau_s=self.tau_s,
                eos=self.eos,
            )
            for i in range(N)
        ]

        if ctrl.mode == "shut_in":
            q = 0.0
        elif fixed_q is not None:
            q = fixed_q
        else:
            WI = self._well_index(Tk[0], Pk[0], y_k[0])
            q = WI * (ctrl.bhp_target - Pk[0])

        self._last_y_new = y_k
        return Pk, Tk, q

    def _solve_pressure(
        self,
        P: list[float],
        T: list[float],
        ctrl: WellControl,
        dt: float,
    ) -> tuple[list[float], list[float], float, float]:
        """
        Solve one pressure/temperature step with an active-set well model.

        The well is first solved against its BHP target.  If the resulting rate
        exceeds ``q_max`` the well is switched to rate control and the pressure
        is re-solved with a fixed mass source.  The well-cell pressure is never
        overwritten, so the reported rate, BHP and reservoir mass are mutually
        consistent.
        """
        if ctrl.mode == "shut_in":
            P_new, T_new, _ = self._solve_mass_step(P, T, ctrl, dt, fixed_q=0.0)
            self.y_h2 = list(self._last_y_new)
            return P_new, T_new, 0.0, P_new[0]

        WI0 = self._well_index(T[0], P[0], self.y_h2[0])
        q_try = WI0 * (ctrl.bhp_target - P[0])
        # A well cannot flow backwards through its intended direction; shut it
        # in if the reservoir has already crossed the BHP target.
        if (ctrl.mode == "inject" and q_try < 0.0) or (ctrl.mode == "withdraw" and q_try > 0.0):
            P_new, T_new, _ = self._solve_mass_step(P, T, ctrl, dt, fixed_q=0.0)
            self.y_h2 = list(self._last_y_new)
            return P_new, T_new, 0.0, P_new[0]

        fixed_q: float | None = None
        if ctrl.q_max > 0.0 and abs(q_try) > ctrl.q_max:
            fixed_q = math.copysign(ctrl.q_max, q_try)

        P_new, T_new, q = self._solve_mass_step(P, T, ctrl, dt, fixed_q)

        if fixed_q is not None:
            # If the BHP constraint is already satisfied at the capped rate,
            # the rate cap is not binding and the BHP-controlled solution holds.
            binding = (ctrl.mode == "inject" and P_new[0] <= ctrl.bhp_target) or (
                ctrl.mode == "withdraw" and P_new[0] >= ctrl.bhp_target
            )
            if not binding:
                P_new, T_new, q = self._solve_mass_step(P, T, ctrl, dt, None)
        elif ctrl.q_max > 0.0 and abs(q) > ctrl.q_max * (1.0 + 1e-6):
            # Post-solve rate-cap check: if mobility changed during the step so
            # that BHP control exceeded q_max, switch to rate control.
            fixed_q = math.copysign(ctrl.q_max, q)
            P_new, T_new, q = self._solve_mass_step(P, T, ctrl, dt, fixed_q)
        elif (ctrl.mode == "inject" and q < 0.0) or (ctrl.mode == "withdraw" and q > 0.0):
            P_new, T_new, q = self._solve_mass_step(P, T, ctrl, dt, fixed_q=0.0)

        self.y_h2 = list(self._last_y_new)
        WI = self._well_index(T_new[0], P_new[0], self.y_h2[0])
        bhp = P_new[0] + q / max(WI, 1e-30)
        return P_new, T_new, q, bhp

    def run_schedule(
        self,
        schedule: list[WellControl],
        dt_days: float = 0.1,
        snap_days: list[float] | None = None,
    ) -> RadialResult:
        P = [self.props.P_init] * self.N
        T = [self.props.T_init] * self.N
        self.y_h2 = [self._y_init] * self.N
        self._last_y_new = list(self.y_h2)
        result = RadialResult()

        cum_inj = 0.0
        cum_wdr = 0.0
        t = 0.0
        snap_targets = sorted(snap_days or [])
        snap_idx = 0
        m_initial = sum(
            self._cell_density(P[i], T[i], self.y_h2[i]) * self.Vp[i] for i in range(self.N)
        )
        m_max = m_initial
        m_min = m_initial
        max_picard_res = 0.0
        min_wdr_y_h2 = self._y_init

        for ctrl in schedule:
            n_steps = max(round(ctrl.duration_d / dt_days), 1)
            dt_act = ctrl.duration_d * 86400.0 / n_steps
            for _ in range(n_steps):
                P_new, T_new, q_well, BHP = self._solve_pressure(P, T, ctrl, dt_act)
                max_picard_res = max(max_picard_res, self._last_picard_residual_Pa)
                t += dt_act

                if q_well >= 0.0:
                    cum_inj += q_well * dt_act
                else:
                    cum_wdr += (-q_well) * dt_act

                P = P_new
                T = T_new

                m_now = sum(
                    self._cell_density(P[i], T[i], self.y_h2[i]) * self.Vp[i] for i in range(self.N)
                )
                m_max = max(m_max, m_now)
                m_min = min(m_min, m_now)

                y_well = self.y_h2[0]
                wobbe_well = self._mix_eos.wobbe_index(self._comp_at(y_well))
                if q_well < -1e-9:
                    min_wdr_y_h2 = min(min_wdr_y_h2, y_well)

                result.times_s.append(t)
                result.BHP_Pa.append(BHP)
                result.q_kg_s.append(q_well)
                result.cum_inj_kg.append(cum_inj)
                result.cum_wdr_kg.append(cum_wdr)
                result.P_avg_Pa.append(sum(P) / len(P))
                result.T_avg_K.append(sum(T) / len(T))
                result.mass_kg.append(m_now)
                result.well_h2_mole_fraction.append(y_well)
                result.well_wobbe_index.append(wobbe_well)

                t_days = t / 86400.0
                while snap_idx < len(snap_targets) and t_days >= snap_targets[snap_idx] - 1e-9:
                    result.snapshots.append(
                        RadialSnapshot(
                            time_s=t,
                            r_m=list(self.r_c),
                            P_Pa=list(P),
                            T_K=list(T),
                            h2_mole_fraction=list(self.y_h2),
                        )
                    )
                    snap_idx += 1

        total_Vp = sum(self.Vp)
        m_final = result.mass_kg[-1] if result.mass_kg else m_initial
        net = cum_inj - cum_wdr
        mb_err = abs((m_final - m_initial) - net) / max(abs(net), 1.0)
        result.metrics = {
            "total_injected_kg": cum_inj,
            "total_withdrawn_kg": cum_wdr,
            "net_stored_kg": net,
            "cycled_gas_kg": cum_inj + cum_wdr,
            "initial_mass_kg": m_initial,
            "final_mass_kg": m_final,
            "mass_balance_error": mb_err,
            "working_gas_kg": m_max - m_min,
            "cushion_gas_kg": m_min,
            "pore_volume_m3": total_Vp,
            "final_BHP_Pa": result.BHP_Pa[-1] if result.BHP_Pa else P[0],
            "final_avg_P_Pa": result.P_avg_Pa[-1] if result.P_avg_Pa else P[0],
            "max_picard_residual_Pa": max_picard_res,
            "converged": 1.0 if max_picard_res < self.picard_tol_Pa else 0.0,
            "initial_h2_purity_mol_pct": self._y_init * 100.0,
            "final_well_h2_purity_mol_pct": (
                result.well_h2_mole_fraction[-1] if result.well_h2_mole_fraction else self._y_init
            )
            * 100.0,
            "min_withdrawn_h2_purity_mol_pct": min_wdr_y_h2 * 100.0,
            "final_well_wobbe_index_MJ_sm3": (
                result.well_wobbe_index[-1]
                if result.well_wobbe_index
                else self._mix_eos.wobbe_index(self._comp_at(self._y_init))
            ),
        }
        return result


H2Reservoir1DModel = UHSRadial1D


def okoroafor_radial() -> tuple[RadialReservoirProps, list[WellControl]]:
    """Canonical 1D radial scenario matching Okoroafor et al. (2022)."""
    props = RadialReservoirProps(
        r_well=0.10,
        r_ext=500.0,
        h=100.0,
        phi=0.20,
        k=mD_to_m2(100.0),
        Ct=5.0e-10,
        P_init=bar_to_pa(80.0),
        T_init=316.0,
        skin=0.0,
        outer_bc="no_flow",
        tau_days=30.0,
    )
    sched = [
        WellControl("cushion_fill", "inject", 90.0, bar_to_pa(108.0), 35.0),
        WellControl("shut_in_1", "shut_in", 30.0, 0.0, 0.0),
        WellControl("withdrawal_1", "withdraw", 60.0, bar_to_pa(55.0), 35.0),
        WellControl("shut_in_2", "shut_in", 30.0, 0.0, 0.0),
    ]
    return props, sched


def scenario_1d_nearwell() -> tuple[RadialReservoirProps, list[WellControl]]:
    return okoroafor_radial()


def rough_radial() -> tuple[RadialReservoirProps, list[WellControl]]:
    """
    Rough-style 1D radial scenario. 150 bar / 369 K, 2700 m depth.
    90-day cushion fill -> 60-day withdrawal -> 10-day shut-in.
    """
    props = RadialReservoirProps(
        r_well=0.15,
        r_ext=500.0,
        h=30.0,
        phi=0.15,
        k=mD_to_m2(100.0),
        Ct=4.0e-10,
        P_init=bar_to_pa(150.0),
        T_init=369.0,
        skin=0.0,
        outer_bc="no_flow",
        tau_days=40.0,
    )
    schedule = [
        WellControl("cushion_90d", "inject", 90.0, bar_to_pa(217.5), 90.0),
        WellControl("withdraw_60d", "withdraw", 60.0, bar_to_pa(97.5), 90.0),
        WellControl("shut_in_10d", "shut_in", 10.0, 0.0, 0.0),
    ]
    return props, schedule


def export_radial_timeseries_csv(result: RadialResult, filepath: str) -> None:
    """Export 1D radial time series to CSV."""
    if not result.times_s:
        return
    os.makedirs(os.path.dirname(os.path.abspath(filepath)), exist_ok=True)
    with open(filepath, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(
            [
                "time_days",
                "bhp_bar",
                "bhp_psi",
                "q_kg_s",
                "cum_inj_tonnes",
                "cum_wdr_tonnes",
                "p_avg_bar",
                "t_avg_C",
            ]
        )
        for t, bhp, q, ci, cw, pa, ta in zip(
            result.times_s,
            result.BHP_Pa,
            result.q_kg_s,
            result.cum_inj_kg,
            result.cum_wdr_kg,
            result.P_avg_Pa,
            result.T_avg_K,
        ):
            w.writerow(
                [
                    f"{t / 86400:.4f}",
                    f"{pa_to_bar(bhp):.4f}",
                    f"{pa_to_psi(bhp):.2f}",
                    f"{q:.4f}",
                    f"{kg_to_tonnes(ci):.4f}",
                    f"{kg_to_tonnes(cw):.4f}",
                    f"{pa_to_bar(pa):.4f}",
                    f"{K_to_C(ta):.2f}",
                ]
            )
