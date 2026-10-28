"""Prototype fix: replace the liquid phi*Ct*Vp storage term with the EOS
gas storage Vp*d(rho)/dP and remove the non-physical cell-pressure clamp.
Compares against the mass-conservative reference at constant BHP."""

import math
import sys

import numpy as np

import pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
from h2multisim.models.one_d import (  # noqa: E402
    UHSRadial1D,
    WellControl,
    build_log_grid,
    okoroafor_radial,
    thomas_solve,
)
from h2multisim.thermo import PengRobinsonH2, viscosity_h2  # noqa: E402
from h2multisim.units import bar_to_pa  # noqa: E402

EOS = PengRobinsonH2()


class Fixed1D(UHSRadial1D):
    """Gas-storage (dV rho / dP) accumulation, rate-consistent well handling."""

    def _solve_pressure(self, P, T, ctrl, dt):
        N = self.N
        T_half = [self._trans(P, T, i) for i in range(N)]
        B, drho = [], []
        for i in range(N):
            dP = 1.0e3
            rho_i = self.eos.density(P[i], T[i])
            cg = (
                self.eos.density(P[i] + dP, T[i]) - self.eos.density(max(P[i] - dP, 1.0), T[i])
            ) / (2.0 * dP)
            B.append(self.Vp[i] * cg / dt)
            drho.append(rho_i)

        a, b, c, d = [0.0] * N, [0.0] * N, [0.0] * N, [0.0] * N
        Tr_outer = (
            self._trans(P + [self.props.P_init], T + [self.props.T_init], N - 1)
            if self.props.outer_bc == "const_pressure"
            else 0.0
        )
        for i in range(N):
            Tl = T_half[i - 1] if i > 0 else 0.0
            Tr = T_half[i] if i < N - 1 else 0.0
            a[i] = -Tl
            c[i] = -Tr
            b[i] = B[i] + Tl + Tr
            d[i] = B[i] * P[i]
        b[N - 1] += Tr_outer
        d[N - 1] += Tr_outer * self.props.P_init
        if ctrl.mode != "shut_in":
            WI = self._well_index(T[0], P[0])
            b[0] += WI
            d[0] += WI * ctrl.bhp_target
        Pn = [max(p, 1.0e3) for p in thomas_solve(a, b, c, d)]

        WI = self._well_index(T[0], Pn[0])
        if ctrl.mode == "inject":
            q = min(WI * max(ctrl.bhp_target - Pn[0], 0.0), ctrl.q_max)
        elif ctrl.mode == "withdraw":
            q = -min(WI * max(Pn[0] - ctrl.bhp_target, 0.0), ctrl.q_max)
        else:
            q = 0.0
        return Pn, q, Pn[0]


props, _ = okoroafor_radial()
props.tau_days = 0.0
N = 40
BHP = bar_to_pa(60.0)
QMAX = 1.0e9
dt_days = 0.2
days = 10.0

for cls, label in ((UHSRadial1D, "shipped"), (Fixed1D, "fixed  ")):
    sim = cls(props, N=N)
    res = sim.run_schedule([WellControl("w", "withdraw", days, BHP, QMAX)], dt_days=dt_days,
                           snap_days=[days])
    m0 = sum(EOS.density(props.P_init, props.T_init) * v for v in sim.Vp)
    m1 = res.metrics["final_mass_kg"]
    P = np.array(res.snapshots[-1].P_Pa)
    print(f"{label}: cum={res.metrics['total_withdrawn_kg']/1e3:9.2f} t  dM={(m1-m0)/1e3:9.2f} t  "
          f"MBerr={((m1-m0)+res.metrics['total_withdrawn_kg'])/max(res.metrics['total_withdrawn_kg'],1e-9)*100:8.2f}%  "
          f"P0={P[0]/1e5:6.2f}  avgP={P.mean()/1e5:6.2f}")

# reference (same WI)
r_c, r_i = build_log_grid(props.r_well, props.r_ext, N)
V = np.array([math.pi * (r_i[i + 1] ** 2 - r_i[i] ** 2) * props.h for i in range(N)])
Vp = props.phi * V
P = np.full(N, props.P_init)
rho = np.array([EOS.density(p, props.T_init) for p in P])
sim0 = UHSRadial1D(props, N=N)
WI = sim0._well_index(props.T_init, props.P_init)
dt = days * 86400 / 100
removed = 0.0


def trans(Pk):
    rk = np.array([EOS.density(p, props.T_init) for p in Pk])
    mk = np.array([viscosity_h2(props.T_init, p) for p in Pk])
    Tf = np.zeros(N - 1)
    for i in range(N - 1):
        rh = 2 * rk[i] * rk[i + 1] / (rk[i] + rk[i + 1])
        mh = 2 * mk[i] * mk[i + 1] / (mk[i] + mk[i + 1])
        Tf[i] = 2 * math.pi * props.k * props.h * rh / (mh * math.log(r_i[i + 2] / r_i[i + 1]))
    return Tf


for _ in range(100):
    Pk = P.copy()
    Tr = trans(Pk)
    for _it in range(80):
        rk = np.array([EOS.density(p, props.T_init) for p in Pk])
        cp = np.array([(EOS.density(p + 1e3, props.T_init) - EOS.density(p - 1e3, props.T_init)) / 2e3
                       for p in Pk])
        A = Vp * cp / dt
        diag = A.copy()
        diag[:-1] += Tr
        diag[1:] += Tr
        diag[0] += WI
        rhs = A * Pk - (rk - rho) * Vp / dt
        rhs[0] += WI * BHP
        a = np.zeros(N); b = diag.copy(); c = np.zeros(N)
        a[1:] = -Tr; c[:-1] = -Tr
        cc = np.zeros(N); dd = np.zeros(N)
        cc[0] = c[0] / b[0]; dd[0] = rhs[0] / b[0]
        for i in range(1, N):
            den = b[i] - a[i] * cc[i - 1]
            cc[i] = c[i] / den
            dd[i] = (rhs[i] - a[i] * dd[i - 1]) / den
        x = np.zeros(N); x[-1] = dd[-1]
        for i in range(N - 2, -1, -1):
            x[i] = dd[i] - cc[i] * x[i + 1]
        ch = np.max(np.abs(x - Pk))
        Pk = 0.5 * Pk + 0.5 * x
        if ch < 1e-2:
            break
    q = WI * (Pk[0] - BHP)
    removed += q * dt
    P = Pk
    rho = np.array([EOS.density(p, props.T_init) for p in P])

print(f"ref    : cum={removed/1e3:9.2f} t  P0={P[0]/1e5:6.2f}  avgP={P.mean()/1e5:6.2f}")
