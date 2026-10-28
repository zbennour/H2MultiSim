"""Constant-BHP deliverability: shipped 1D model vs a mass-conservative
finite-volume reference on the same grid, same PR-EOS, same well index."""

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
)
from h2multisim.thermo import PengRobinsonH2, viscosity_h2  # noqa: E402
from h2multisim.units import bar_to_pa  # noqa: E402

EOS = PengRobinsonH2()
props, _ = okoroafor_radial()
props.tau_days = 0.0
N = 40
BHP = bar_to_pa(60.0)

sim = UHSRadial1D(props, N=N)
WI = sim._well_index(props.T_init, props.P_init)
print(f"well index used by code: {WI:.3e} kg/(Pa s)")

res = sim.run_schedule([WellControl("w", "withdraw", 10.0, BHP, 1.0e9)], dt_days=0.2,
                       snap_days=[10.0])
m0 = sum(EOS.density(props.P_init, props.T_init) * v for v in sim.Vp)
m1 = res.metrics["final_mass_kg"]
Pc = np.array(res.snapshots[-1].P_Pa)
print(f"shipped code: cum withdrawn = {res.metrics['total_withdrawn_kg']/1e3:10.2f} t, "
      f"P0 = {Pc[0]/1e5:6.2f} bar, avg P = {Pc.mean()/1e5:6.2f} bar, "
      f"dM = {(m1-m0)/1e3:9.2f} t, MB err = {((m1-m0)+res.metrics['total_withdrawn_kg'])/res.metrics['total_withdrawn_kg']*100:6.1f}%")

# ---------------- reference ----------------
r_c, r_i = build_log_grid(props.r_well, props.r_ext, N)
V = np.array([math.pi * (r_i[i + 1] ** 2 - r_i[i] ** 2) * props.h for i in range(N)])
Vp = props.phi * V
P = np.full(N, props.P_init)
rho = np.array([EOS.density(p, props.T_init) for p in P])
dt = 10 * 86400 / 100
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

print(f"reference   : cum withdrawn = {removed/1e3:10.2f} t, P0 = {P[0]/1e5:6.2f} bar, "
      f"avg P = {P.mean()/1e5:6.2f} bar")
print(f"ratio (reference / code) = {removed/max(res.metrics['total_withdrawn_kg'],1e-9):.1f}x")

print("\nprofile:")
print(f"{'r [m]':>8} {'ref [bar]':>10} {'code [bar]':>11}")
for i in (0, 2, 5, 10, 20, 30, 39):
    print(f"{r_c[i]:8.2f} {P[i]/1e5:10.2f} {Pc[i]/1e5:11.2f}")

# storage-coefficient ratio
rho0 = EOS.density(props.P_init, props.T_init)
cg = (EOS.density(props.P_init + 1e3, props.T_init) - EOS.density(props.P_init - 1e3, props.T_init)) / (2e3 * rho0)
print(f"\ntrue gas compressibility cg  = {cg:.3e} 1/Pa")
print(f"code storage term  phi^2*Ct = {props.phi**2*props.Ct:.3e} 1/Pa")
print(f"stiffness ratio (physical / code) = {props.phi*cg/(props.phi**2*props.Ct):.0f}x")
