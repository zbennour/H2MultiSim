"""Thermophysical benchmark: H2MultiSim vs NIST-reference (CoolProp HEOS)."""

import sys

import CoolProp.CoolProp as CP

import pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
from h2multisim.thermo import (  # noqa: E402
    PengRobinsonH2,
    cp_ideal_h2,
    gamma_real_h2,
    joule_thomson_coefficient,
    viscosity_h2,
)

FLUID = "Hydrogen"
eos = PengRobinsonH2()

def ref(prop, T, P):
    return CP.PropsSI(prop, "T", T, "P", P, FLUID)


def table(name, cases, fn_model, fn_ref, unit="", rel=True):
    print(f"\n## {name}")
    print(f"{'T[K]':>6} {'P[bar]':>7} {'model':>12} {'REF':>12} {'err%':>8}")
    errs = []
    for T, Pbar in cases:
        P = Pbar * 1e5
        try:
            m = fn_model(T, P)
            r = fn_ref(T, P)
        except Exception as exc:  # noqa: BLE001
            print(f"{T:6.1f} {Pbar:7.1f}  FAILED {exc}")
            continue
        e = (m - r) / r * 100.0 if rel else (m - r) * 100.0
        errs.append(abs(e))
        print(f"{T:6.1f} {Pbar:7.1f} {m:12.5g} {r:12.5g} {e:8.3f}")
    if errs:
        print(
            f"  -> mean|err| {sum(errs)/len(errs):.3f}%   max|err| {max(errs):.3f}%   (n={len(errs)})"
        )
    return errs


cases = [
    (273.15, 1), (293.15, 1), (293.15, 50), (293.15, 80), (293.15, 100),
    (293.15, 150), (293.15, 200), (293.15, 250), (293.15, 300), (293.15, 350),
    (316.0, 80), (316.0, 120), (333.0, 82), (350.0, 100), (350.0, 250),
    (400.0, 100), (400.0, 250), (400.0, 350), (500.0, 200), (700.0, 350),
]

print("=" * 78)
print("A. EOS / thermophysics vs CoolProp HEOS (Leachman 2009 reference EOS)")
print("=" * 78)

table("Z-factor", cases, lambda T, P: eos.z_factor(P, T), lambda T, P: ref("Z", T, P))
table("Density [kg/m3]", cases, lambda T, P: eos.density(P, T), lambda T, P: ref("Dmass", T, P))
table(
    "Ideal-gas cp (Shomate) [J/mol/K]",
    [(T, 1) for T in (250, 293.15, 300, 350, 400, 500, 600, 700, 800, 1000)],
    lambda T, P: cp_ideal_h2(T),
    lambda T, P: ref("Cpmass", T, P) * 2.01594e-3,
)
table(
    "Real-gas cp [J/kg/K] (not implemented; Cp0 only)",
    [(316.0, 80), (316.0, 120), (350.0, 250)],
    lambda T, P: cp_ideal_h2(T) / 2.01594e-3,
    lambda T, P: ref("Cpmass", T, P),
)
table(
    "Viscosity [uPa.s]",
    cases,
    lambda T, P: viscosity_h2(T, P) * 1e6,
    lambda T, P: ref("VISCOSITY", T, P) * 1e6,
)
table(
    "gamma = cp/cv [-]",
    [(316.0, 80), (316.0, 120), (350.0, 250), (293.15, 1), (400.0, 100)],
    lambda T, P: gamma_real_h2(T, P, eos),
    lambda T, P: ref("Cpmass", T, P) / ref("Cvmass", T, P),
)
table(
    "Joule-Thomson [K/bar]",
    [(293.15, 1), (293.15, 100), (316.0, 80), (350.0, 250), (400.0, 300)],
    lambda T, P: joule_thomson_coefficient(T, P, eos) * 1e5,
    lambda T, P: ref("Joule-Thomson coefficient", T, P) * 1e5,
)

print("\n\n## B. The repository's own 'NIST' Z table vs CoolProp")
repo = [
    (293.15, 101325.0, 1.00059), (293.15, 5.0e6, 1.02999), (293.15, 10.0e6, 1.06081),
    (293.15, 20.0e6, 1.12473), (293.15, 25.0e6, 1.15751), (350.00, 10.0e6, 1.05338),
    (350.00, 25.0e6, 1.13579), (400.00, 10.0e6, 1.04789), (400.00, 25.0e6, 1.12074),
    (316.00, 8.0e6, 1.04594),
]
print(f"{'T[K]':>6} {'P[bar]':>7} {'repo-NIST':>10} {'CoolProp':>10} {'repo err%':>10} {'model err%':>11}")
for T, P, Zr in repo:
    Zref = ref("Z", T, P)
    Zmod = eos.z_factor(P, T)
    print(
        f"{T:6.1f} {P/1e5:7.1f} {Zr:10.5f} {Zref:10.5f} "
        f"{(Zr-Zref)/Zref*100:10.3f} {(Zmod-Zref)/Zref*100:11.3f}"
    )
