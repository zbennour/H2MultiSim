"""
h2multisim.verify
=================
Self-verification of the installed thermodynamic core against reference values
tabulated from CoolProp 8.0.0 (Leachman et al. 2009 normal-H2 EOS; reference
EOS for CH4/CO2/N2).  No runtime dependency on CoolProp.

Each check compares one property with the acceptance tolerance documented in
LIMITATIONS.md.  ``h2multisim verify`` prints the table and exits non-zero if
any check fails.
"""

from __future__ import annotations

from dataclasses import dataclass

from h2multisim.thermo import (
    PengRobinsonH2,
    PengRobinsonMixture,
    adiabatic_temperature,
    cp_real_h2,
    joule_thomson_coefficient,
    viscosity_h2,
)

# (T [K], P [bar], Z, viscosity [µPa·s], Cp [J/(mol·K)])  — CoolProp 8.0.0, normal H2
H2_STATE_REFS: list[tuple[float, float, float, float, float]] = [
    (250.0, 10.0, 1.00643, 7.888, 28.4058),
    (250.0, 100.0, 1.06726, 8.0569, 29.0317),
    (250.0, 200.0, 1.14047, 8.3668, 29.4872),
    (250.0, 300.0, 1.21728, 8.7554, 29.7757),
    (300.0, 10.0, 1.00586, 8.9439, 28.9022),
    (300.0, 100.0, 1.05985, 9.0576, 29.3269),
    (300.0, 200.0, 1.1225, 9.2823, 29.6585),
    (300.0, 300.0, 1.18696, 9.5768, 29.8859),
    (350.0, 10.0, 1.00528, 9.9477, 29.1248),
    (350.0, 100.0, 1.05338, 10.0259, 29.4274),
    (350.0, 200.0, 1.10808, 10.1952, 29.6759),
    (350.0, 300.0, 1.16366, 10.4268, 29.8555),
    (400.0, 10.0, 1.00477, 10.91, 29.2137),
    (400.0, 100.0, 1.04789, 10.9639, 29.4368),
    (400.0, 200.0, 1.09634, 11.0949, 29.6272),
    (400.0, 300.0, 1.14519, 11.2822, 29.7703),
]

# (T0 [K], P0 [bar], P1 [bar], isentropic ΔT [K])
ISENTROPE_REFS: list[tuple[float, float, float, float]] = [
    (293.15, 40.0, 120.0, 108.392),
    (316.0, 80.0, 108.0, 28.396),
    (293.15, 120.0, 40.0, -81.159),
]

# (species, P [bar], density [kg/m³], viscosity [µPa·s]) at 320 K
SPECIES_REFS: list[tuple[str, float, float, float]] = [
    ("CH4", 100.0, 67.614, 14.126),
    ("CH4", 200.0, 138.234, 18.306),
    ("CO2", 100.0, 448.277, 31.758),
    ("CO2", 200.0, 802.331, 72.476),
    ("N2", 100.0, 103.649, 20.657),
    ("N2", 200.0, 197.001, 23.51),
]

JT_REF_K_PER_BAR: tuple[float, float, float] = (293.15, 100.0, -0.03466)

# Acceptance tolerances (relative), matching LIMITATIONS.md
TOL = {
    "Z": 0.0085,  # 280-450 K; 1.2 % below 280 K (constant Peneloux shift)
    "Z_cold": 0.012,
    "viscosity": 0.005,
    "Cp": 0.01,
    "isentropic_dT": 0.03,
    "JT_100bar": 0.10,
    "species_density": 0.045,
    "species_viscosity": 0.02,
}


@dataclass
class Check:
    name: str
    condition: str
    reference: float
    value: float
    tolerance: float

    @property
    def error(self) -> float:
        return (self.value - self.reference) / self.reference

    @property
    def passed(self) -> bool:
        return abs(self.error) <= self.tolerance


def run_checks() -> list[Check]:
    """Evaluate every reference check with the Peneloux-translated PR EOS."""
    eos = PengRobinsonH2(peneloux=True)
    checks: list[Check] = []
    for T, Pb, Z_ref, mu_ref, cp_ref in H2_STATE_REFS:
        P = Pb * 1e5
        cond = f"H2 {T:.0f} K, {Pb:.0f} bar"
        z_tol = TOL["Z"] if T >= 280.0 else TOL["Z_cold"]
        checks.append(Check("Z", cond, Z_ref, eos.z_factor(P, T), z_tol))
        checks.append(
            Check("viscosity", cond, mu_ref, viscosity_h2(T, P, eos=eos) * 1e6, TOL["viscosity"])
        )
        checks.append(Check("Cp", cond, cp_ref, cp_real_h2(T, P, eos), TOL["Cp"]))
    for T0, P0, P1, dT_ref in ISENTROPE_REFS:
        dT = adiabatic_temperature(T0, P0 * 1e5, P1 * 1e5, eos) - T0
        cond = f"H2 {T0:.0f} K, {P0:.0f}→{P1:.0f} bar"
        checks.append(Check("isentropic_dT", cond, dT_ref, dT, TOL["isentropic_dT"]))
    T, Pb, jt_ref = JT_REF_K_PER_BAR
    jt = joule_thomson_coefficient(T, Pb * 1e5) * 1e5
    checks.append(Check("JT_100bar", f"H2 {T:.0f} K, {Pb:.0f} bar", jt_ref, jt, TOL["JT_100bar"]))
    for sp, Pb, rho_ref, mu_ref in SPECIES_REFS:
        mix = PengRobinsonMixture({sp: 1.0}, peneloux=True)
        cond = f"{sp} 320 K, {Pb:.0f} bar"
        checks.append(
            Check(
                "species_density",
                cond,
                rho_ref,
                mix.density(Pb * 1e5, 320.0),
                TOL["species_density"],
            )
        )
        checks.append(
            Check(
                "species_viscosity",
                cond,
                mu_ref,
                mix.viscosity(Pb * 1e5, 320.0) * 1e6,
                TOL["species_viscosity"],
            )
        )
    return checks


def format_report(checks: list[Check]) -> str:
    """Human-readable verification table."""
    lines = [
        f"{'property':<18} {'condition':<28} {'reference':>11} {'H2MultiSim':>11} {'error %':>8} {'tol %':>6}  ",
        "-" * 92,
    ]
    for c in checks:
        lines.append(
            f"{c.name:<18} {c.condition:<28} {c.reference:>11.5g} {c.value:>11.5g} "
            f"{100 * c.error:>+8.2f} {100 * c.tolerance:>6.2f}  {'PASS' if c.passed else 'FAIL'}"
        )
    n_fail = sum(not c.passed for c in checks)
    lines.append("-" * 92)
    lines.append(f"{len(checks) - n_fail}/{len(checks)} checks passed (reference: CoolProp 8.0.0)")
    return "\n".join(lines)
