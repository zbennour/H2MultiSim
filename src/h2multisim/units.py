"""
h2multisim.units
================
Unit conversions, dimension constants, and pure-Python numerical utilities.
Standard library only — zero external dependencies.
"""

from __future__ import annotations

import math
import random
from typing import Callable


# ── Pressure conversions ─────────────────────────────────────────────────────
def pa_to_bar(P: float) -> float:
    """Convert Pascal to bar."""
    return P / 1e5


def bar_to_pa(P: float) -> float:
    """Convert bar to Pascal."""
    return P * 1e5


def pa_to_MPa(P: float) -> float:
    """Convert Pascal to MegaPascal."""
    return P / 1e6


def MPa_to_pa(P: float) -> float:
    """Convert MegaPascal to Pascal."""
    return P * 1e6


def pa_to_psi(P: float) -> float:
    """Convert Pascal to pounds per square inch (psi)."""
    return P / 6894.757293168


def psi_to_pa(P: float) -> float:
    """Convert pounds per square inch (psi) to Pascal."""
    return P * 6894.757293168


# ── Temperature conversions ──────────────────────────────────────────────────
def K_to_C(T: float) -> float:
    """Convert Kelvin to degrees Celsius."""
    return T - 273.15


def C_to_K(T: float) -> float:
    """Convert degrees Celsius to Kelvin."""
    return T + 273.15


def K_to_F(T: float) -> float:
    """Convert Kelvin to degrees Fahrenheit."""
    return (T - 273.15) * 9.0 / 5.0 + 32.0


def F_to_K(T: float) -> float:
    """Convert degrees Fahrenheit to Kelvin."""
    return (T - 32.0) * 5.0 / 9.0 + 273.15


# ── Mass conversions ─────────────────────────────────────────────────────────
def kg_to_tonnes(m: float) -> float:
    """Convert kilograms to metric tonnes."""
    return m / 1000.0


def tonnes_to_kg(m: float) -> float:
    """Convert metric tonnes to kilograms."""
    return m * 1000.0


def kg_to_short_tons(m: float) -> float:
    """Convert kilograms to US short tons (2,000 lb)."""
    return m / 907.18474


def short_tons_to_kg(m: float) -> float:
    """Convert US short tons to kilograms."""
    return m * 907.18474


def kg_to_lb(m: float) -> float:
    """Convert kilograms to pounds (lb)."""
    return m * 2.2046226218487757


def lb_to_kg(m: float) -> float:
    """Convert pounds (lb) to kilograms."""
    return m / 2.2046226218487757


# ── Volume & Permeability conversions ────────────────────────────────────────
def m3_to_cubic_feet(V: float) -> float:
    """Convert cubic metres to cubic feet."""
    return V * 35.31466672148859


def cubic_feet_to_m3(V: float) -> float:
    """Convert cubic feet to cubic metres."""
    return V / 35.31466672148859


def millidarcy_to_m2(k: float) -> float:
    """Convert permeability from milliDarcy (mD) to square metres (m²)."""
    return k * 9.869233e-16


def m2_to_millidarcy(k: float) -> float:
    """Convert permeability from square metres (m²) to milliDarcy (mD)."""
    return k / 9.869233e-16


# Standard shorthand aliases
mD_to_m2 = millidarcy_to_m2
m2_to_mD = m2_to_millidarcy
pa_to_mpa = pa_to_MPa
mpa_to_pa = MPa_to_pa
k_to_c = K_to_C
c_to_k = C_to_K
k_to_f = K_to_F
f_to_k = F_to_K


# ── Pure-Python Numerical Utilities ──────────────────────────────────────────
def cbrt(x: float) -> float:
    """Cube root handling negative arguments accurately."""
    if x >= 0.0:
        return x ** (1.0 / 3.0)
    return -((-x) ** (1.0 / 3.0))


def cubic_real_roots(p: float, q: float, r: float) -> list[float]:
    """
    Solve monic cubic polynomial: x³ + p·x² + q·x + r = 0.

    Uses analytical Cardano / Viète trigonometric method.
    Returns a sorted list of all real roots without iteration.
    """
    p3 = p / 3.0
    A = q - p3 * p
    B = r + p3 * (2.0 * p3**2 - q)
    D = (B / 2.0) ** 2 + (A / 3.0) ** 3  # Discriminant

    if D > 1e-14:
        sq = math.sqrt(D)
        u = cbrt(-B / 2.0 + sq)
        v = cbrt(-B / 2.0 - sq)
        return [u + v - p3]

    if D < -1e-14:
        m = 2.0 * math.sqrt(max(-A / 3.0, 0.0))
        denom = A * m
        arg = (3.0 * B / denom) if abs(denom) > 1e-30 else 0.0
        arg = max(-1.0, min(1.0, arg))
        theta = math.acos(arg) / 3.0
        return sorted([m * math.cos(theta - 2.0 * math.pi * k / 3.0) - p3 for k in range(3)])

    u = cbrt(-B / 2.0)
    return sorted([2.0 * u - p3, -u - p3])


def runge_kutta_4(
    func: Callable[[float, list[float]], list[float]], t: float, y: list[float], dt: float
) -> list[float]:
    """
    One step of classic 4th-order Runge-Kutta integrator.
    func(t, y) -> dy/dt
    """
    k1 = func(t, y)
    k2 = func(t + dt / 2.0, [y[i] + dt / 2.0 * k1[i] for i in range(len(y))])
    k3 = func(t + dt / 2.0, [y[i] + dt / 2.0 * k2[i] for i in range(len(y))])
    k4 = func(t + dt, [y[i] + dt * k3[i] for i in range(len(y))])
    return [y[i] + dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) for i in range(len(y))]


def box_muller(rng: random.Random | None = None) -> float:
    """Generate standard normal variate N(0, 1) using Box-Muller transform."""
    r = rng if rng is not None else random
    u1 = max(r.random(), 1e-15)
    u2 = r.random()
    return math.sqrt(-2.0 * math.log(u1)) * math.cos(2.0 * math.pi * u2)


def lognormal_sample(mean: float, cv: float, rng: random.Random | None = None) -> float:
    """
    Sample from lognormal distribution with specified mean and coefficient of variation (cv = std/mean).
    Ensures strictly positive physical parameters (permeability, porosity, thickness).
    """
    sigma = math.sqrt(math.log(1.0 + cv**2))
    mu = math.log(mean) - 0.5 * sigma**2
    z = box_muller(rng)
    return math.exp(mu + sigma * z)


def percentile(sorted_values: list[float], p: float) -> float:
    """Compute percentile p in [0, 100] from pre-sorted list of floats."""
    if not sorted_values:
        return float("nan")
    p = max(0.0, min(100.0, p))
    idx = (len(sorted_values) - 1) * p / 100.0
    lo = int(idx)
    hi = min(lo + 1, len(sorted_values) - 1)
    return sorted_values[lo] + (sorted_values[hi] - sorted_values[lo]) * (idx - lo)
