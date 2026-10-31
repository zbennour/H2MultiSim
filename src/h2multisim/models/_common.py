"""Shared helpers for the h2multisim solver tiers."""

from __future__ import annotations

import math

WELL_MODES: tuple[str, ...] = ("inject", "withdraw", "shut_in")


def validate_mode(mode: str, where: str) -> str:
    """Return ``mode`` if it is a known well mode, otherwise raise ``ValueError``."""
    if mode not in WELL_MODES:
        raise ValueError(
            f"{where}: unknown well mode {mode!r}. Choose from {', '.join(WELL_MODES)}"
        )
    return mode


def check_inputs(
    where: str,
    *,
    positive: dict[str, float] | None = None,
    non_negative: dict[str, float] | None = None,
    fraction: dict[str, float] | None = None,
) -> None:
    """
    Validate physical inputs, raising ``ValueError`` with the offending field.

    ``positive``: must be finite and > 0; ``non_negative``: finite and >= 0;
    ``fraction``: finite and in (0, 1].
    """
    for name, val in (positive or {}).items():
        if not (math.isfinite(val) and val > 0.0):
            raise ValueError(f"{where}: {name} must be > 0, got {val!r}")
    for name, val in (non_negative or {}).items():
        if not (math.isfinite(val) and val >= 0.0):
            raise ValueError(f"{where}: {name} must be >= 0, got {val!r}")
    for name, val in (fraction or {}).items():
        if not (math.isfinite(val) and 0.0 < val <= 1.0):
            raise ValueError(f"{where}: {name} must be in (0, 1], got {val!r}")


def pcg_solve(
    diag: list[float],
    nbrs: list[list[tuple[int, float]]],
    b: list[float],
    x0: list[float],
    tol: float,
    max_iter: int,
) -> tuple[list[float], bool]:
    """
    Jacobi-preconditioned conjugate gradients for the symmetric positive-definite
    finite-volume pressure system ``diag[i]·x[i] − Σ_(j,T) T·x[j] = b[i]``.

    Converged when the Jacobi-scaled residual ``max |r_i| / diag_i`` (the pressure
    correction a point-Jacobi sweep would make, in the units of ``x``) is below
    ``tol``.  Returns ``(x, converged)``.
    """
    n = len(b)
    x = list(x0)
    inv_d = [1.0 / max(d, 1e-300) for d in diag]

    def matvec(v: list[float]) -> list[float]:
        out = [0.0] * n
        for i in range(n):
            acc = diag[i] * v[i]
            for j, t in nbrs[i]:
                acc -= t * v[j]
            out[i] = acc
        return out

    ax = matvec(x)
    r = [b[i] - ax[i] for i in range(n)]
    if max(abs(r[i]) * inv_d[i] for i in range(n)) < tol:
        return x, True
    z = [r[i] * inv_d[i] for i in range(n)]
    p = list(z)
    rz = sum(r[i] * z[i] for i in range(n))
    for _ in range(max_iter):
        ap = matvec(p)
        pap = sum(p[i] * ap[i] for i in range(n))
        if pap <= 0.0 or not math.isfinite(pap):
            return x, False
        alpha = rz / pap
        for i in range(n):
            x[i] += alpha * p[i]
            r[i] -= alpha * ap[i]
        if max(abs(r[i]) * inv_d[i] for i in range(n)) < tol:
            return x, True
        z = [r[i] * inv_d[i] for i in range(n)]
        rz_new = sum(r[i] * z[i] for i in range(n))
        beta = rz_new / rz
        rz = rz_new
        p = [z[i] + beta * p[i] for i in range(n)]
    return x, False


class AndersonAccelerator:
    """
    Anderson acceleration (type II, depth ``m``) for a fixed-point map x = g(x)
    on flat float lists.  ``step(x, gx)`` returns the next iterate.  Used to
    speed up the linearly convergent pressure–composition Picard coupling.
    """

    def __init__(self, m: int = 3) -> None:
        self.m = m
        self._dg: list[list[float]] = []  # differences of g(x)
        self._df: list[list[float]] = []  # differences of residual f = g(x) - x
        self._prev: tuple[list[float], list[float]] | None = None

    def step(self, x: list[float], gx: list[float]) -> list[float]:
        f = [g - xi for g, xi in zip(gx, x)]
        if self._prev is not None:
            g_prev, f_prev = self._prev
            self._dg.append([a - b for a, b in zip(gx, g_prev)])
            self._df.append([a - b for a, b in zip(f, f_prev)])
            if len(self._df) > self.m:
                self._dg.pop(0)
                self._df.pop(0)
        self._prev = (gx, f)
        k = len(self._df)
        if k == 0:
            return list(gx)
        # Least squares min ||f − ΔF·γ|| via the k×k normal equations (k <= m, small).
        A = [
            [sum(a * b for a, b in zip(self._df[r], self._df[c])) for c in range(k)]
            for r in range(k)
        ]
        rhs = [sum(a * b for a, b in zip(self._df[r], f)) for r in range(k)]
        for r in range(k):
            A[r][r] *= 1.0 + 1e-10
        gamma = _solve_small(A, rhs)
        if gamma is None:
            return list(gx)
        out = list(gx)
        for c in range(k):
            gc = gamma[c]
            dgc = self._dg[c]
            for i in range(len(out)):
                out[i] -= gc * dgc[i]
        return out


def _solve_small(A: list[list[float]], b: list[float]) -> list[float] | None:
    """Gaussian elimination with partial pivoting for a small dense system."""
    n = len(b)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(M[r][col]))
        if abs(M[piv][col]) < 1e-300:
            return None
        M[col], M[piv] = M[piv], M[col]
        for r in range(col + 1, n):
            fac = M[r][col] / M[col][col]
            for c in range(col, n + 1):
                M[r][c] -= fac * M[col][c]
    x = [0.0] * n
    for r in range(n - 1, -1, -1):
        x[r] = (M[r][n] - sum(M[r][c] * x[c] for c in range(r + 1, n))) / M[r][r]
    return x
