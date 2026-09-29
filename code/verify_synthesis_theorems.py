#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
verify_synthesis_theorems.py
============================
Independent numerical checks of the statements that the integrating manuscript proves or recomputes itself
(the theorems of Papers I-VII are verified by the scripts archived with each paper).

    [1] Theorem 1.1(1)  h(delta) = 2 - sqrt(1 + 2 delta) - sqrt(1 - 2 delta) = delta^2 + 5/4 delta^4 + 21/8 delta^6 + ...
    [2] Theorem 1.1(2)  pair divergence D_k = (1/2) ln((1 - rho^2/4) / sqrt(1 - rho^2)) = rho^2/8 + 7 rho^4/64 + O(rho^6);
                        partial sums of D_k with rho_k = c k^-alpha diverge iff alpha <= 1/2
    [3] Proposition 1.2 pocket-local tilts: Bhattacharyya coefficient of the full ensembles equals that of the pocket
                        marginals (Gaussian example with a quadratic pocket energy)
    [4] Theorem 3.3     Szego rate s = 0.0094968218 by quadrature; D_cov(N) = s N + c_0 with c_0 = -0.0093715854 from exact
                        Toeplitz determinants; D_mean, D_cov and CKI at N = 501

Output: results/synthesis_verification.txt; exit code 0 if every check passes.
"""

from __future__ import annotations

import math
import os
import sys

import numpy as np
from scipy.integrate import quad
from scipy.linalg import toeplitz

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "results", "synthesis_verification.txt")
os.makedirs(os.path.dirname(OUT), exist_ok=True)
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

lines: list[str] = []
ok_all = True


def out(s: str = "") -> None:
    print(s)
    lines.append(s)


def check(name: str, cond: bool) -> None:
    global ok_all
    ok_all &= bool(cond)
    out(f"  {name}: {'PASS' if cond else 'FAIL'}")


def logdet(M: np.ndarray) -> float:
    s, v = np.linalg.slogdet(M)
    assert s > 0
    return float(v)


def bhattacharyya(mA, CA, mB, CB) -> float:
    Cb = 0.5 * (CA + CB)
    d = mB - mA
    return float(d @ np.linalg.solve(Cb, d) / 8 + 0.5 * logdet(Cb) - 0.25 * (logdet(CA) + logdet(CB)))


out("=" * 78)
out("Verification of the theorems proved or recomputed in the integrating manuscript")
out("=" * 78)

# ---------------- [1] ----------------
out("\n[1] Theorem 1.1(1): exact local index at p = 1/2 (50-digit arithmetic; the next coefficient is 429/64)")
import mpmath  # noqa: E402
mpmath.mp.dps = 50
ratios = []
for dl in ("0.001", "0.01", "0.05", "0.1"):
    d = mpmath.mpf(dl)
    exact = 2 - mpmath.sqrt(1 + 2 * d) - mpmath.sqrt(1 - 2 * d)
    series = d ** 2 + mpmath.mpf(5) / 4 * d ** 4 + mpmath.mpf(21) / 8 * d ** 6
    ratio = (exact - series) / d ** 8
    ratios.append(float(ratio))
    out(f"  delta = {dl:<6}: exact {mpmath.nstr(exact, 16)}, three-term series {mpmath.nstr(series, 16)}, "
        f"(exact - series) / delta^8 = {mpmath.nstr(ratio, 10)}")
check("remainder / delta^8 -> 429/64 = 6.703125 as delta -> 0", abs(ratios[0] - 429 / 64) < 1e-4)

# ---------------- [2] ----------------
out("\n[2] Theorem 1.1(2): divergence of a correlated pair with unit variances")
for rho in (0.05, 0.1, 0.3, 0.6):
    D = 0.5 * math.log((1 - rho ** 2 / 4) / math.sqrt(1 - rho ** 2))
    num = bhattacharyya(np.zeros(2), np.eye(2), np.zeros(2), np.array([[1, rho], [rho, 1]]))
    lead = rho ** 2 / 8
    out(f"  rho = {rho:<5}: closed form {D:.12f}, direct Bhattacharyya {num:.12f}, rho^2/8 {lead:.12f}, "
        f"rho^2/8 + 7 rho^4/64 {lead + 7 * rho ** 4 / 64:.12f}")
    check(f"closed form equals the direct computation (rho = {rho})", abs(D - num) < 1e-13)
c = 0.5
out("  partial sums S(K) = sum_{k<=K} D_k with rho_k = 0.5 k^-alpha:")
for alpha in (0.25, 0.5, 0.75):
    ks = np.arange(1, 10 ** 6 + 1, dtype=float)
    rho = c * ks ** -alpha
    Dk = 0.5 * np.log((1 - rho ** 2 / 4) / np.sqrt(1 - rho ** 2))
    S = np.cumsum(Dk)
    vals = [S[10 ** e - 1] for e in (2, 4, 6)]
    incr = vals[2] - vals[1]
    out(f"    alpha = {alpha}: S(1e2) = {vals[0]:.6f}, S(1e4) = {vals[1]:.6f}, S(1e6) = {vals[2]:.6f}; increment 1e4 -> 1e6 = {incr:.6f}")
    if alpha == 0.5:
        pred = c ** 2 / 8 * math.log(100)
        check(f"alpha = 1/2: logarithmic growth, increment {incr:.6f} vs (c^2/8) ln 100 = {pred:.6f}", abs(incr / pred - 1) < 1e-3)
    if alpha == 0.75:
        tail = c ** 2 / 8 * 2 * (1e4 ** -0.5)
        check(f"alpha = 3/4: convergent, increment {incr:.2e} vs tail bound (c^2/8) * 2 / sqrt(1e4) = {tail:.2e}", incr < 1.1 * tail)
    if alpha == 0.25:
        check("alpha = 1/4: power growth (S(1e6) / S(1e4) close to 10^(2 * 0.5) = 10)", 8 < vals[2] / vals[1] < 11)

# ---------------- [3] ----------------
out("\n[3] Proposition 1.2: pocket-local tilts of a common ensemble")
rng = np.random.default_rng(20260929)
n, F = 12, [0, 1, 2]
B = rng.normal(size=(n, n))
C0 = B @ B.T / n + 0.5 * np.eye(n)
Q0 = np.linalg.inv(C0)


def tilt(a, b):
    """dP_L/dP_0 proportional to exp(-U(x_F)), U = 1/2 x_F^T a x_F - b^T x_F: again Gaussian."""
    Q = Q0.copy()
    Q[np.ix_(F, F)] += a
    h = np.zeros(n)
    h[F] = b
    C = np.linalg.inv(Q)
    return C @ h, C


a1, a2 = np.diag([0.8, 0.2, 0.5]), np.array([[0.3, 0.1, 0.0], [0.1, 0.9, 0.2], [0.0, 0.2, 0.4]])
m1, C1 = tilt(a1, np.array([0.5, -0.2, 0.1]))
m2, C2 = tilt(a2, np.array([-0.3, 0.4, 0.2]))
full = bhattacharyya(m1, C1, m2, C2)
pocket = bhattacharyya(m1[F], C1[np.ix_(F, F)], m2[F], C2[np.ix_(F, F)])
out(f"  n = {n}, pocket F = {F}: D_B(full) = {full:.15f}, D_B(pocket marginals) = {pocket:.15f}")
check("full-ensemble divergence equals the pocket-marginal divergence", abs(full - pocket) < 1e-12)

# ---------------- [4] ----------------
out("\n[4] Theorem 3.3: Szego rate, constant c_0 and the values at N = 501")
fr = lambda r, t: (1 - r * r) / (1 - 2 * r * np.cos(t) + r * r)  # noqa: E731
g = lambda t: 0.5 * np.log((fr(.35, t) + fr(.55, t)) / 2) - 0.25 * np.log(fr(.35, t)) - 0.25 * np.log(fr(.55, t))  # noqa: E731
s = quad(g, -np.pi, np.pi, limit=400, epsabs=1e-15, epsrel=1e-13)[0] / (2 * np.pi)
out(f"  Szego rate by quadrature: s = {s:.12f}")
check("s = 0.0094968218 to ten digits", abs(s - 0.0094968218) < 5e-11)
consts = {}
for N in (50, 100, 200, 501, 1000):
    CA, CB = toeplitz(0.35 ** np.arange(N)), toeplitz(0.55 ** np.arange(N))
    Dc = 0.5 * logdet(0.5 * (CA + CB)) - 0.25 * (logdet(CA) + logdet(CB))
    consts[N] = Dc - s * N
    out(f"  N = {N:>4}: D_cov = {Dc:.10f}, D_cov - s N = {consts[N]:.10f}")
check("D_cov(N) - s N = c_0 = -0.0093715854 for N >= 50", all(abs(v + 0.0093715854) < 1e-9 for v in consts.values()))
N = 501
i = np.arange(1, N + 1, dtype=float)
dmu = 0.12 * i ** -0.75
CA, CB = toeplitz(0.35 ** np.arange(N)), toeplitz(0.55 ** np.arange(N))
Cb = 0.5 * (CA + CB)
Dm = float(dmu @ np.linalg.solve(Cb, dmu) / 8)
Dc = 0.5 * logdet(Cb) - 0.25 * (logdet(CA) + logdet(CB))
share = Dc / (Dm + Dc)
out(f"  N = 501: D_mean = {Dm:.7f}, D_cov = {Dc:.5f}, CKI = {Dm + Dc:.5f}, covariance share = {100 * share:.3f}%")
check("D_mean(501) = 0.0024408, D_cov(501) = 4.74854, CKI = 4.75098, share 99.949%",
      abs(Dm - 0.0024408) < 5e-8 and abs(Dc - 4.74854) < 5e-6 and abs(Dm + Dc - 4.75098) < 5e-6 and abs(100 * share - 99.949) < 5e-4)

out(f"\nOVERALL: {'PASS' if ok_all else 'FAIL'}")
with open(OUT, "w", encoding="utf-8") as f:
    f.write("\n".join(lines) + "\n")
sys.exit(0 if ok_all else 1)
