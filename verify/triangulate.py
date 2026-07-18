"""triangulate.py — independent cross-check via the ZERO-CONTRAST LIMIT (non-gating INFO).

A grating whose inclusion permittivity equals its host permittivity is not a grating at
all — it is a uniform planar multilayer. In that limit the PPML RCWA result MUST reduce to
the ordinary thin-film transfer-matrix (TMM) result. We implement an INDEPENDENT TM
(p-polarization) Abeles transfer-matrix solver here (no PPML code reused) and compare it to
`RTA_1d_tm` evaluated with epsxA==epsxB, epszA==epszB (zero contrast).

Agreement of two independent formulations (RCWA reduced vs first-principles TMM) is a
strong trust signal. Disagreement would be a finding, NOT a tolerance adjustment.

Run: .venv/bin/python verify/triangulate.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "python_src"))
from ppml_1d_tm import RTA_1d_tm  # noqa: E402
from ppml_1d_tm.general import sqrt_whittaker  # noqa: E402


def tmm_tm(eps_layers, d_layers, epssup, epssub, k0, kpar):
    """Independent TM (p-pol) transfer-matrix R, T for a planar multilayer.

    eps_layers, d_layers: the L INTERNAL layers (super/sub handled separately).
    Uses the same branch convention as PPML (sqrt_whittaker, Im(kz) >= 0) so the two are
    compared on equal footing. Admittance for TM: y = kz / eps.
    """
    def kz(eps):
        return sqrt_whittaker(np.array([eps * k0**2 - kpar**2]))[0]

    kz0 = kz(epssup); y0 = kz0 / epssup
    kzs = kz(epssub); ys = kzs / epssub

    M = np.eye(2, dtype=complex)
    for eps, d in zip(eps_layers, d_layers):
        kzl = kz(eps); y = kzl / eps
        beta = kzl * d
        c, s = np.cos(beta), np.sin(beta)
        Ml = np.array([[c, 1j * s / y],
                       [1j * y * s, c]], dtype=complex)
        M = M @ Ml

    m11, m12, m21, m22 = M[0, 0], M[0, 1], M[1, 0], M[1, 1]
    denom = (y0 * m11 + y0 * ys * m12 + m21 + ys * m22)
    r = (y0 * m11 + y0 * ys * m12 - m21 - ys * m22) / denom
    t = 2 * y0 / denom
    R = abs(r) ** 2
    T = np.real(ys) / np.real(y0) * abs(t) ** 2   # power transmittance (real strate)
    return R, T


def _cases():
    return [
        # (name, a, eps_internal, d_internal, epssup, epssub, k0, theta_deg)
        ("single dielectric slab", 1.2, [6.25], [0.6], 1.0, 1.0, 2 * np.pi / 1.0, 20.0),
        ("two-layer stack", 1.5, [12.0, 2.25], [0.4, 0.7], 1.0, 1.0, 2 * np.pi / 1.1, 30.0),
        ("immersion, real sub", 1.0, [9.0], [0.5], 2.25, 4.0, 2 * np.pi / 1.0, 15.0),
        ("subwavelength", 0.4, [11.0], [0.3], 1.0, 1.0, 2 * np.pi / 2.0, 12.0),
    ]


def run():
    print("  ZERO-CONTRAST TRIANGULATION — PPML RCWA (zero contrast) vs independent TMM")
    print("  case                         R_ppml      R_tmm       T_ppml      T_tmm      max|Δ|")
    worst = 0.0
    for name, a, eps_int, d_int, esup, esub, k0, th in _cases():
        kpar = k0 * np.sqrt(esup) * np.sin(th * np.pi / 180)
        L = len(eps_int)
        d = [1.0] + list(d_int) + [1.0]
        # zero contrast: A == B, so f is irrelevant; grating reduces to uniform layers
        RR, TT, _ = RTA_1d_tm(a, L, esup, esub, list(eps_int), list(eps_int),
                              list(eps_int), list(eps_int), [0] * (L + 1),
                              [0.5] * L, d, 8, k0, kpar)
        R_t, T_t = tmm_tm(eps_int, d_int, esup, esub, k0, kpar)
        dmax = max(abs(RR.real - R_t), abs(TT.real - T_t))
        worst = max(worst, dmax)
        print(f"  {name:26s}  {RR.real:.8f}  {R_t:.8f}  {TT.real:.8f}  {T_t:.8f}  {dmax:.2e}")
    print(f"  worst |Δ| across cases: {worst:.2e}  "
          f"({'AGREE (independent formulations match)' if worst < 1e-9 else 'FINDING — investigate'})")
    return worst


if __name__ == "__main__":
    import sys as _sys
    _sys.exit(0 if run() < 1e-9 else 1)
