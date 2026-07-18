"""convergence.py — truncation-order (N) convergence study for the PPML 1-D TM port.

Non-gating standalone artifact (domain addendum, Phase 5). Sweeps the number of Fourier
harmonics N = 2*halfnpw+1 for a fixed grating, records the physical outputs R/T vs N,
and estimates the converged limit by Richardson/Aitken extrapolation with a
discretization-error estimate |value(N) - limit|. Runs the PORT only (no MATLAB needed).

This doubles as evidence that the truncation order is a first-class fixture axis: the
MATLAB-vs-Python equivalence is ALWAYS gated at identical N; this study characterises the
physics of N-convergence itself, which is common to both implementations.

Run: .venv/bin/python verify/convergence.py
"""
from __future__ import annotations

import datetime as _dt
import os
import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "python_src"))
from ppml_1d_tm import RTA_1d_tm  # noqa: E402


def _structure():
    """A lossless dielectric grating (high/low index) — TM convergence in N is clean and
    monotone here, so Richardson/Aitken extrapolation is well-behaved. (Metallic TM
    gratings converge far more slowly and oscillate — a known RCWA fact, not a port
    issue; see VERIFICATION_NOTES.)"""
    DA, DB = 12.0, 2.25            # Si-like / SiO2-like, real
    return dict(
        a=1.2, L=1,
        f=[0.5],
        d=[1.0, 0.5, 1.0],
        epsxA=[DA], epszA=[DA], epsxB=[DB], epszB=[DB],
        sigma=[0, 0], epssup=1.0, epssub=1.0,
        k0=2 * np.pi / 1.0, theta=15.0)


def aitken(seq):
    """Aitken Δ² extrapolation of the tail of a scalar sequence."""
    x = np.asarray(seq, float)
    if x.size < 3:
        return x[-1]
    a, b, c = x[-3], x[-2], x[-1]
    denom = (c - 2 * b + a)
    if abs(denom) < 1e-300:
        return c
    return a - (b - a) ** 2 / denom


def run():
    S = _structure()
    kpar = S["k0"] * np.sin(S["theta"] * np.pi / 180)
    ladder = [2, 3, 4, 6, 8, 10, 12, 15, 18, 20, 25, 30]
    Ns, Rs, Ts = [], [], []
    for hp in ladder:
        RR, TT, _ = RTA_1d_tm(S["a"], S["L"], S["epssup"], S["epssub"],
                              S["epsxA"], S["epszA"], S["epsxB"], S["epszB"],
                              S["sigma"], S["f"], S["d"], hp, S["k0"], kpar)
        Ns.append(2 * hp + 1); Rs.append(RR.real); Ts.append(TT.real)

    R_lim = aitken(Rs); T_lim = aitken(Ts)
    R_err = [abs(r - R_lim) for r in Rs]
    T_err = [abs(t - T_lim) for t in Ts]

    print("  N      R              |R-R∞|         T              |T-T∞|")
    for N, r, re, t, te in zip(Ns, Rs, R_err, Ts, T_err):
        print(f"  {N:3d}   {r:.10f}   {re:.2e}   {t:.10f}   {te:.2e}")
    print(f"  Richardson/Aitken limit: R∞={R_lim:.10f}  T∞={T_lim:.10f}")
    print(f"  discretization error at N={Ns[-1]}: R {R_err[-1]:.2e}, T {T_err[-1]:.2e}")

    figures_dir = str(ROOT / "reports" / "figures")
    os.makedirs(figures_dir, exist_ok=True)
    date = _dt.date.today().isoformat()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    ax1.plot(Ns, Rs, "o-", label="R(N)")
    ax1.plot(Ns, Ts, "s-", label="T(N)")
    ax1.axhline(R_lim, color="C0", ls=":", lw=1)
    ax1.axhline(T_lim, color="C1", ls=":", lw=1)
    ax1.set_xlabel("N = 2·halfnpw + 1"); ax1.set_ylabel("R, T")
    ax1.set_title("Physical outputs vs truncation order N")
    ax1.legend(); ax1.grid(True, alpha=0.25)

    ax2.loglog(Ns, np.maximum(R_err, 1e-18), "o-", label="|R(N) − R∞|")
    ax2.loglog(Ns, np.maximum(T_err, 1e-18), "s-", label="|T(N) − T∞|")
    ax2.set_xlabel("N"); ax2.set_ylabel("discretization error")
    ax2.set_title("Convergence to the Richardson/Aitken limit")
    ax2.legend(); ax2.grid(True, which="both", alpha=0.25)
    fig.text(0.99, 0.01, f"PPML 1-D TM convergence study · {date}",
             ha="right", va="bottom", fontsize=7, color="#888")
    out = os.path.join(figures_dir, "convergence.png")
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out}")
    return R_lim, T_lim, R_err[-1], T_err[-1]


if __name__ == "__main__":
    run()
