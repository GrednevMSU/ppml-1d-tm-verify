r"""PPML v3.0 — 1-D TM shared kernel (Python port).

Faithful port of the MATLAB `general/` helpers. Numerical semantics are preserved over
Pythonic style: MATLAB `\` (mldivide) and `/` (mrdivide) are reproduced exactly, the
sqrt branch convention is bit-for-bit, and the hardcoded vacuum impedance `376.730`
(ORIGINAL_CODE_FINDINGS F1) is replicated — do NOT substitute scipy.constants.

Original: PPML v3.0, S. Zanotto (BSD). See matlab_src/general/*.m.
"""
from __future__ import annotations

import numpy as np

# ORIGINAL_CODE_FINDINGS F1: original hardcodes the vacuum impedance truncated to
# 376.730 Ohm (physical value 376.730313668...). Replicated bug-for-bug; a fix is
# gated behind the `fix_vacuum_impedance` flag in tm_1d (default OFF).
Z0_TRUNCATED = 376.730
Z0_FULL = 376.73031366857  # only used when fix_vacuum_impedance=True (DEVIATION)


# Singular-matrix policy (NAME_MAP D2). At an EXACT singular point (e.g. a Rayleigh/Wood
# cutoff where a mode kz -> 0) MATLAB `\` warns and returns Inf/NaN-laden output; numpy
# `solve` raises LinAlgError. The two ORIGINAL engines already DISAGREE there (MATLAB
# returns NaN, Octave returns finite garbage), so no port choice can be bit-for-bit.
#   'raise'     (default) -> let LinAlgError propagate. Safer for the user; the outcome
#                            model records it as DEVIATION(flag=on_singular) at that
#                            exact-cutoff fixture (near_anomaly, non-gating).
#   'replicate'           -> best-effort min-norm lstsq (NOT bit-exact to MATLAB;
#                            documented). Lets the pipeline produce a value if desired.
_SOLVE_MODE = "raise"


def set_solve_mode(mode: str) -> str:
    """Set the singular-matrix policy ('raise'|'replicate'); returns the previous mode."""
    global _SOLVE_MODE
    if mode not in ("raise", "replicate"):
        raise ValueError("on_singular must be 'raise' or 'replicate'")
    prev, _SOLVE_MODE = _SOLVE_MODE, mode
    return prev


def mldivide(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """MATLAB ``A\\B`` — solve A X = B. Mirrors left-division of the original."""
    if _SOLVE_MODE == "replicate":
        try:
            return np.linalg.solve(A, B)
        except np.linalg.LinAlgError:
            return np.linalg.lstsq(A, B, rcond=None)[0]
    return np.linalg.solve(A, B)


def mrdivide(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """MATLAB ``A/B`` — solve X B = A, i.e. ``A @ inv(B)`` without forming the inverse."""
    if _SOLVE_MODE == "replicate":
        try:
            return np.linalg.solve(B.T, A.T).T
        except np.linalg.LinAlgError:
            return np.linalg.lstsq(B.T, A.T, rcond=None)[0].T
    return np.linalg.solve(B.T, A.T).T


def sqrt_whittaker(qq):
    """Branch-fixed complex square root (port of general/sqrt_whittaker.m).

    Chooses the root in the upper half-plane (Im > 0) so evanescent modes decay with
    the code's e^{-i w t} / e^{+i q z} convention; the real-axis tie is broken toward
    Re > 0. Exact 4-case sign logic of the original, vectorized elementwise.
    """
    q = np.sqrt(np.asarray(qq, dtype=complex))
    # MATLAB:
    #   if Re(q) >= 0:  if Im(q) <  0:  q = -q
    #   else (Re<0):    if Im(q) <= 0:  q = -q
    flip = np.where(np.real(q) >= 0.0, np.imag(q) < 0.0, np.imag(q) <= 0.0)
    return np.where(flip, -q, q)


def smpropag_fw_cond(S1, S2, p1, p2, A1, A2, f1, f2, sigma):
    """Forward S-matrix recursion across one conducting interface + layer.

    Port of general/smpropag_fw_cond.m. `f1`,`f2` are the layer phase VECTORS
    (exp(i q d) diagonals); `sigma` is the scalar interface term sigma*Z0/k0.
    """
    A1A2 = mldivide(A1, A2)
    p1A2 = mldivide(p1, A2)
    p1p2 = mldivide(p1, p2)

    I1 = ( A1A2 + sigma * p1A2 + p1p2) / 2.0
    I2 = (-A1A2 - sigma * p1A2 + p1p2) / 2.0
    I3 = (-A1A2 + sigma * p1A2 + p1p2) / 2.0
    I4 = ( A1A2 - sigma * p1A2 + p1p2) / 2.0

    Df1 = np.diag(f1)
    M = I1 - Df1 @ S2 @ I3
    S1r = mldivide(M, Df1 @ S1)
    S2r = mldivide(M, Df1 @ S2 @ I4 - I2) @ np.diag(f2)
    return S1r, S2r


def smpropag_bw_cond(S3, S4, p1, p2, A1, A2, f1, f2, sigma):
    """Backward S-matrix recursion across one conducting interface + layer.

    Port of general/smpropag_bw_cond.m.
    """
    A1A2 = mldivide(A1, A2)
    p1A2 = mldivide(p1, A2)
    p1p2 = mldivide(p1, p2)

    J1 = ( A1A2 - sigma * p1A2 + p1p2) / 2.0
    J2 = (-A1A2 + sigma * p1A2 + p1p2) / 2.0
    J3 = (-A1A2 - sigma * p1A2 + p1p2) / 2.0
    J4 = ( A1A2 + sigma * p1A2 + p1p2) / 2.0

    Df1 = np.diag(f1)
    M = J4 - Df1 @ S3 @ J2
    S3r = mldivide(M, Df1 @ S3 @ J1 - J3) @ np.diag(f2)
    S4r = mldivide(M, Df1 @ S4)
    return S3r, S4r
