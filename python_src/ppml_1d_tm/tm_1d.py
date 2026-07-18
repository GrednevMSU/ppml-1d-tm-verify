"""PPML v3.0 — 1-D TM public functions (Python port).

Ports of matlab_src/1d_tm/{SM_1d_tm,RTA_1d_tm,field_1d_tm}.m. Numerical semantics
preserved bug-for-bug (ORIGINAL_CODE_FINDINGS F1/F2/F3). Original comments in the .m
files remain authoritative; line references below point into them.

Conventions (VERIFICATION_NOTES §6): time e^{-i w t}; loss Im(eps) > 0; sqrt branch
Im(q) >= 0 (evanescent decay). Li inverse rule for TM at the epsx build (SM_1d_tm.m:115).

Gating fixes (both default OFF -> replicate original; ON -> DEVIATION in reports):
  fix_vacuum_impedance : use full-precision Z0 instead of 376.730 (F1).
  guard_real_epssup    : raise on complex epssup instead of silently returning
                         nonphysical values (F2). NB the original RTA guard at
                         RTA.m:183 is a telescoping tautology (F3) and never fires.
"""
from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass

import numpy as np

from .general import (
    Z0_FULL,
    Z0_TRUNCATED,
    mldivide,
    mrdivide,
    set_solve_mode,
    smpropag_bw_cond,
    smpropag_fw_cond,
    sqrt_whittaker,
)


@contextmanager
def _solve_policy(on_singular: str):
    """Scope the module singular-matrix policy for one public call (NAME_MAP D2)."""
    prev = set_solve_mode(on_singular)
    try:
        yield
    finally:
        set_solve_mode(prev)


class PPMLEnergyError(RuntimeError):
    """Port of the original's ``error('RTA:EnNotCons', ...)`` (RTA_1d_tm.m:183).

    Carries the MATLAB identifier so the outcome model (Phase 5) can compare error
    ids across engines. NOTE (F3): the checked quantity is identically 1 by
    telescoping, so this fires only on NaN/Inf — it is NOT a real energy check.
    """

    identifier = "RTA:EnNotCons"


def _as_vec(x) -> np.ndarray:
    """Coerce a scalar / 0-d / list to a 1-D complex array (handles L=0 empties)."""
    a = np.atleast_1d(np.asarray(x, dtype=complex))
    return a


@dataclass
class _Stack:
    """Eigen-solutions + assembled S-matrices for one structure (shared core)."""

    n: np.ndarray
    npw: int
    kx: np.ndarray
    q: list          # each (npw,npw) diagonal
    phi: list
    A: list
    epsx: list
    etaz: list
    S1: list
    S2: list
    S3: list
    S4: list


def _build_stack(a, L, epssup, epssub, epsxA, epszA, epsxB, epszB, sigma, f, d,
                 halfnpw, k0, kpar, z0) -> _Stack:
    """Shared eigen-decomposition + Redheffer forward/backward recursion.

    Mirrors lines 76-148 of every 1-D TM function (identical in RTA/SM/field). Factored
    to a single helper so the three ports cannot drift — a maintainability refactor with
    NO change to the arithmetic.
    """
    n = np.arange(-halfnpw, halfnpw + 1)
    npw = n.size
    kx = (2.0 * np.pi / a) * n + kpar

    epsxA = _as_vec(epsxA); epszA = _as_vec(epszA)
    epsxB = _as_vec(epsxB); epszB = _as_vec(epszB)
    sigma = _as_vec(sigma)
    f = _as_vec(f)
    d = np.asarray(d, dtype=float)

    I = np.eye(npw)
    q = [None] * (L + 2)
    phi = [None] * (L + 2)
    A = [None] * (L + 2)
    epsx = [None] * (L + 2)
    etaz = [None] * (L + 2)

    # --- superstrate (index 0) and substrate (index L+1) --------------------
    q[0] = np.diag(sqrt_whittaker(epssup * k0**2 - kx * kx))
    phi[0] = I.copy()
    A[0] = mrdivide(np.diag(k0**2 - kx * kx / epssup), q[0])
    epsx[0] = I / epssup
    etaz[0] = I * epssup

    q[L + 1] = np.diag(sqrt_whittaker(epssub * k0**2 - kx * kx))
    phi[L + 1] = I.copy()
    A[L + 1] = mrdivide(np.diag(k0**2 - kx * kx / epssub), q[L + 1])
    epsx[L + 1] = I / epssub
    etaz[L + 1] = I * epssub

    # --- internal layers 1..L ----------------------------------------------
    for lay in range(L):                       # MATLAB l = 1..L -> layer index lay+1
        li = lay + 1
        dn = n[:, None] - n[None, :]
        with np.errstate(divide="ignore", invalid="ignore"):
            F = np.sin(np.pi * f[lay] * dn) / (np.pi * dn)
        np.fill_diagonal(F, f[lay])            # removable sinc singularity -> f(l)

        epsx[li] = (1.0 / epsxB[lay] - 1.0 / epsxA[lay]) * F + (1.0 / epsxA[lay]) * I
        etaz[li] = (epszB[lay] - epszA[lay]) * F + epszA[lay] * I

        Kz = k0**2 * I - np.diag(kx) @ mldivide(etaz[li], np.diag(kx))
        w, V = np.linalg.eig(mldivide(epsx[li], Kz))
        q[li] = np.diag(sqrt_whittaker(w))
        phi[li] = V
        A[li] = mrdivide(Kz @ V, q[li])

    # layer phase vectors exp(i q d)
    expqd = [np.exp(1j * np.diag(q[l]) * d[l]) for l in range(L + 2)]

    # --- forward recursion --------------------------------------------------
    Z = np.zeros((npw, npw), dtype=complex)
    S1 = [Z.copy() for _ in range(L + 2)]
    S2 = [Z.copy() for _ in range(L + 2)]
    S1[0] = I.copy()
    for l in range(L + 1):
        S1[l + 1], S2[l + 1] = smpropag_fw_cond(
            S1[l], S2[l], phi[l], phi[l + 1], A[l], A[l + 1],
            expqd[l], expqd[l + 1], sigma[l] * z0 / k0)

    # --- backward recursion -------------------------------------------------
    S3 = [Z.copy() for _ in range(L + 2)]
    S4 = [Z.copy() for _ in range(L + 2)]
    S4[L + 1] = I.copy()
    for p in range(L + 1, 0, -1):
        S3[p - 1], S4[p - 1] = smpropag_bw_cond(
            S3[p], S4[p], phi[p], phi[p - 1], A[p], A[p - 1],
            expqd[p], expqd[p - 1], sigma[p - 1] * z0 / k0)

    return _Stack(n, npw, kx, q, phi, A, epsx, etaz, S1, S2, S3, S4)


def _z0(fix_vacuum_impedance: bool) -> float:
    return Z0_FULL if fix_vacuum_impedance else Z0_TRUNCATED


def _check_epssup(epssup, guard_real_epssup: bool):
    if guard_real_epssup and np.imag(epssup) != 0.0:
        raise ValueError("epssup must be real (guard_real_epssup=ON, DEVIATION)")


def SM_1d_tm(a, L, epssup, epssub, epsxA, epszA, epsxB, epszB, sigma, f, d,
             halfnpw, k0, kpar, *, fix_vacuum_impedance=False, guard_real_epssup=False,
             on_singular="raise"):
    """S-matrix 0-order coefficients (rl, rr, tlr, trl). Port of SM_1d_tm.m."""
    _check_epssup(epssup, guard_real_epssup)
    with _solve_policy(on_singular):
        st = _build_stack(a, L, epssup, epssub, epsxA, epszA, epsxB, epszB, sigma, f, d,
                          halfnpw, k0, kpar, _z0(fix_vacuum_impedance))
        h = halfnpw
        q0 = st.q[0][h, h]
        qN = st.q[L + 1][h, h]
        d = np.asarray(d, dtype=float)
        rl = st.S3[0][h, h] * np.exp(-1j * q0 * d[0])
        rr = st.S2[L + 1][h, h] * np.exp(-1j * qN * d[L + 1])
        tlr = st.S1[L + 1][h, h] * np.exp(-1j * q0 * d[0])
        trl = st.S4[0][h, h] * np.exp(-1j * qN * d[L + 1])
    return rl, rr, tlr, trl


def RTA_1d_tm(a, L, epssup, epssub, epsxA, epszA, epsxB, epszB, sigma, f, d,
              halfnpw, k0, kpar, *, fix_vacuum_impedance=False, guard_real_epssup=False,
              on_singular="raise"):
    """Reflectance, Transmittance, per-layer Absorbance. Port of RTA_1d_tm.m.

    Returns (RR, TT, AA) with AA a length-L array. Raises PPMLEnergyError only if the
    (tautological, F3) guard's residual exceeds 1e-5 (i.e. on NaN/Inf).
    """
    _check_epssup(epssup, guard_real_epssup)
    with _solve_policy(on_singular):
        st = _build_stack(a, L, epssup, epssub, epsxA, epszA, epsxB, epszB, sigma, f, d,
                          halfnpw, k0, kpar, _z0(fix_vacuum_impedance))
        npw, h = st.npw, halfnpw
        d = np.asarray(d, dtype=float)
        I = np.eye(npw)

        a_s = np.zeros((npw, L + 2), dtype=complex)
        b_s = np.zeros((npw, L + 2), dtype=complex)
        a_s[h, 0] = np.exp(-1j * st.q[0][h, h] * d[0]) * np.sqrt(epssup) / k0
        incflux = st.q[0][h, h] / k0**2

        for l in range(L + 2):
            a_s[:, l] = mldivide(I - st.S2[l] @ st.S3[l], st.S1[l] @ a_s[:, 0])
            b_s[:, l] = mldivide(I - st.S3[l] @ st.S2[l], st.S3[l] @ st.S1[l] @ a_s[:, 0])

        # dtype=complex: MATLAB `flux=zeros(...)` is real but AUTO-PROMOTES because it is
        # divided by incflux (complex when epssup is complex, F2). numpy does NOT promote
        # — a real array would silently drop the imaginary part. Bug-for-bug fidelity
        # requires complex (caught by out_of_domain_epssup_complex_RTA: real flux gave
        # real RR 0.4018, MATLAB gives 0.4018+0.0938i).
        flux = np.zeros(L + 1, dtype=complex)
        for l in range(1, L + 2):                 # MATLAB l = 2..L+2
            expd = np.exp(1j * np.diag(st.q[l]) * d[l])
            hy = st.phi[l] @ (a_s[:, l] + np.diag(expd) @ b_s[:, l])
            ex = st.A[l] @ (a_s[:, l] - np.diag(expd) @ b_s[:, l])
            flux[l - 1] = np.real(hy @ np.conj(ex)) / incflux

        RR = 1.0 - flux[0]
        TT = flux[L]
        AA = np.array([flux[l] - flux[l + 1] for l in range(L)])

    if abs(np.sum(AA) + RR + TT - 1.0) > 1e-5:   # F3: tautology, fires only on NaN/Inf
        raise PPMLEnergyError("Energy not conserved")
    return RR, TT, AA


def field_1d_tm(a, L, epssup, epssub, epsxA, epszA, epsxB, epszB, sigma, f, d,
                halfnpw, k0, kpar, nx, nz, *,
                fix_vacuum_impedance=False, guard_real_epssup=False, on_singular="raise"):
    """Real-space fields Ex, Ez and unit-cell-averaged Poynting Sz. Port of field_1d_tm.m.

    Returns (x, z, Ex, Ez, Sz). Ex/Ez are (nx, sum(nz)); Sz is (sum(nz),).
    """
    _check_epssup(epssup, guard_real_epssup)
    with _solve_policy(on_singular):
        st = _build_stack(a, L, epssup, epssub, epsxA, epszA, epsxB, epszB, sigma, f, d,
                          halfnpw, k0, kpar, _z0(fix_vacuum_impedance))
        npw, h = st.npw, halfnpw
        d = np.asarray(d, dtype=float)
        nz = np.atleast_1d(np.asarray(nz, dtype=int))
        I = np.eye(npw)
        kx = st.kx

        a_s = np.zeros((npw, L + 2), dtype=complex)
        b_s = np.zeros((npw, L + 2), dtype=complex)
        a_s[h, 0] = np.exp(-1j * st.q[0][h, h] * d[0]) * np.sqrt(epssup) / k0
        incflux = st.q[0][h, h] / k0**2

        for l in range(L + 2):
            a_s[:, l] = mldivide(I - st.S2[l] @ st.S3[l], st.S1[l] @ a_s[:, 0])
            b_s[:, l] = mldivide(I - st.S3[l] @ st.S2[l], st.S3[l] @ st.S1[l] @ a_s[:, 0])

        ntot = int(np.sum(nz))
        Ez = np.zeros((nx, ntot), dtype=complex)
        Ex = np.zeros((nx, ntot), dtype=complex)
        Sz = np.zeros(ntot, dtype=complex)  # complex: MATLAB auto-promotes via /incflux (F2)
        x = np.linspace(-a / 2.0, a / 2.0, nx)

        z0 = 0.0
        nz0 = 0
        z = np.zeros(0)
        qdiag = [np.diag(st.q[l]) for l in range(L + 2)]
        for l in range(L + 2):
            zz = np.linspace(0.0, d[l], nz[l])
            z = np.concatenate([z, zz + z0])
            for iz in range(nz[l]):
                e_a = np.exp(1j * qdiag[l] * zz[iz])
                e_b = np.exp(1j * qdiag[l] * (d[l] - zz[iz]))
                hy = st.phi[l] @ (np.diag(e_a) @ a_s[:, l] + np.diag(e_b) @ b_s[:, l])
                ex = st.A[l] @ (np.diag(e_a) @ a_s[:, l] - np.diag(e_b) @ b_s[:, l])
                ez = -mldivide(st.etaz[l], np.diag(kx)) @ hy
                for ix in range(nx):
                    phase = np.exp(1j * kx * x[ix])
                    Ez[ix, iz + nz0] = phase @ ez
                    Ex[ix, iz + nz0] = phase @ ex
                Sz[iz + nz0] = np.real(np.conj(hy) @ ex) / incflux
            z0 = z0 + d[l]
            nz0 = nz0 + int(nz[l])

    return x, z, Ex, Ez, Sz
