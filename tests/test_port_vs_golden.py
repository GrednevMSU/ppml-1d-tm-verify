"""Unit tests: Python port vs the committed MATLAB golden reference_outputs/.

No MATLAB/Octave needed — runs against reference_outputs/matlab/*.mat that are
committed to the repo. Mirrors the gate logic of verify/compare.py (OR pass logic,
budget-derived tolerances) but as standalone pytest cases, one per fixture.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "python_src"))
from ppml_1d_tm import RTA_1d_tm, SM_1d_tm, field_1d_tm  # noqa: E402

OUTPUTS = {"RTA": ["RR", "TT", "AA"], "SM": ["rl", "rr", "tlr", "trl"],
           "field": ["Ex", "Ez", "Sz"]}
# budget-derived tolerances (see verify/tolerances.yaml)
TOL = {"RTA": (1e-11, 1e-10), "SM": (1e-11, 1e-10), "field": (1e-9, 1e-11)}
NEAR_ANOMALY = (1e-3, 1e-4)   # loose INFO tier

_MANIFEST = json.loads((ROOT / "fixtures" / "manifest.json").read_text())
_GOLDEN = ROOT / "reference_outputs" / "matlab"


def _phys_fixtures():
    out = []
    for f in _MANIFEST["fixtures"]:
        if f["expect"] != "ok":
            continue
        if not (_GOLDEN / f"{f['id']}.mat").exists():
            continue
        out.append(f)
    return out


PHYS = _phys_fixtures()


def _vec(x):
    return np.atleast_1d(np.asarray(x, dtype=complex)).ravel()


def _run(fx):
    args = (fx.a, int(fx.L), complex(fx.epssup), complex(fx.epssub),
            fx.epsxA, fx.epszA, fx.epsxB, fx.epszB, fx.sigma, fx.f, fx.d,
            int(fx.halfnpw), float(fx.k0), float(fx.kpar))
    if fx.func == "RTA":
        RR, TT, AA = RTA_1d_tm(*args)
        return {"RR": RR, "TT": TT, "AA": AA}
    if fx.func == "SM":
        rl, rr, tlr, trl = SM_1d_tm(*args)
        return {"rl": rl, "rr": rr, "tlr": tlr, "trl": trl}
    x, z, Ex, Ez, Sz = field_1d_tm(*args, int(fx.nx), fx.nz)
    return {"Ex": Ex, "Ez": Ez, "Sz": Sz}


def test_corpus_nonempty():
    assert len(PHYS) > 200, "expected the full frozen corpus of physical fixtures"


def test_manifest_matches_committed_fixtures():
    assert _MANIFEST["corpus_version"] >= 3
    assert _MANIFEST["n_fixtures"] == len(_MANIFEST["fixtures"])


@pytest.mark.parametrize("f", PHYS, ids=[f["id"] for f in PHYS])
def test_port_matches_golden(f):
    fid, func = f["id"], f["function"]
    tol_rel, tol_abs = NEAR_ANOMALY if f["near_anomaly"] else TOL[func]
    ref = loadmat(_GOLDEN / f"{fid}.mat", struct_as_record=False, squeeze_me=True)["ref"]
    if str(ref.outcome) != "ok":
        pytest.skip(f"{fid}: golden outcome {ref.outcome}")
    fx = loadmat(ROOT / "fixtures" / f"{fid}.mat",
                 struct_as_record=False, squeeze_me=True)["fx"]
    try:
        out = _run(fx)
    except np.linalg.LinAlgError:
        # D2: at an EXACT Wood cutoff the system is singular; the port raises by design
        # (on_singular='raise') where the reference engines themselves disagree
        # (MATLAB NaN vs Octave finite). Documented DEVIATION(flag=on_singular), non-gating.
        if f["near_anomaly"]:
            pytest.xfail(f"{fid}: singular at exact cutoff — DEVIATION(flag=on_singular), see D2")
        raise
    for k in OUTPUTS[func]:
        pv = _vec(out[k]); rv = _vec(getattr(ref, k))
        assert pv.size == rv.size, f"{fid}.{k}: shape {pv.size} != {rv.size}"
        if pv.size == 0:
            continue
        aerr = np.abs(pv - rv)
        scale = np.abs(rv)
        with np.errstate(divide="ignore", invalid="ignore"):
            rerr = np.where(scale > 0, aerr / scale, 0.0)
        el_pass = (aerr <= tol_abs) | (rerr <= tol_rel)
        assert np.all(el_pass), (
            f"{fid}.{k}: worst abs {aerr.max():.2e} (tol {tol_abs:g}), "
            f"worst rel {rerr[scale > tol_abs].max() if np.any(scale > tol_abs) else 0:.2e} "
            f"(tol {tol_rel:g})")


def test_vacuum_impedance_replicated_bug_for_bug():
    """F1: the port must hardcode 376.730, not the full-precision Z0."""
    from ppml_1d_tm.general import Z0_TRUNCATED
    assert Z0_TRUNCATED == 376.730


def test_sqrt_whittaker_upper_half_plane():
    """Branch convention: root chosen with Im >= 0 (evanescent decay)."""
    from ppml_1d_tm import sqrt_whittaker
    q = sqrt_whittaker(np.array([-4.0, 1.0, -1j, 3 + 4j]))
    # negative real -> +2i (Im>0); positive real -> real; etc.
    assert np.all(np.imag(q) >= -1e-15)


def test_complex_epssup_is_nonphysical_not_raising():
    """F2/F3: complex epssup returns (nonphysical) values; RTA does NOT raise."""
    f = next(x for x in _MANIFEST["fixtures"]
             if x["id"] == "out_of_domain_epssup_complex_RTA")
    fx = loadmat(ROOT / "fixtures" / f"{f['id']}.mat",
                 struct_as_record=False, squeeze_me=True)["fx"]
    RR, TT, AA = RTA_1d_tm(fx.a, int(fx.L), complex(fx.epssup), complex(fx.epssub),
                           fx.epsxA, fx.epszA, fx.epsxB, fx.epszB, fx.sigma, fx.f, fx.d,
                           int(fx.halfnpw), float(fx.k0), float(fx.kpar))
    assert abs(RR.imag) > 1e-3   # nonphysical: reflectance has a nonzero imaginary part
