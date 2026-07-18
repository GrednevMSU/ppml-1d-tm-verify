"""Non-gating INFO checks as pytest cases: zero-contrast triangulation + N-convergence."""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "verify"))


def test_zero_contrast_matches_independent_tmm():
    """A zero-contrast grating must reduce to the first-principles thin-film TMM."""
    from triangulate import run
    worst = run()
    assert worst < 1e-9, f"RCWA zero-contrast vs TMM disagree by {worst:.2e}"


def test_n_convergence_reaches_limit():
    """R/T converge to the Richardson/Aitken limit as N grows (dielectric grating)."""
    from convergence import run
    R_lim, T_lim, R_err, T_err = run()
    assert R_err < 1e-3 and T_err < 1e-3, "did not converge to ~1e-3 by N=61"
    assert abs((R_lim + T_lim) - 1.0) < 1e-6, "lossless: R+T must be 1 at the limit"
