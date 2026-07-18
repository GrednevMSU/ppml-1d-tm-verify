"""benchmark.py — honest wall-clock timing of the Python port.

Speed is NOT a goal of this project (correctness/equivalence is). This benchmark exists so
any speed claim is measured, not asserted. It times RTA_1d_tm over a representative
structure at a ladder of truncation orders N, on this machine. Pair with bench_matlab.m
(run in your MATLAB) to fill the ratio column in BENCHMARK.md.

Run: .venv/bin/python bench/benchmark.py
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "python_src"))
from ppml_1d_tm import RTA_1d_tm  # noqa: E402


def _args(halfnpw):
    epsAu = -4000 + 300j
    epsbg = 10.05
    k0 = 150 * (2 * np.pi / 1240)
    return (3.2, 3, 1.0, 1.0,
            [epsbg, epsbg, epsbg], [epsbg, epsbg, epsbg],
            [epsbg, epsAu, epsbg], [epsbg, epsAu, epsbg],
            [0, 0, 0, 0], [0.0, 0.73, 0.0], [0.0, 0.5, 0.05, 1.3, 0.0],
            halfnpw, k0, k0 * np.sin(0.1 * np.pi / 180))


def run():
    print("  N     halfnpw   calls   mean [ms]   median [ms]")
    rows = []
    for hp in [5, 10, 20, 30]:
        args = _args(hp)
        RTA_1d_tm(*args)  # warm up
        reps = 200 if hp <= 20 else 80
        ts = []
        for _ in range(reps):
            t0 = time.perf_counter()
            RTA_1d_tm(*args)
            ts.append((time.perf_counter() - t0) * 1e3)
        ts = np.array(ts)
        rows.append((2 * hp + 1, hp, reps, ts.mean(), np.median(ts)))
        print(f"  {2*hp+1:3d}    {hp:3d}      {reps:4d}    {ts.mean():8.3f}    {np.median(ts):8.3f}")
    return rows


if __name__ == "__main__":
    import platform
    print(f"  platform: {platform.platform()}  python {platform.python_version()} numpy {np.__version__}")
    run()
