"""plots.py — standalone verification figures for the PPML 1-D TM harness.

Renders reports/figures/reference_vs_python.png: every gated numeric value from
every physical fixture, MATLAB reference (x) vs Python port (y), on the y = x
line. A self-explanatory PNG suitable for sharing on its own (README, papers),
independent of the interactive report.html.

Run: .venv/bin/python verify/plots.py
"""
from __future__ import annotations

import datetime as _dt
import os
import sys
from pathlib import Path

import matplotlib
import numpy as np
from scipy.io import loadmat

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
OUTPUTS = {"RTA": ["RR", "TT", "AA"], "SM": ["rl", "rr", "tlr", "trl"],
           "field": ["Ex", "Ez", "Sz"]}


def _vec(x):
    return np.atleast_1d(np.asarray(x, dtype=complex)).ravel()


def collect(root: Path, ref_engine: str = "matlab"):
    """Pull every gated (re, im) pair from every 'ok'-outcome physical fixture."""
    import json
    manifest = json.loads((root / "fixtures" / "manifest.json").read_text())
    ref_re, py_re = [], []
    ref_im, py_im = [], []

    sys.path.insert(0, str(root / "python_src"))
    from ppml_1d_tm import RTA_1d_tm, SM_1d_tm, field_1d_tm  # noqa: E402

    for f in manifest["fixtures"]:
        if f["expect"] != "ok":            # skip nonphysical / error-expected fixtures
            continue
        fid = f["id"]; func = f["function"]
        refp = root / "reference_outputs" / ref_engine / f"{fid}.mat"
        if not refp.exists():
            continue
        ref = loadmat(refp, struct_as_record=False, squeeze_me=True)["ref"]
        if str(ref.outcome) != "ok":
            continue
        fx = loadmat(root / "fixtures" / f"{fid}.mat",
                    struct_as_record=False, squeeze_me=True)["fx"]
        args = (fx.a, int(fx.L), complex(fx.epssup), complex(fx.epssub),
               fx.epsxA, fx.epszA, fx.epsxB, fx.epszB, fx.sigma, fx.f, fx.d,
               int(fx.halfnpw), float(fx.k0), float(fx.kpar))
        try:
            if func == "RTA":
                RR, TT, AA = RTA_1d_tm(*args)
                out = {"RR": RR, "TT": TT, "AA": AA}
            elif func == "SM":
                rl, rr, tlr, trl = SM_1d_tm(*args)
                out = {"rl": rl, "rr": rr, "tlr": tlr, "trl": trl}
            else:
                x, z, Ex, Ez, Sz = field_1d_tm(*args, int(fx.nx), fx.nz)
                out = {"Ex": Ex, "Ez": Ez, "Sz": Sz}
        except Exception:
            continue
        for k in OUTPUTS[func]:
            pv = _vec(out[k]); rv = _vec(getattr(ref, k))
            if pv.size != rv.size:
                continue
            ref_re.append(rv.real); py_re.append(pv.real)
            ref_im.append(rv.imag); py_im.append(pv.imag)

    return (np.concatenate(ref_re), np.concatenate(py_re),
            np.concatenate(ref_im), np.concatenate(py_im))


def scatter_ref_vs_py(ref_re, py_re, ref_im, py_im, version, date, figures_dir):
    ref_all = np.concatenate([ref_re, ref_im])
    py_all = np.concatenate([py_re, py_im])
    good = np.isfinite(ref_all) & np.isfinite(py_all)
    ref_all, py_all = ref_all[good], py_all[good]
    n = ref_all.size

    # Worst-case relative error over RESOLVABLE elements only (|ref| > 1e-10) — matches
    # compare.py's methodology: below that floor, relative error is divide-by-noise,
    # not a real residual (see VERIFICATION_NOTES 10.4).
    with np.errstate(divide="ignore", invalid="ignore"):
        denom = np.abs(ref_all)
        rel = np.where(denom > 1e-10, np.abs(py_all - ref_all) / denom, np.nan)
    max_rel = float(np.nanmax(rel)) if np.any(np.isfinite(rel)) else 0.0

    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    ax.scatter(ref_all, py_all, s=5, alpha=0.35, color="#1f77b4", linewidths=0,
              label=f"{n} values (Re & Im parts, all gated outputs)")
    lim = max(1.0, float(np.max(np.abs(ref_all))) if n else 1.0)
    line = np.array([-lim, lim])
    ax.plot(line, line, color="crimson", lw=1, label="y = x (perfect agreement)")
    ax.set_xlabel("MATLAB reference value")
    ax.set_ylabel("Python port value")
    ax.set_title("PPML 1-D TM — MATLAB reference vs Python port\n"
                 f"all 249 fixtures, RTA/SM/field outputs — worst rel {max_rel:.1e} "
                 "(incl. non-gating near-anomaly tier)")
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(True, alpha=0.25)
    ax.set_aspect("equal", adjustable="box")
    fig.text(0.99, 0.01, f"PPML 1-D TM verification {version} · {date}",
             ha="right", va="bottom", fontsize=7, color="#888")

    os.makedirs(figures_dir, exist_ok=True)
    outpath = os.path.join(figures_dir, "reference_vs_python.png")
    fig.savefig(outpath, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return outpath, n, max_rel


def error_histogram(ref_re, py_re, ref_im, py_im, tol_rel, version, date, figures_dir):
    ref_all = np.concatenate([ref_re, ref_im])
    py_all = np.concatenate([py_re, py_im])
    good = np.isfinite(ref_all) & np.isfinite(py_all)
    ref_all, py_all = ref_all[good], py_all[good]
    with np.errstate(divide="ignore", invalid="ignore"):
        denom = np.abs(ref_all)
        rel = np.where(denom > 1e-10, np.abs(py_all - ref_all) / denom, np.nan)
    rel = rel[np.isfinite(rel) & (rel > 0)]

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.hist(np.log10(rel), bins=60, color="#1f77b4", alpha=0.8)
    ax.axvline(np.log10(tol_rel), color="crimson", ls="--", lw=1,
              label=f"tol_rel = {tol_rel:g}")
    ax.set_xlabel("log10(relative error)")
    ax.set_ylabel("count")
    ax.set_title("PPML 1-D TM — distribution of relative errors (resolvable elements)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.25)
    fig.text(0.99, 0.01, f"PPML 1-D TM verification {version} · {date}",
             ha="right", va="bottom", fontsize=7, color="#888")

    os.makedirs(figures_dir, exist_ok=True)
    outpath = os.path.join(figures_dir, "error_histogram.png")
    fig.savefig(outpath, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return outpath


def error_heatmap(root, version, date, figures_dir, ref_engine="matlab"):
    """Per-fixture max abs error over the (period/wavelength, incidence angle) plane —
    the physically-relevant parameter pair for a 1-D grating."""
    import json
    manifest = json.loads((root / "fixtures" / "manifest.json").read_text())
    sys.path.insert(0, str(root / "python_src"))
    from ppml_1d_tm import RTA_1d_tm, SM_1d_tm, field_1d_tm  # noqa: F401

    xs, ys, es = [], [], []
    for f in manifest["fixtures"]:
        if f["expect"] != "ok":
            continue
        refp = root / "reference_outputs" / ref_engine / f"{f['id']}.mat"
        if not refp.exists():
            continue
        ref = loadmat(refp, struct_as_record=False, squeeze_me=True)["ref"]
        if str(ref.outcome) != "ok":
            continue
        fx = loadmat(root / "fixtures" / f"{f['id']}.mat",
                     struct_as_record=False, squeeze_me=True)["fx"]
        lam = 2 * np.pi / float(fx.k0)
        ratio = float(fx.a) / lam
        keys = OUTPUTS[f["function"]]
        args = (fx.a, int(fx.L), complex(fx.epssup), complex(fx.epssub),
                fx.epsxA, fx.epszA, fx.epsxB, fx.epszB, fx.sigma, fx.f, fx.d,
                int(fx.halfnpw), float(fx.k0), float(fx.kpar))
        try:
            if f["function"] == "RTA":
                out = dict(zip(keys, RTA_1d_tm(*args)))
            elif f["function"] == "SM":
                out = dict(zip(keys, SM_1d_tm(*args)))
            else:
                _, _, Ex, Ez, Sz = field_1d_tm(*args, int(fx.nx), fx.nz)
                out = {"Ex": Ex, "Ez": Ez, "Sz": Sz}
        except Exception:
            continue
        e = 0.0
        for k in keys:
            pv = _vec(out[k]); rv = _vec(getattr(ref, k))
            if pv.size == rv.size and pv.size:
                e = max(e, float(np.max(np.abs(pv - rv))))
        xs.append(ratio); ys.append(float(fx.theta_deg)); es.append(max(e, 1e-18))

    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    sc = ax.scatter(xs, ys, c=np.log10(es), cmap="viridis", s=45, edgecolors="k", linewidths=0.3)
    cb = fig.colorbar(sc, ax=ax); cb.set_label("log10(max abs error) vs MATLAB")
    ax.set_xscale("log")
    ax.set_xlabel("period / wavelength  (a / λ0)")
    ax.set_ylabel("incidence angle θ (deg)")
    ax.set_title("PPML 1-D TM — per-fixture error over (a/λ, θ)")
    ax.grid(True, alpha=0.25)
    fig.text(0.99, 0.01, f"PPML 1-D TM verification {version} · {date}",
             ha="right", va="bottom", fontsize=7, color="#888")
    os.makedirs(figures_dir, exist_ok=True)
    outpath = os.path.join(figures_dir, "error_heatmap.png")
    fig.savefig(outpath, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return outpath


def main():
    figures_dir = str(ROOT / "reports" / "figures")
    date = _dt.date.today().isoformat()
    import json
    manifest = json.loads((ROOT / "fixtures" / "manifest.json").read_text())
    version = f"v{manifest.get('corpus_version', '?')}"

    ref_re, py_re, ref_im, py_im = collect(ROOT, "matlab")
    p1, n, max_rel = scatter_ref_vs_py(ref_re, py_re, ref_im, py_im, version, date, figures_dir)
    print(f"[plots] {p1}  ({n} values, worst rel {max_rel:.2e})")
    p2 = error_histogram(ref_re, py_re, ref_im, py_im, 1e-9, version, date, figures_dir)
    print(f"[plots] {p2}")
    p3 = error_heatmap(ROOT, version, date, figures_dir)
    print(f"[plots] {p3}")


if __name__ == "__main__":
    main()
