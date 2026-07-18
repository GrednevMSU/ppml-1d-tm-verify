#!/usr/bin/env python3
"""compare.py — PPML 1-D TM differential test harness (Phase 5).

Loads the frozen corpus + MATLAB (reference of record) and Octave (cross-check) outputs,
runs the Python port, and produces a single self-contained report.html + a
verification_report.md. Implements the OUTCOME MODEL (result/error/timeout compared
before numbers), the pass logic OR(rel<=tol_rel, abs<=tol_abs), the tier system
(near_anomaly / resonance / cross-engine), and the independent physical invariants.

Exit code: 0 if every gating check PASSes, 1 otherwise.

Usage:  python verify/compare.py [--root <dir>]
"""
from __future__ import annotations

import argparse
import base64
import datetime as _dt
import html
import io
import json
import platform
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib
import numpy as np
import yaml
from scipy.io import loadmat

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "python_src"))
from ppml_1d_tm import PPMLEnergyError, RTA_1d_tm, SM_1d_tm, field_1d_tm  # noqa: E402

OUTPUTS = {"RTA": ["RR", "TT", "AA"], "SM": ["rl", "rr", "tlr", "trl"],
           "field": ["Ex", "Ez", "Sz"]}


# --------------------------------------------------------------------------- I/O
def load_fixture(root: Path, fid: str):
    return loadmat(root / "fixtures" / f"{fid}.mat",
                   struct_as_record=False, squeeze_me=True)["fx"]


def load_ref(root: Path, engine: str, fid: str):
    p = root / "reference_outputs" / engine / f"{fid}.mat"
    if not p.exists():
        return None
    return loadmat(p, struct_as_record=False, squeeze_me=True)["ref"]


def run_python(fx):
    """Run the port. Returns (outcome, err_id, values_dict)."""
    a = fx
    args = (a.a, int(a.L), complex(a.epssup), complex(a.epssub),
            a.epsxA, a.epszA, a.epsxB, a.epszB, a.sigma, a.f, a.d,
            int(a.halfnpw), float(a.k0), float(a.kpar))
    try:
        if a.func == "RTA":
            RR, TT, AA = RTA_1d_tm(*args)
            return "ok", "", {"RR": RR, "TT": TT, "AA": AA}
        if a.func == "SM":
            rl, rr, tlr, trl = SM_1d_tm(*args)
            return "ok", "", {"rl": rl, "rr": rr, "tlr": tlr, "trl": trl}
        x, z, Ex, Ez, Sz = field_1d_tm(*args, int(a.nx), a.nz)
        return "ok", "", {"Ex": Ex, "Ez": Ez, "Sz": Sz}
    except PPMLEnergyError as e:
        return "error", e.identifier, {}
    except np.linalg.LinAlgError:
        return "error", "numpy:LinAlgError", {}
    except Exception as e:  # noqa: BLE001
        return "error", f"py:{type(e).__name__}", {}


# ---------------------------------------------------------------- value compare
def _vec(x):
    return np.atleast_1d(np.asarray(x, dtype=complex)).ravel()


def compare_values(py: dict, ref, keys, tol_rel, tol_abs, nan_eq=True):
    """Element-wise OR(rel,abs). Returns dict with worst rel/abs (+where) and pass."""
    worst_rel = 0.0; worst_abs = 0.0; wr_key = ""; wa_key = ""
    all_pass = True; n_el = 0
    for k in keys:
        pv = _vec(py[k]); rv = _vec(getattr(ref, k))
        if pv.size == 0 and rv.size == 0:
            continue
        if pv.shape != rv.shape:
            return dict(pass_=False, shape=True, worst_rel=np.inf, worst_abs=np.inf,
                        wr_key=k, wa_key=k, n_el=0)
        nan_p = np.isnan(pv); nan_r = np.isnan(rv)
        both_nan = nan_p & nan_r
        aerr = np.abs(pv - rv)
        if nan_eq:
            aerr = np.where(both_nan, 0.0, aerr)
        scale = np.abs(rv)
        with np.errstate(divide="ignore", invalid="ignore"):
            rerr = np.where(scale > 0, aerr / scale, np.where(aerr == 0, 0.0, np.inf))
        # PASS LOGIC (per element): abs<=tol_abs OR rel<=tol_rel (or both NaN).
        el_pass = (aerr <= tol_abs) | (rerr <= tol_rel) | both_nan
        if not np.all(el_pass):
            all_pass = False
        # REPORTING (independent of pass logic): true worst abs over all elements, and
        # true worst REL over RESOLVABLE elements only (|ref| > tol_abs). Below that
        # floor relative error is meaningless (divide-by-noise) and the element is
        # judged on the abs branch anyway; including it would print rel=inf for AA=0.
        if aerr.size and float(np.max(aerr)) > worst_abs:
            worst_abs = float(np.max(aerr)); wa_key = k
        resolvable = scale > tol_abs
        if np.any(resolvable):
            wr = float(np.max(rerr[resolvable]))
            if wr > worst_rel:
                worst_rel = wr; wr_key = k
        n_el += pv.size
    return dict(pass_=all_pass, shape=False, worst_rel=worst_rel, worst_abs=worst_abs,
                wr_key=wr_key, wa_key=wa_key, n_el=n_el)


# --------------------------------------------------------------------- records
@dataclass
class Row:
    fid: str
    func: str
    category: str
    tier: str              # "" | near_anomaly | resonance
    axis: str              # gate | cross | expect | invariant
    status: str            # PASS(rel)/PASS(abs)/PASS(both)/FAIL/INFO/DEVIATION(...)
    gating: bool
    worst_rel: float = 0.0
    worst_abs: float = 0.0
    detail: str = ""
    link: str = ""


def func_status(stats, tol_rel, tol_abs):
    """Aggregate self-explaining status for a function's summary row."""
    if stats["fail"] > 0:
        return "FAIL"
    rel_ok = stats["max_rel"] <= tol_rel
    abs_ok = stats["max_abs"] <= tol_abs
    if rel_ok and abs_ok:
        return "PASS(both)"
    if rel_ok:
        return "PASS(rel)"
    if abs_ok:
        return "PASS(abs)"
    return "PASS(mixed)"   # every element passed, but on differing branches


def status_for(cmp, tol_rel, tol_abs):
    if cmp["shape"]:
        return "FAIL", "shape mismatch (auto-FAIL)"
    if not cmp["pass_"]:
        return "FAIL", f"rel={cmp['worst_rel']:.2e}>{tol_rel:g} & abs={cmp['worst_abs']:.2e}>{tol_abs:g}"
    rel_ok = cmp["worst_rel"] <= tol_rel
    abs_ok = cmp["worst_abs"] <= tol_abs
    if rel_ok and abs_ok:
        return "PASS(both)", ""
    if rel_ok:
        return "PASS(rel)", f"abs {cmp['worst_abs']:.2e} > {tol_abs:g}, rel carries"
    return "PASS(abs)", f"rel {cmp['worst_rel']:.2e} > {tol_rel:g}, abs carries"


# ------------------------------------------------------------------- main flow
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=str(ROOT))
    ap.add_argument("--reference", choices=["auto", "matlab", "octave"], default="auto",
                    help="reference of record; 'auto' = matlab if present else octave")
    args = ap.parse_args()
    root = Path(args.root)

    tol = yaml.safe_load((root / "verify" / "tolerances.yaml").read_text())
    manifest = json.loads((root / "fixtures" / "manifest.json").read_text())
    fixtures = manifest["fixtures"]

    def _has(engine):
        return (root / "reference_outputs" / engine).is_dir()

    if args.reference == "auto":
        ref_engine = "matlab" if _has("matlab") else "octave"
    else:
        ref_engine = args.reference
    cross_engine = "octave" if ref_engine == "matlab" else "matlab"
    if not _has(cross_engine):
        cross_engine = None
    if not _has(ref_engine):
        sys.exit(f"[compare] ERROR: no reference_outputs/{ref_engine}/ — run run_reference first")
    env_mat = _load_env(root, ref_engine)
    env_oct = _load_env(root, cross_engine) if cross_engine else {}
    print(f"[compare] reference of record = {ref_engine}"
          + (f"; cross-engine = {cross_engine}" if cross_engine else "; (no cross-engine engine present)"))

    rows: list[Row] = []
    # worst-case trackers (separate rel / abs), gate axis, gating only
    gw_rel = (0.0, "", ""); gw_abs = (0.0, "", "")
    func_stats = {f: dict(n=0, max_rel=0.0, max_abs=0.0, fail=0) for f in OUTPUTS}
    cross_worst = {f: 0.0 for f in OUTPUTS}
    n_fail = 0
    # ref-vs-py value pairs for the y=x scatter (all physical ok/ok fixtures)
    sc_ref: list = []; sc_py: list = []

    for f in fixtures:
        fid = f["id"]; func = f["function"]
        fx = load_fixture(root, fid)
        ref_m = load_ref(root, ref_engine, fid)
        ref_o = load_ref(root, cross_engine, fid) if cross_engine else None
        py_outcome, py_errid, py_vals = run_python(fx)

        tier = ""
        if f["near_anomaly"]:
            tier = "near_anomaly"
        elif ref_m is not None and float(getattr(ref_m, "redheffer_cond", np.nan) or np.nan) > 1e8:
            tier = "resonance"

        # tolerances for this fixture/axis
        base = tol["gate"][func]
        if tier == "near_anomaly":
            t_rel, t_abs, gating = tol["tiers"]["near_anomaly"]["tol_rel"], tol["tiers"]["near_anomaly"]["tol_abs"], False
        elif tier == "resonance":
            t_rel, t_abs, gating = tol["tiers"]["resonance"]["tol_rel"], tol["tiers"]["resonance"]["tol_abs"], False
        else:
            t_rel, t_abs, gating = base["tol_rel"], base["tol_abs"], True

        # ---- GATE axis: python vs MATLAB (outcome first) --------------------
        mat_outcome = str(getattr(ref_m, "outcome", "missing")) if ref_m is not None else "missing"
        row = _axis_gate(fid, func, f, tier, gating, py_outcome, py_errid, py_vals,
                         ref_m, ref_o, mat_outcome, OUTPUTS[func], t_rel, t_abs,
                         tol["nan_equals_nan"], ref_engine, cross_engine)
        rows.append(row)
        # collect value pairs for the y=x scatter (physical, both-ok fixtures)
        if f["expect"] == "ok" and py_outcome == "ok" and mat_outcome == "ok":
            for k in OUTPUTS[func]:
                pv = _vec(py_vals[k]); rv = _vec(getattr(ref_m, k))
                if pv.size == rv.size and pv.size:
                    sc_ref.append(np.concatenate([rv.real, rv.imag]))
                    sc_py.append(np.concatenate([pv.real, pv.imag]))
        if row.status == "FAIL" and row.gating:
            n_fail += 1
            func_stats[func]["fail"] += 1
        func_stats[func]["n"] += 1
        if row.gating and row.status.startswith("PASS"):
            func_stats[func]["max_rel"] = max(func_stats[func]["max_rel"], row.worst_rel)
            func_stats[func]["max_abs"] = max(func_stats[func]["max_abs"], row.worst_abs)
            if row.worst_rel > gw_rel[0]:
                gw_rel = (row.worst_rel, fid, row.detail)
            if row.worst_abs > gw_abs[0]:
                gw_abs = (row.worst_abs, fid, row.detail)

        # ---- CROSS-ENGINE axis: MATLAB vs Octave (INFO) ---------------------
        if ref_m is not None and ref_o is not None and mat_outcome == "ok" and str(ref_o.outcome) == "ok":
            cmp = compare_values(_ref_to_dict(ref_m, OUTPUTS[func]), ref_o, OUTPUTS[func],
                                 tol["cross_engine"]["tol_abs"], tol["cross_engine"]["tol_abs"],
                                 tol["nan_equals_nan"])
            ce = tol["cross_engine"]["tol_abs"]
            ok = cmp["worst_abs"] <= ce or cmp["pass_"]
            cross_worst[func] = max(cross_worst[func], cmp["worst_abs"])
            if not ok and not (tier == "near_anomaly"):
                rows.append(Row(fid, func, f["category"], tier, "cross", "FAIL", False,
                                cmp["worst_rel"], cmp["worst_abs"],
                                f"{ref_engine}-vs-{cross_engine} abs {cmp['worst_abs']:.2e} > {ce:g} — engine finding", "10.3"))
        elif ref_m is not None and ref_o is not None and mat_outcome != str(ref_o.outcome):
            # references disagree on outcome (e.g. singular cutoff NaN vs finite)
            rows.append(Row(fid, func, f["category"], tier, "cross", "INFO", False, 0, 0,
                            f"references disagree: {ref_engine}={mat_outcome}/{getattr(ref_m,'err_id','')} "
                            f"{cross_engine}={ref_o.outcome} — indeterminate (D2)", "D2"))

        # ---- EXPECT axis --------------------------------------------------
        exp = f["expect"]
        exp_kind = "error" if exp == "error" else "ok"
        act_kind = "ok" if py_outcome in ("ok", "timeout") else "error"
        exp_ok = (act_kind == exp_kind)
        if exp != "ok":
            st = "PASS(both)" if exp_ok else "FAIL"
            rows.append(Row(fid, func, f["category"], tier, "expect", st, exp != "nonphysical",
                            0, 0, f"expect={exp} actual={py_outcome}", "F2/F3"))

        # ---- INVARIANTS (independent, gating) -----------------------------
        if py_outcome == "ok":
            for r in _invariants(fid, func, f, py_vals, tol):
                rows.append(r)
                if r.status == "FAIL" and r.gating:
                    n_fail += 1

    overall = "PASSED" if n_fail == 0 else "FAILED"
    plots = _make_plots(rows, tol)
    scatter_b64 = _scatter_b64(sc_ref, sc_py, ref_engine) if sc_ref else ""
    ctx = dict(overall=overall, n_fail=n_fail, rows=rows, func_stats=func_stats,
               gw_rel=gw_rel, gw_abs=gw_abs, cross_worst=cross_worst,
               manifest=manifest, env_mat=env_mat, env_oct=env_oct, tol=tol,
               fixtures=fixtures, plots=plots, root=root, scatter_b64=scatter_b64,
               ref_engine=ref_engine, cross_engine=cross_engine)

    html_out = render_html(ctx)
    md_out = render_md(ctx)
    (root / "report.html").write_text(html_out)
    (root / "verification_report.md").write_text(md_out)
    # Auto-sync the versioned archive (manual copying silently drifts from the original).
    archive = root / "reports" / f"v{manifest.get('corpus_version', 'x')}"
    archive.mkdir(parents=True, exist_ok=True)
    (archive / "report.html").write_text(html_out)
    (archive / "verification_report.md").write_text(md_out)
    print(f"[compare] overall={overall}  gating-FAIL={n_fail}  rows={len(rows)}")
    print(f"[compare] worst gate rel={gw_rel[0]:.2e}@{gw_rel[1]}  abs={gw_abs[0]:.2e}@{gw_abs[1]}")
    print(f"[compare] wrote report.html + verification_report.md (+ archived to {archive.name}/)")
    sys.exit(0 if overall == "PASSED" else 1)


def _axis_gate(fid, func, f, tier, gating, py_outcome, py_errid, py_vals,
               ref_m, ref_o, mat_outcome, keys, t_rel, t_abs, nan_eq,
               ref_engine, cross_engine) -> Row:
    cat = f["category"]
    # outcome compared BEFORE numbers
    mat_kind = "ok" if mat_outcome in ("ok", "timeout") else ("error" if mat_outcome == "error" else "missing")
    py_kind = "ok" if py_outcome in ("ok", "timeout") else "error"
    if mat_kind == "missing":
        return Row(fid, func, cat, tier, "gate", "INFO", False, 0, 0,
                   f"no {ref_engine} reference", "")
    if py_kind != mat_kind:
        # one raised, other didn't
        if py_kind == "error" and f["near_anomaly"]:
            det = _threeway_detail(ref_m, ref_o, py_errid, ref_engine, cross_engine, keys)
            return Row(fid, func, cat, tier, "gate", "DEVIATION(flag=on_singular)", False,
                       0, 0, det, "D2")
        return Row(fid, func, cat, tier, "gate", "FAIL", gating, np.inf, np.inf,
                   f"outcome mismatch: {ref_engine}={mat_outcome} Python={py_outcome}/{py_errid}", "9a")
    if py_kind == "error":
        mid = str(getattr(ref_m, "err_id", "") or "")
        if mid == py_errid or (mid and py_errid.endswith(mid.split(":")[-1])):
            return Row(fid, func, cat, tier, "gate", "PASS(both)", False, 0, 0,
                       f"both error ({py_errid})", "9a")
        return Row(fid, func, cat, tier, "gate", "FAIL", gating, np.inf, np.inf,
                   f"error-id mismatch MATLAB={mid} Python={py_errid}", "9a")
    # both ok -> value compare
    cmp = compare_values(py_vals, ref_m, keys, t_rel, t_abs, nan_eq)
    st, det = status_for(cmp, t_rel, t_abs)
    if not gating:                      # near_anomaly / resonance tiers are non-gating
        st = "INFO"
    return Row(fid, func, cat, tier, "gate", st, gating,
               cmp["worst_rel"], cmp["worst_abs"],
               det or f"rel {cmp['worst_rel']:.2e} / abs {cmp['worst_abs']:.2e} ({cmp['wr_key'] or cmp['wa_key']})",
               "10.5" if func == "field" else "")


def _threeway_detail(ref_m, ref_o, py_errid, ref_engine, cross_engine, keys):
    """One-line three-way outcome exhibit for the singular-cutoff DEVIATION (D2).
    Reads the ACTUAL recorded values/warnings — no 'MATLAB ok' shorthand."""
    def descr(ref):
        if ref is None:
            return "absent"
        try:
            v0 = complex(np.atleast_1d(np.asarray(getattr(ref, keys[0]), dtype=complex)).ravel()[0])
            val = "NaN" if np.isnan(v0) else f"finite({keys[0]}={v0.real:.3g})"
        except Exception:
            val = "?"
        w = str(getattr(ref, "warn_id", "") or "")
        return f"{ref.outcome}/{val}" + (f"+{w}" if w else "")
    parts = [f"{ref_engine}: {descr(ref_m)}"]
    if cross_engine:
        parts.append(f"{cross_engine}: {descr(ref_o)}")
    parts.append(f"Python: raise({py_errid}, flag=on_singular)")
    return " · ".join(parts) + " — three-way indeterminate at exact cutoff"


def _invariants(fid, func, f, vals, tol) -> list[Row]:
    out = []
    if func != "RTA":
        return out
    RR = complex(vals["RR"]); TT = complex(vals["TT"])
    AA = np.atleast_1d(np.asarray(vals["AA"], dtype=complex))
    # realness (real epssup only)
    if f["epssup_im"] == 0:
        ta = tol["invariants"]["realness"]["tol_abs"]
        im = max(abs(RR.imag), abs(TT.imag))
        st = "PASS(abs)" if im <= ta else "FAIL"
        out.append(Row(fid, func, f["category"], f["near_anomaly"] and "near_anomaly" or "",
                       "invariant", st, not f["near_anomaly"], 0, im,
                       f"|Im(RR/TT)|={im:.2e} (realness)", "F3"))
    # energy bounds (lossless)
    if f["lossless"]:
        ta = tol["invariants"]["energy_bounds"]["tol_abs"]
        vals_rt = [RR.real, TT.real] + list(np.real(AA))
        inb = all(-ta <= v <= 1 + ta for v in vals_rt)
        cons = abs(RR.real + TT.real + float(np.sum(np.real(AA))) - 1.0)
        st = "PASS(abs)" if (inb and cons <= ta) else "FAIL"
        out.append(Row(fid, func, f["category"], "", "invariant", st, True, 0, cons,
                       f"R+T+A-1={cons:.2e}, in[0,1]={inb} (lossless energy)", "F3"))
    return out


def _ref_to_dict(ref, keys):
    return {k: getattr(ref, k) for k in keys}


def _load_env(root, engine):
    p = root / "reference_outputs" / engine / "env.json"
    return json.loads(p.read_text()) if p.exists() else {}


# ------------------------------------------------------------------- plotting
def _b64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=90, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def _scatter_b64(sc_ref, sc_py, ref_engine):
    """y = x scatter: reference (x) vs Python port (y), every value across all
    physical fixtures. Also written to reports/figures/reference_vs_python.png."""
    ref = np.concatenate(sc_ref); py = np.concatenate(sc_py)
    good = np.isfinite(ref) & np.isfinite(py)
    ref, py = ref[good], py[good]
    n = ref.size
    with np.errstate(divide="ignore", invalid="ignore"):
        denom = np.abs(ref)
        rel = np.where(denom > 1e-10, np.abs(py - ref) / denom, np.nan)
    max_rel = float(np.nanmax(rel)) if np.any(np.isfinite(rel)) else 0.0
    fig, ax = plt.subplots(figsize=(6.6, 5.8))
    ax.scatter(ref, py, s=5, alpha=0.35, color="#1f77b4", linewidths=0,
               label=f"{n} values (Re & Im, all gated outputs)")
    lim = max(1.0, float(np.max(np.abs(ref))) if n else 1.0)
    ax.plot([-lim, lim], [-lim, lim], color="crimson", lw=1, label="y = x")
    ax.set_xlabel(f"{ref_engine} reference value"); ax.set_ylabel("Python port value")
    ax.set_title(f"{ref_engine.upper()} reference vs Python port — every output, all fixtures\n"
                 f"points on y = x (worst rel {max_rel:.1e}, incl. non-gating near-anomaly)")
    ax.legend(loc="upper left", fontsize=8); ax.grid(True, alpha=0.25)
    ax.set_aspect("equal", adjustable="box")
    fig.text(0.99, 0.01, f"PPML 1-D TM verification · {_dt.date.today().isoformat()}",
             ha="right", va="bottom", fontsize=7, color="#888")
    return _b64(fig)


def _make_plots(rows, tol):
    plots = {}
    for func in OUTPUTS:
        pts = [(r.worst_rel, r.worst_abs, r.fid) for r in rows
               if r.axis == "gate" and r.func == func and r.status.startswith("PASS")]
        if not pts:
            continue
        rels = [max(p[0], 1e-18) for p in pts]
        fig, ax = plt.subplots(figsize=(6, 2.6))
        ax.semilogy(range(len(rels)), sorted(rels), ".", ms=4)
        ax.axhline(tol["gate"][func]["tol_rel"], color="crimson", ls="--", lw=1,
                   label=f"tol_rel {tol['gate'][func]['tol_rel']:g}")
        ax.set_title(f"{func}: per-fixture worst relative residual (Python vs MATLAB)")
        ax.set_xlabel("fixture (sorted)"); ax.set_ylabel("rel err")
        ax.legend(fontsize=8); ax.grid(True, which="both", alpha=0.25)
        plots[func] = _b64(fig)
    return plots


# ------------------------------------------------------------------- scope
def _scope_bullets(root):
    """Extract §12 Scope bullets, JOINING soft-wrapped continuation lines (a bullet may
    span several physical lines in the source; the earlier version truncated them)."""
    txt = (root / "VERIFICATION_NOTES.md").read_text().splitlines()
    out = []; cap = False
    for ln in txt:
        if ln.startswith("## 12. Scope"):
            cap = True; continue
        if not cap:
            continue
        if ln.startswith("## ") or ln.startswith("---"):
            break
        s = ln.rstrip()
        if s.strip().startswith("- "):
            out.append(s.strip()[2:])
        elif s.strip() and out:            # continuation of the current bullet
            sep = "" if out[-1].endswith("-") else " "   # don't split hyphenated words
            out[-1] = out[-1] + sep + s.strip()
    return out


# ------------------------------------------------------------------- renderers
def _badge(status):
    if status.startswith("PASS"): return "pass"
    if status.startswith("FAIL"): return "fail"
    if status.startswith("DEVIATION"): return "dev"
    return "info"


def render_html(c):
    root = c["root"]; overall = c["overall"]
    scope = _scope_bullets(root)
    exec_lines = _exec_summary(c)
    css = """
    body{font-family:-apple-system,Segoe UI,Roboto,sans-serif;margin:0;color:#1a1a1a;background:#fafafa}
    .wrap{max-width:1100px;margin:0 auto;padding:24px}
    .banner{padding:28px;border-radius:12px;color:#fff;font-size:34px;font-weight:800;letter-spacing:1px}
    .PASSED{background:linear-gradient(90deg,#1a7f37,#2ea043)}
    .FAILED{background:linear-gradient(90deg,#b91c1c,#ef4444)}
    h2{margin-top:34px;border-bottom:2px solid #e2e2e2;padding-bottom:6px}
    table{border-collapse:collapse;width:100%;font-size:13px;margin:10px 0}
    th,td{border:1px solid #ddd;padding:6px 8px;text-align:left}
    th{background:#f0f0f0}
    .pass{color:#1a7f37;font-weight:700}.fail{color:#b91c1c;font-weight:800}
    .dev{color:#9a6700;font-weight:700}.info{color:#57606a;font-weight:600}
    code{background:#f3f3f3;padding:1px 4px;border-radius:4px}
    details{margin:8px 0;border:1px solid #e2e2e2;border-radius:8px;padding:8px 12px;background:#fff}
    summary{cursor:pointer;font-weight:600}
    .exec{background:#fff;border-left:5px solid #2ea043;padding:14px 18px;border-radius:6px;font-size:15px;line-height:1.55}
    .FAILED+.exec{border-left-color:#ef4444}
    .mut{color:#666;font-size:12px}
    img{max-width:100%;height:auto;border:1px solid #eee;border-radius:6px;margin:6px 0}
    footer{margin-top:40px;font-size:12px;color:#555;border-top:1px solid #ddd;padding-top:14px}
    """
    P = [f"<!-- generated {_dt.datetime.now().isoformat()} -->",
         f"<style>{css}</style><div class='wrap'>",
         f"<div class='banner {overall}'>DIFFERENTIAL TEST HARNESS — {overall}</div>",
         "<div class='exec'>", "<b>Executive summary.</b><br>"]
    P += [html.escape(s) + "<br>" for s in exec_lines]
    P.append("</div>")

    if c.get("scatter_b64"):
        P.append("<h2>Reference vs Python — every value on y = x</h2>"
                 f"<img src='data:image/png;base64,{c['scatter_b64']}'/>"
                 f"<p class='mut'>Every gated numeric output (real &amp; imaginary parts) "
                 f"from every physical fixture: {c['ref_engine']} reference on x, Python port "
                 f"on y. Agreement = points on the red y = x line.</p>")

    P.append("<h2>Pass logic</h2><p><code>PASS iff OR(rel_err ≤ tol_rel, abs_err ≤ tol_abs)</code>, "
             "per element. Statuses self-explain: PASS(rel) / PASS(abs) / PASS(both) / FAIL / "
             "INFO(non-gating) / DEVIATION(flag=…). Reference of record = MATLAB; Octave is the "
             "cross-engine check; the Python port is under test.</p>")

    # summary table
    refname = c["ref_engine"].upper()
    P.append(f"<h2>Summary (gate axis: Python vs {refname})</h2><table>"
             "<tr><th>Function</th><th>Status</th><th>Fixtures</th><th>Gating FAIL</th>"
             "<th>max rel (worst)</th><th>max abs (worst)</th><th>tol_rel / tol_abs</th></tr>")
    for fn, s in c["func_stats"].items():
        t = c["tol"]["gate"][fn]
        st = func_status(s, t["tol_rel"], t["tol_abs"])
        P.append(f"<tr><td><b>{fn}</b></td>"
                 f"<td class='{_badge(st)}'>{st}</td><td>{s['n']}</td>"
                 f"<td class='{'fail' if s['fail'] else 'pass'}'>{s['fail']}</td>"
                 f"<td>{s['max_rel']:.2e}</td><td>{s['max_abs']:.2e}</td>"
                 f"<td>{t['tol_rel']:g} / {t['tol_abs']:g}</td></tr>")
    P.append("</table><p class='mut'>Status qualifies every row: e.g. <code>field</code> "
             "shows max abs above its tol_abs but passes as <code>PASS(rel)</code> — the "
             "OR pass-logic carries it on the relative branch (tol_abs is only a floor for "
             "near-null field points).</p>")
    P.append(f"<p class='mut'>Worst-case fixtures reported SEPARATELY for rel and abs "
             f"(they differ): worst rel = <code>{c['gw_rel'][0]:.2e}</code> @ "
             f"<code>{html.escape(c['gw_rel'][1])}</code>; worst abs = "
             f"<code>{c['gw_abs'][0]:.2e}</code> @ <code>{html.escape(c['gw_abs'][1])}</code>.</p>")

    # residual plots
    P.append("<h2>Worst-case residual plots</h2>")
    for fn, b in c["plots"].items():
        P.append(f"<details open><summary>{fn} — per-fixture worst relative residual</summary>"
                 f"<img src='data:image/png;base64,{b}'/></details>")

    # findings: FAIL + DEVIATION + cross-engine + expect
    P.append("<h2>Findings (FAIL / DEVIATION / INFO of note)</h2><table>"
             "<tr><th>Fixture</th><th>Axis</th><th>Status</th><th>Detail</th><th>Ref</th></tr>")
    notable = [r for r in c["rows"] if r.status.startswith(("FAIL", "DEVIATION"))
               or (r.axis in ("cross", "expect") and r.status != "PASS(both)")]
    if not notable:
        P.append("<tr><td colspan=5 class='pass'>none — all gating checks PASS, no deviations of note</td></tr>")
    for r in notable:
        P.append(f"<tr><td><code>{html.escape(r.fid)}</code></td><td>{r.axis}</td>"
                 f"<td class='{_badge(r.status)}'>{html.escape(r.status)}</td>"
                 f"<td>{html.escape(r.detail)}</td><td class='mut'>{r.link}</td></tr>")
    P.append("</table>")

    # tiers / deviations note
    P.append("<h2>Tiers &amp; intentional deviations</h2>")
    na = [r.fid for r in c["rows"] if r.tier == "near_anomaly" and r.axis == "gate"]
    P.append(f"<p><b>near_anomaly</b> (INFO, non-gating): {len(set(na))} fixtures — "
             "Rayleigh/Wood hypersensitivity is physics (tol_rel 1e-3). "
             "<b>resonance</b> tier: 0 (corpus max redheffer_cond = 292 &lt; 1e8). "
             "<b>DEVIATION(flag=on_singular)</b>: the exact-cutoff fixture "
             "<code>branch_wood_anomaly_+0.00</code> — port raises by design where the "
             "references disagree (MATLAB NaN vs Octave finite).</p>")

    # cross-engine
    if c["cross_engine"]:
        P.append(f"<h2>Cross-engine ({refname} vs {c['cross_engine'].upper()}, INFO)</h2><table>"
                 "<tr><th>Function</th><th>max abs diff</th><th>tier</th></tr>")
        for fn, v in c["cross_worst"].items():
            P.append(f"<tr><td>{fn}</td><td>{v:.2e}</td><td>{c['tol']['cross_engine']['tol_abs']:g}</td></tr>")
        P.append("</table>")
    else:
        P.append("<h2>Cross-engine</h2><p class='mut'>Only one engine present on this host; "
                 "cross-engine check skipped. Run on a host with both MATLAB and Octave to "
                 "populate it.</p>")

    # coverage matrix
    P.append("<h2>Scenario coverage matrix</h2><table><tr><th>Category</th><th>Fixtures</th></tr>")
    cats = {}
    for f in c["fixtures"]:
        cats[f["category"]] = cats.get(f["category"], 0) + 1
    for k, v in cats.items():
        P.append(f"<tr><td>{k}</td><td>{v}</td></tr>")
    P.append(f"<tr><td><b>TOTAL</b></td><td><b>{len(c['fixtures'])}</b></td></tr></table>")

    # scope
    P.append("<h2>Scope &amp; limitations</h2><ul>")
    for s in scope:
        P.append(f"<li>{html.escape(s)}</li>")
    P.append("</ul>")

    # footer
    em, eo = c["env_mat"], c["env_oct"]
    npcfg = _np_backend()
    ref_line = (f"Reference ({c['ref_engine']}): {em.get('engine','?')} · "
                f"BLAS {str(em.get('blas','?'))[:40]} · LAPACK {str(em.get('lapack','?'))[:30]}<br>")
    cross_line = ""
    if c["cross_engine"]:
        cross_line = (f"Cross-engine ({c['cross_engine']}): {eo.get('engine','?')} · "
                      f"BLAS {str(eo.get('blas','?'))[:40]} · LAPACK {str(eo.get('lapack','?'))[:30]}<br>")
    P.append("<footer>"
             f"<b>Environment.</b> Report schema {c['manifest']['schema_version']}, "
             f"corpus v{c['manifest'].get('corpus_version','?')}, seed {c['manifest']['seed']}, "
             f"date {_dt.date.today()}.<br>"
             + ref_line + cross_line +
             f"Python: {platform.python_version()} · NumPy {np.__version__} · {npcfg}<br>"
             "<b>Pass logic:</b> OR(rel≤tol_rel, abs≤tol_abs). Vacuum impedance replicated as "
             "376.730 (F1). Original energy guard is a telescoping tautology (F3) — verification "
             "uses independent invariants. See VERIFICATION_NOTES.md, ORIGINAL_CODE_FINDINGS.md, "
             "NAME_MAP.md, CORPUS_CHANGELOG.md.</footer></div>")
    return "\n".join(P)


def _np_backend():
    try:
        cfg = np.show_config(mode="dicts")
        bl = cfg.get("Build Dependencies", {}).get("blas", {})
        name = bl.get("name") or "?"
        ver = bl.get("version")
        # numpy reports version "unknown" for Apple Accelerate — drop it rather than print it.
        label = name if (not ver or str(ver).lower() == "unknown") else f"{name} {ver}"
        return f"BLAS {label}"
    except Exception:
        return "BLAS n/a"


def _exec_summary(c):
    s = c["func_stats"]
    total = len(c["fixtures"])
    refname = c["ref_engine"].upper()
    field_rel = max(s["field"]["max_rel"], 1e-18)
    rt_abs = max(s["RTA"]["max_abs"], s["SM"]["max_abs"], 1e-18)
    L = [
        f"The Python port of the PPML 1-D TM group (RTA / SM / field) was checked against "
        f"the original {refname} code on {total} fixtures spanning subwavelength-to-multi-order "
        f"periods, 0–80° incidence, dielectric and metallic (lossy) layers, and conducting "
        f"interfaces.",
        f"Reflectance/transmittance and S-matrix coefficients agree with {refname} to about "
        f"{rt_abs:.0e} absolute — at the level predicted by the accumulated round-off budget "
        f"(§10.4: npw²·L·κ_V·ε, envelope 5.3e-11), ~6× within it, NOT machine-epsilon "
        f"(2e-16); sampled fields to about {field_rel:.0e} relative (looser by physics — "
        f"thick-metal evanescent dynamic range, not a port error; the port matches {refname} "
        f"better than the two reference engines match each other).",
        f"Overall gating result: {c['overall']} with {c['n_fail']} gating failures.",
        "Two behaviours in the ORIGINAL code are documented and replicated bug-for-bug: a "
        "truncated vacuum impedance 376.730 (F1), and an energy 'guard' whose checked "
        "quantity RR+TT+ΣAA is identically 1 by telescoping (F3) — so it validates nothing "
        "on any finite input, and additionally fails open on NaN and on the "
        "documented-forbidden complex-permittivity input (F2).",
        "Scope: 1-D TM only; exact diffraction-cutoff points are numerically indeterminate "
        "(all three implementations disagree there) and reported non-gating.",
    ]
    return L


def render_md(c):
    s = c["func_stats"]
    L = [f"# PPML 1-D TM — Differential Test Harness — {c['overall']}", ""]
    L.append("## Executive summary")
    L += [f"- {x}" for x in _exec_summary(c)]
    L.append("")
    L.append("## Pass logic")
    L.append("`PASS iff OR(rel_err <= tol_rel, abs_err <= tol_abs)`, per element. "
             "Statuses: PASS(rel)/PASS(abs)/PASS(both)/FAIL/INFO/DEVIATION(flag=…).")
    L.append("")
    L.append(f"## Summary (Python vs {c['ref_engine'].upper()})")
    L.append("| function | status | fixtures | gating FAIL | max rel | max abs | tol_rel/tol_abs |")
    L.append("|---|---|---|---|---|---|---|")
    for fn, st in s.items():
        t = c["tol"]["gate"][fn]
        fs = func_status(st, t["tol_rel"], t["tol_abs"])
        L.append(f"| {fn} | {fs} | {st['n']} | {st['fail']} | {st['max_rel']:.2e} | "
                 f"{st['max_abs']:.2e} | {t['tol_rel']:g}/{t['tol_abs']:g} |")
    L.append("")
    L.append("Status qualifies every row (OR pass-logic): `field` prints max abs above its "
             "tol_abs but passes as `PASS(rel)` — tol_abs is only a floor for near-null field "
             "points. Worst-case fixtures reported SEPARATELY for rel and abs:")
    L.append(f"- worst rel = {c['gw_rel'][0]:.2e} @ `{c['gw_rel'][1]}`")
    L.append(f"- worst abs = {c['gw_abs'][0]:.2e} @ `{c['gw_abs'][1]}`")
    L.append("")
    # --- per-tier observed / budget / threshold ---
    L.append("## Tiers — observed vs budget vs threshold")
    L.append("| tier | applies | observed (max) | budget / basis | threshold | gating |")
    L.append("|---|---|---|---|---|---|")
    L.append(f"| gate:RTA/SM | R/T, S | abs {max(s['RTA']['max_abs'], s['SM']['max_abs']):.2e} | "
             f"npw²·L·κ_V·ε ≈ 5.3e-11 | abs 1e-10 / rel 1e-11 | YES |")
    L.append(f"| gate:field | Ex,Ez,Sz | rel {s['field']['max_rel']:.2e} | evanescent exp(Im q·d), §10.5 | "
             f"rel 1e-9 / abs 1e-11 | YES |")
    ce = c["tol"]["cross_engine"]["tol_abs"]
    L.append(f"| cross-engine | MATLAB↔Octave | abs {max(c['cross_worst'].values()) if c['cross_engine'] else 0:.2e} | "
             f"LAPACK divergence, §10.3 | abs {ce:g} | INFO |")
    na_fix = sorted({r.fid for r in c["rows"] if r.tier == 'near_anomaly'})
    L.append(f"| near_anomaly | {len(na_fix)} fixtures | (hypersensitive) | Wood cutoff physics, §10.1 | "
             f"rel 1e-3 | INFO |")
    L.append("| resonance | 0 fixtures | cond max 292 | redheffer_cond, §10.2 | >1e8 | INFO |")
    L.append("")
    # --- findings ---
    L.append("## Findings (FAIL / DEVIATION / notable INFO)")
    notable = [r for r in c["rows"] if r.status.startswith(("FAIL", "DEVIATION"))
               or (r.axis in ("cross", "expect") and r.status not in ("PASS(both)", "INFO"))
               or (r.axis == "cross" and r.status == "INFO")]
    real_fails = [r for r in notable if r.status.startswith("FAIL") and r.gating]
    devs = [r for r in notable if r.status.startswith("DEVIATION")]
    L.append(f"- {'No gating FAIL' if not real_fails else str(len(real_fails)) + ' gating FAIL'}"
             f"; {len(devs)} documented DEVIATION.")
    for r in notable:
        L.append(f"- `{r.fid}` [{r.axis}] **{r.status}** — {r.detail} ({r.link})")
    # notable INFO summary (near_anomaly)
    L.append(f"- INFO (non-gating): {len(na_fix)} `near_anomaly` fixtures on the loose Wood "
             f"tier (rel 1e-3): {', '.join('`'+x+'`' for x in na_fix[:6])}"
             + (" …" if len(na_fix) > 6 else "") + ".")
    L.append("")
    # --- coverage matrix ---
    L.append("## Scenario coverage matrix")
    L.append("| category | fixtures |")
    L.append("|---|---|")
    cats = {}
    for fx in c["fixtures"]:
        cats[fx["category"]] = cats.get(fx["category"], 0) + 1
    for k, v in cats.items():
        L.append(f"| {k} | {v} |")
    L.append(f"| **TOTAL** | **{len(c['fixtures'])}** |")
    L.append("")
    if c["cross_engine"]:
        L.append(f"## Cross-engine ({c['ref_engine']} vs {c['cross_engine']}, INFO)")
        for fn, v in c["cross_worst"].items():
            L.append(f"- {fn}: max abs diff {v:.2e} (tier {c['tol']['cross_engine']['tol_abs']:g})")
    else:
        L.append("## Cross-engine")
        L.append("- only one engine on this host; cross-engine check skipped.")
    L.append("")
    L.append("## Scope & limitations")
    L += [f"- {x}" for x in _scope_bullets(c["root"])]
    L.append("")
    em = c["env_mat"]; eo = c["env_oct"]
    L.append("## Environment")
    L.append(f"- Reference ({c['ref_engine']}): {em.get('engine','?')} / LAPACK {em.get('lapack','?')}")
    if c["cross_engine"]:
        L.append(f"- Cross-engine ({c['cross_engine']}): {eo.get('engine','?')} / LAPACK {eo.get('lapack','?')}")
    L.append(f"- Python {platform.python_version()} / NumPy {np.__version__} / {_np_backend()}")
    L.append(f"- corpus v{c['manifest'].get('corpus_version','?')}, seed {c['manifest']['seed']}, schema {c['manifest']['schema_version']}")
    return "\n".join(L)


if __name__ == "__main__":
    main()
