# NAME_MAP.md — MATLAB → Python renames & semantic mappings (1-D TM port)

Python package: `python_src/ppml_1d_tm/`. NumPy/SciPy only.

## Function / file map
| MATLAB (`matlab_src/`) | Python (`ppml_1d_tm/`) | Notes |
|---|---|---|
| `1d_tm/SM_1d_tm.m` | `tm_1d.SM_1d_tm` | signature identical order; returns tuple `(rl,rr,tlr,trl)` |
| `1d_tm/RTA_1d_tm.m` | `tm_1d.RTA_1d_tm` | returns `(RR, TT, AA)`, `AA` a length-L `ndarray` |
| `1d_tm/field_1d_tm.m` | `tm_1d.field_1d_tm` | returns `(x, z, Ex, Ez, Sz)` |
| `general/sqrt_whittaker.m` | `general.sqrt_whittaker` | vectorized; identical 4-case branch |
| `general/smpropag_fw_cond.m` | `general.smpropag_fw_cond` | — |
| `general/smpropag_bw_cond.m` | `general.smpropag_bw_cond` | — |
| (inline core, dup in all 3) | `tm_1d._build_stack` | shared eigen+Redheffer core; refactor only, no arithmetic change |
| `error('RTA:EnNotCons',…)` | `tm_1d.PPMLEnergyError` (`.identifier`) | carries the MATLAB id for outcome comparison |

## Operator / idiom map
| MATLAB | Python | Note |
|---|---|---|
| `A\B` (mldivide) | `general.mldivide(A,B)` = `np.linalg.solve(A,B)` | see singular caveat below |
| `A/B` (mrdivide) | `general.mrdivide(A,B)` = `solve(B.T, A.T).T` | avoids forming `inv` |
| `diag(v)` (v vector) | `np.diag(v)` | build diagonal |
| `diag(M)` (M matrix) | `np.diag(M)` | extract diagonal |
| `X'` (ctranspose) | `np.conj(X).T` / `np.conj(x)` | conjugated — RTA flux, field Sz |
| `X.'` (transpose) | `X.T` | non-conjugated — the `hy.'`/`ex.'` in RTA:169 |
| `halfnpw+1` (1-based center) | `halfnpw` (0-based) | zero-order index |
| `n = -halfnpw:halfnpw` | `np.arange(-halfnpw, halfnpw+1)` | length `npw` |
| `sin(pi*f*dn)/(pi*dn)`, diag=`f` | `np.sin(...)/(pi*dn)` + `fill_diagonal(F,f)` | hand-coded Fourier kernel (no toolbox `sinc`) |
| layer axis `(:,:,l)` 1..L+2 | Python list index `l` 0..L+1 | 0-based layers |

## Semantic differences that REQUIRED explicit handling (fidelity-critical)

### D1 — complex auto-promotion of `flux` / `Sz` (bug-for-bug, F2)
MATLAB `flux = zeros(...)` / `Sz = zeros(...)` are real, but AUTO-PROMOTE to complex on
assignment because they are divided by `incflux = q0/k0²`, which is **complex when
`epssup` is complex**. NumPy does NOT auto-promote — a real array silently truncates the
imaginary part. **Fix:** declare `flux` and `Sz` as `dtype=complex`. Caught by
`out_of_domain_epssup_complex_RTA` (real `flux` gave `RR=0.4018`; MATLAB gives
`0.4018+0.0938i`). After the fix, all three epssup-complex fixtures match MATLAB to
< 1.2e-15.

### D2 — singular matrix at an EXACT cutoff: THREE-WAY divergence + explicit policy
At an EXACT Rayleigh/Wood cutoff a mode wavevector `q → 0` and the linear system is
genuinely singular. The three implementations DISAGREE — confirmed on the one fixture
that hits it, `branch_wood_anomaly_+0.00` (dth=0, `near_anomaly=true`):

| impl | outcome | value | warning captured |
|---|---|---|---|
| MATLAB R2025a | `ok` | `RR = NaN` | `MATLAB:illConditionedMatrix` |
| Octave 11.3.0 | `ok` | `RR = 0.7686` (finite garbage) | `Octave:singular-matrix` |
| Python (default) | **error** | — (`LinAlgError`) | — |

The two REFERENCE engines already disagree (NaN vs finite), so **no port choice is
bit-for-bit** — the point is fundamentally indeterminate. `run_reference.m` now records
`lastwarn` (`warn_id`/`warn_msg`) as part of the outcome so this is visible in the data.

**Explicit port policy (flag `on_singular`, threaded through all three public funcs):**
- `'raise'` (**default**) — let `LinAlgError` propagate. Safer for the user than
  returning silent NaN/garbage. This is an intentional deviation from the original's
  limp-along, so its single, self-explaining status on the Python-vs-MATLAB axis is
  **`DEVIATION(flag=on_singular)`** — NOT "INFO", NOT a bare PASS.
- `'replicate'` — best-effort min-norm `lstsq` (returns `0.7686`, matching Octave's
  finite value, NOT MATLAB's NaN). Documented as NOT bit-exact; provided for pipelines
  that must produce a value.
On the MATLAB-vs-Octave (cross-engine) axis the same fixture is **INFO (indeterminate
singular point)** — the references' own NaN-vs-finite disagreement, not a port matter.
Real usage does not sit exactly on a cutoff; the off-cutoff Wood fixtures (`+0.15`,
`-0.50`) port cleanly (< 1e-11). (VERIFICATION_NOTES R2, R8; Scope.)

## Replicated original behaviors (see ORIGINAL_CODE_FINDINGS.md)
- **F1** vacuum impedance `376.730` — hardcoded in `general.Z0_TRUNCATED`; fix behind
  `fix_vacuum_impedance=True` (default OFF → DEVIATION when ON).
- **F2** complex `epssup` → silent nonphysical (all three) — replicated; input guard
  behind `guard_real_epssup=True` (default OFF → DEVIATION when ON).
- **F3** RTA energy guard is a telescoping tautology — replicated as-is (never fires
  for finite inputs); verification uses independent invariants, not this guard.

## Full-corpus self-check (Python vs MATLAB R2025a, all 249 fixtures)
`RTA max 8.3e-12 · SM max 1.9e-13 · field max 4.5e-10` (abs); 248/249 produce values,
1 (`branch_wood_anomaly_+0.00`) raises per D2. Formal gating is Phase 5.
