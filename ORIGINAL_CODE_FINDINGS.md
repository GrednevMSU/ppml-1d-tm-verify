# ORIGINAL_CODE_FINDINGS.md — PPML v3.0, 1-D TM group

Behaviors in the ORIGINAL MATLAB code that a faithful port must decide how to handle.
Default policy: **replicate bug-for-bug**; any fix is gated behind an explicit flag
(default OFF) and marked `DEVIATION` in the report when enabled. Started in Phase 1/2;
extended during Phase 4.

---

## F1 — Vacuum impedance truncated to `376.730` Ω
- **Where:** `SM_1d_tm.m:135,147`, `RTA_1d_tm.m:135,147`, `field_1d_tm.m:126,138`
  (`sigma(l)*376.730/k0`).
- **What:** the free-space impedance Z₀ is hardcoded as `376.730`. The physical value
  is `376.730313668…` Ω (= μ₀c). The literal drops after the 3rd decimal.
- **Impact:** only affects fixtures with `sigma≠0` (conducting interfaces). Relative
  error in Z₀ ≈ `8.3e-7`, which propagates linearly into the interface term. For our
  `sigma_nonzero` fixtures this is well ABOVE the `rel 1e-11` gate, so substituting a
  more precise constant would make the port MISMATCH the original.
- **Decision:** **replicate exactly** — Python hardcodes `376.730`. Do NOT use
  `scipy.constants` Z₀. Fix (use full-precision Z₀) available behind flag
  `fix_vacuum_impedance` (default OFF) → `DEVIATION(flag=fix_vacuum_impedance)`.
- **Status:** replicate. Verified by the `sigma_nonzero` fixtures (50 of them).

## F2 — Forbidden `epssup` complex input silently returns nonphysical (all 3 funcs)
- **Where:** `RTA_1d_tm.m`, `SM_1d_tm.m`, `field_1d_tm.m` — none validate `epssup`.
- **What:** the docstrings state `epssup` must be real ("otherwise the incident waves
  are ill-defined"), but nothing validates it. **EMPIRICALLY CONFIRMED** (Phase-3
  `run_reference.m`, Octave 11.3.0): all three functions **run to completion and return
  a nonphysical result** — `RR`,`TT` come out **complex** (e.g. `RR=0.4018+0.0938i`)
  where a physical reflectance must be real. No warning, no error, from any of them.
- **Hypothesis refuted:** an earlier draft (and the review) hypothesized RTA would
  `throw RTA:EnNotCons` via its energy guard. The OUTCOME MODEL flagged the mismatch:
  RTA returns `ok`. Reason = finding **F3** (the guard is vacuous). This is a worked
  example of hypothesis-vs-data — the harness caught it instead of hard-coding a wrong
  expectation.
- **Detectable signature:** `|Im(RR)|` (and `|Im(TT)|`) ≈ 0.094 for this input vs
  ≈ 0 (machine eps) for every physical fixture. The independent realness invariant
  (VERIFICATION_NOTES §9) uses `|Im(RR)|`, `|Im(TT)|` for real-`epssup` fixtures —
  precisely because the code's own guard does NOT catch this (F3).
- **Decision:** **replicate bug-for-bug** — port returns the same complex values.
  Fixtures `out_of_domain_epssup_complex_{RTA,SM,field}` all pin `expect=nonphysical`
  (`ok` outcome, physically invalid). Optional input guard for all three behind flag
  `guard_real_epssup` (default OFF) → `DEVIATION(flag=guard_real_epssup)` when enabled.
- **Status:** replicate. Confirmed on Octave; MATLAB expected identical (the mechanism
  is pure algebra, F3) — cross-checked when the MATLAB reference runs.

## F3 — RTA energy-conservation "guard" is a telescoping tautology (vacuous check)
- **Where:** `RTA_1d_tm.m:183`, `if abs(sum(AA)+RR+TT-1) > 1e-5, error('RTA:EnNotCons')`.
- **What:** the checked quantity is identically 1 **by construction**, independent of
  physics or inputs:
  `RR = 1-flux(1)`, `TT = flux(L+1)`, `AA(l) = flux(l)-flux(l+1)` (`RTA.m:172-179`), so
  `RR + TT + Σ AA = 1 - flux(1) + flux(L+1) + (flux(1) - flux(L+1)) ≡ 1` (telescoping).
  Confirmed: for the complex-`epssup` fixture the residual is **exactly 0.000e+00**
  even though the result is nonphysical.
- **Impact:** the guard does **not** validate energy conservation and gives false
  assurance. It cannot catch nonphysical inputs (F2), factorization errors, or a
  mis-ported formula that still telescopes. Worse, it does not even fire on `NaN`:
  Phase-3 confirmed that on the exact-Wood-cutoff fixture MATLAB RTA returns `RR=NaN`
  yet completes with outcome `ok` — because `abs(NaN) > 1e-5` evaluates to **false**.
  So the guard catches neither nonphysical inputs nor NaN — it is effectively dead code
  for finite-and-NaN inputs alike, firing only on `Inf` that survives the telescoping.
- **Decision:** replicate the guard as-is (bug-for-bug — it never fires for finite
  inputs anyway), but **do NOT rely on it** for verification. The harness's INDEPENDENT
  energy/physicality invariants (recompute R+T+A from an independent flux path;
  `|Im(RR)|,|Im(TT)| ≈ 0` for real `epssup`; `0 ≤ R,T,A ≤ 1` for lossless) are what
  actually guard energy — VERIFICATION_NOTES §9. This is exactly the "invariant the two
  ports agreeing would miss" that the corpus spec calls for.
- **Status:** replicate guard; verification leans on independent invariants, not on it.

---

### Watchlist (candidate findings to confirm in Phase 3/4, not yet actionable)
- `sqrt_whittaker` real-axis tie handling (`>=` / `<=`) — a mode exactly on the real or
  imaginary axis could pick a different sign than a naive numpy port; unit-tested on a
  grid straddling both axes.
- Cross-engine (MATLAB vs Octave) LAPACK divergence on `eig`-heavy chains — expected,
  handled by the pre-declared cross-engine tolerance tier (VERIFICATION_NOTES §10.3),
  logged as an engine finding, never a port fix.
