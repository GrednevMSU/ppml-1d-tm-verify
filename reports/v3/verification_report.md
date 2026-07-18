# PPML 1-D TM — Differential Test Harness — PASSED

## Executive summary
- The Python port of the PPML 1-D TM group (RTA / SM / field) was checked against the original MATLAB code on 249 fixtures spanning subwavelength-to-multi-order periods, 0–80° incidence, dielectric and metallic (lossy) layers, and conducting interfaces.
- Reflectance/transmittance and S-matrix coefficients agree with MATLAB to about 8e-12 absolute — at the level predicted by the accumulated round-off budget (§10.4: npw²·L·κ_V·ε, envelope 5.3e-11), ~6× within it, NOT machine-epsilon (2e-16); sampled fields to about 3e-10 relative (looser by physics — thick-metal evanescent dynamic range, not a port error; the port matches MATLAB better than the two reference engines match each other).
- Overall gating result: PASSED with 0 gating failures.
- Two behaviours in the ORIGINAL code are documented and replicated bug-for-bug: a truncated vacuum impedance 376.730 (F1), and an energy 'guard' whose checked quantity RR+TT+ΣAA is identically 1 by telescoping (F3) — so it validates nothing on any finite input, and additionally fails open on NaN and on the documented-forbidden complex-permittivity input (F2).
- Scope: 1-D TM only; exact diffraction-cutoff points are numerically indeterminate (all three implementations disagree there) and reported non-gating.

## Pass logic
`PASS iff OR(rel_err <= tol_rel, abs_err <= tol_abs)`, per element. Statuses: PASS(rel)/PASS(abs)/PASS(both)/FAIL/INFO/DEVIATION(flag=…).

## Summary (Python vs MATLAB)
| function | status | fixtures | gating FAIL | max rel | max abs | tol_rel/tol_abs |
|---|---|---|---|---|---|---|
| RTA | PASS(abs) | 164 | 0 | 1.81e-10 | 8.32e-12 | 1e-11/1e-10 |
| SM | PASS(both) | 53 | 0 | 4.74e-12 | 1.90e-13 | 1e-11/1e-10 |
| field | PASS(rel) | 32 | 0 | 2.89e-10 | 4.55e-10 | 1e-09/1e-11 |

Status qualifies every row (OR pass-logic): `field` prints max abs above its tol_abs but passes as `PASS(rel)` — tol_abs is only a floor for near-null field points. Worst-case fixtures reported SEPARATELY for rel and abs:
- worst rel = 2.89e-10 @ `rand_085`
- worst abs = 4.55e-10 @ `rand_085`

## Tiers — observed vs budget vs threshold
| tier | applies | observed (max) | budget / basis | threshold | gating |
|---|---|---|---|---|---|
| gate:RTA/SM | R/T, S | abs 8.32e-12 | npw²·L·κ_V·ε ≈ 5.3e-11 | abs 1e-10 / rel 1e-11 | YES |
| gate:field | Ex,Ez,Sz | rel 2.89e-10 | evanescent exp(Im q·d), §10.5 | rel 1e-9 / abs 1e-11 | YES |
| cross-engine | MATLAB↔Octave | abs 6.31e-10 | LAPACK divergence, §10.3 | abs 5e-09 | INFO |
| near_anomaly | 9 fixtures | (hypersensitive) | Wood cutoff physics, §10.1 | rel 1e-3 | INFO |
| resonance | 0 fixtures | cond max 292 | redheffer_cond, §10.2 | >1e8 | INFO |

## Findings (FAIL / DEVIATION / notable INFO)
- No gating FAIL; 1 documented DEVIATION.
- `branch_wood_anomaly_+0.00` [gate] **DEVIATION(flag=on_singular)** — matlab: ok/NaN+MATLAB:illConditionedMatrix · octave: ok/finite(RR=0.769)+Octave:singular-matrix · Python: raise(numpy:LinAlgError, flag=on_singular) — three-way indeterminate at exact cutoff (D2)
- INFO (non-gating): 9 `near_anomaly` fixtures on the loose Wood tier (rel 1e-3): `branch_wood_anomaly_+0.00`, `branch_wood_anomaly_+0.15`, `rand_012`, `rand_031`, `rand_059`, `rand_071` ….

## Scenario coverage matrix
| category | fixtures |
|---|---|
| scenario | 6 |
| branch | 15 |
| conducting_interface | 5 |
| out_of_domain | 3 |
| randomized | 220 |
| **TOTAL** | **249** |

## Cross-engine (matlab vs octave, INFO)
- RTA: max abs diff 5.04e-12 (tier 5e-09)
- SM: max abs diff 2.83e-13 (tier 5e-09)
- field: max abs diff 6.31e-10 (tier 5e-09)

## Scope & limitations
- IN: `RTA_1d_tm`, `SM_1d_tm`, `field_1d_tm` + `sqrt_whittaker`, `smpropag_{fw,bw}_cond`.
- OUT: `epar_1d` (biaxial 1-D), `2d_rect`, `2d_Lshape`, `fill_F`, `fill_new_F`, `ZSM_2d_*` — declared out of scope; require toolbox `sinc` shim if later added.
- Reference on this host: MATLAB (record) + Octave 11.3.0 (cross-check). MATLAB-vs-Octave discrepancies are gated at a **pre-declared cross-engine tier** (§10.3, `rel 5e-9`) and above it become findings — never a silent loosening of the Python gate.
- `epssup` complex is out-of-domain (docstring); handled via the outcome model (§9a) with per-function expectations, and logged as finding F2 (guard asymmetry).
- A fixture run is result/error/timeout (§9a); guard-straddling near cutoffs is expected and classified, not a crash.
- At an EXACT Rayleigh/Wood cutoff the linear system is singular and all three impls diverge (MATLAB `RR=NaN`, Octave finite garbage, Python raises by default) — a fundamentally indeterminate point (D2). Python-vs-MATLAB status is `DEVIATION(flag=on_singular)`; cross-engine status is INFO (the references disagree). Not a port defect. Real usage does not sit exactly on a cutoff.

## Environment
- Reference (matlab): MATLAB 25.1 / LAPACK NAG Performance Components (NPC) Release 1.2.1, supporting Linear Algebra PACKage (LAPACK 3.9.1)
- Cross-engine (octave): Octave 11.3.0 / LAPACK Linear Algebra PACKage Version 3.12.0
- Python 3.14.6 / NumPy 2.5.1 / BLAS accelerate
- corpus v3, seed 20260717, schema 1