# Changelog

All notable changes to the PPML 1-D TM verification harness. Versioning tracks the
harness; the frozen test corpus has its own version in `CORPUS_CHANGELOG.md`
(currently corpus v3), stamped into `fixtures/manifest.json`.

## [3.0.0] — 2026-07-18

First public release. Verified MATLAB→Python migration of the PPML v3.0 1-D TM function
group (`RTA_1d_tm`, `SM_1d_tm`, `field_1d_tm` + kernel), delivered as a differential test
harness. **Result: PASSED** on MATLAB R2025a and GNU Octave 11.3.0.

### Added
- Six-phase pipeline: code analysis (`VERIFICATION_NOTES.md`), frozen 249-fixture corpus
  (SHA-256 locked, `CORPUS_CHANGELOG.md`), dual-engine reference runner (`run_reference.m`,
  MATLAB + Octave), NumPy/SciPy port (`python_src/ppml_1d_tm/`), budget-derived comparison
  engine (`verify/compare.py`, `verify/tolerances.yaml`), and one-command client
  verification (`verify.sh` / `verify.bat`).
- Outcome model: every run is result/error/timeout, compared before numbers; `lastwarn`
  captured as part of the outcome.
- Independent physical invariants (realness `|Im(RR)|≈0`, lossless energy bounds) that do
  NOT rely on the original's vacuous energy guard.
- Convergence study vs truncation order N (`verify/convergence.py`, Richardson-style
  discretization-error estimate) — non-gating.
- Triangulation via the zero-contrast limit vs an independent thin-film TMM
  (`verify/triangulate.py`) — non-gating INFO.
- Standalone figures (`verify/plots.py`): reference-vs-python scatter, error histogram,
  error heatmap. Self-contained `report.html` + `verification_report.md`.
- Honest benchmark (`bench/`, `BENCHMARK.md`) — speed is measured, not claimed.
- pytest suite (250 cases) against the committed MATLAB golden; ruff lint; quickstart
  notebook; 4-job CI (verify / pytest matrix / lint / notebook); Dockerfile for
  MATLAB-free reproduction.

### Documented original-code findings (replicated bug-for-bug)
- **F1** — vacuum impedance truncated to `376.730`.
- **F2** — complex `epssup` silently returns nonphysical results in all three functions.
- **F3** — the RTA energy "guard" is a telescoping identity; it validates nothing on any
  finite input and fails open on NaN.

### Known deviations (flags, default OFF)
- `on_singular='raise'` — at an exact Wood cutoff the port raises where the two reference
  engines themselves disagree (MATLAB NaN vs Octave finite garbage) — DEVIATION, non-gating.
- `fix_vacuum_impedance`, `guard_real_epssup` — optional fixes to F1/F2, DEVIATION when ON.

### Scope
1-D TM group only. Biaxial `epar_1d` and 2-D groups are out of scope.
