# Methodology

How this harness decides that the Python port reproduces the original MATLAB, and why you
can trust the verdict without trusting us. The authoritative, exhaustive version lives in
[`../VERIFICATION_NOTES.md`](../VERIFICATION_NOTES.md); this is the orientation map.

## 1. The original is the source of truth
`matlab_src/` is a verbatim, never-edited copy of the PPML v3.0 1-D TM code. Everything
else is measured against it. Bugs in the original are reproduced bug-for-bug and
documented ([`../ORIGINAL_CODE_FINDINGS.md`](../ORIGINAL_CODE_FINDINGS.md)); any fix is
behind an explicit flag (default OFF) and shows up as `DEVIATION` in the report, never as a
silent `PASS`/`FAIL`.

## 2. A frozen, audit-locked corpus
249 input fixtures (`fixtures/*.mat`) are generated once from a disclosed seed and frozen;
each is SHA-256'd into `fixtures/manifest.json`. Every change to the corpus is a logged
event ([`../CORPUS_CHANGELOG.md`](../CORPUS_CHANGELOG.md)) with a bumped `corpus_version`.
The corpus spans the documented usage scenarios plus edge/branch cases and a randomized
stress sweep; a coverage matrix shows scenario → fixtures.

## 3. Outcome model — outcome before numbers
Every run is a first-class value: `ok` / `error(id)` / `timeout`, with warnings captured.
Two runs are compared on *outcome* first (same error id → PASS; one raises, one doesn't →
finding), and only then on *values*. This is what caught that the original's energy guard
never fires (finding F3) and that an exact diffraction cutoff is a three-way-indeterminate
point across MATLAB / Octave / Python.

## 4. Pass logic and tolerances from a budget
An output passes iff `rel_err ≤ tol_rel` **OR** `abs_err ≤ tol_abs`, per element. The OR is
deliberate: near-zero outputs pass on the abs branch where relative error is
divide-by-noise. Tolerances are **derived from an error budget** (`npw²·L·κ_V·ε`,
evanescent dynamic range for fields) and only then compared to the observed residual —
never fitted a hair above the worst case. Fields get a looser, budget-justified tier
because intra-layer field reconstruction spans `exp(Im(q)·d)`.

## 5. Three comparison axes
- **Gate:** Python vs MATLAB (reference of record). This is what PASS/FAIL is about.
- **Cross-engine:** MATLAB vs Octave (different LAPACK) — a pre-declared INFO tier; a
  residual above it is an engine finding, never a reason to loosen the Python gate.
- **Expect:** each out-of-domain fixture's actual outcome vs its declared expectation.

## 6. Independent triangulation (non-gating INFO)
- **Zero-contrast limit** (`verify/triangulate.py`): a grating with equal inclusion and
  host permittivity is a uniform multilayer; the RCWA result must reduce to an independent
  thin-film TMM. It does, to ~1e-16.
- **Convergence study** (`verify/convergence.py`): physical outputs vs truncation order N,
  with a Richardson/Aitken limit and a discretization-error estimate. Characterises the
  physics of N-convergence (common to both implementations); equivalence is always gated
  at *identical* N.

## 7. Reproduce it yourself
`./verify.sh` (or `verify.bat`) regenerates the reference on your machine and rebuilds the
report; the committed golden + `Dockerfile` let you rebuild the report with no MATLAB at
all. CI runs the whole pipeline on Octave + a pytest matrix + ruff + the notebook.
