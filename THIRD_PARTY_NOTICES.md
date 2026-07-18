# THIRD_PARTY_NOTICES

This repository contains, verbatim and unmodified, third-party source code. Its
provenance and license are recorded here so downstream users know exactly what they are
redistributing.

---

## PPML v3.0 (original MATLAB code)

- **Location in this repo:** `matlab_src/` (audit copy of the 1-D TM function group and
  its shared kernel; the full upstream package is at `../PPML_v3/`).
- **Author:** Simone Zanotto, with contributions by Gaia Da Prato.
- **Copyright:** © 2016–2022 Simone Zanotto. All rights reserved.
- **License:** BSD 3-Clause (see `LICENSE`, copied verbatim from the upstream
  `license.txt`).
- **What it is:** an RCWA / Fourier Modal Method implementation in the scattering-matrix
  formalism (Whittaker–Culshaw / Liscidini et al.), with Li's correct factorization rules
  for TM polarization.
- **Files included (1-D TM scope):**
  - `matlab_src/1d_tm/RTA_1d_tm.m`
  - `matlab_src/1d_tm/SM_1d_tm.m`
  - `matlab_src/1d_tm/field_1d_tm.m`
  - `matlab_src/general/sqrt_whittaker.m`
  - `matlab_src/general/smpropag_fw_cond.m`
  - `matlab_src/general/smpropag_bw_cond.m`
- **Modifications:** NONE. The MATLAB files are the source of truth and are never edited;
  the verification harness only reads and runs them.

### Academic citation
The upstream authors request that results obtained through this code, appearing in an
academic publication, cite the original PPML source (and, where applicable, the papers
referenced in the individual `examples/*.m` headers, e.g. APL 103, 091110 (2013);
Sci. Rep. 6 (2016)).

---

## This verification harness (derived work)

The Python port (`python_src/`), test corpus, comparison engine, figures, and
client-experience tooling in this repository are a derived work built around the original
PPML code. They are distributed under the **same BSD 3-Clause terms** as the original (see
`LICENSE`), and do **not** alter the original algorithm. The port replicates two
documented bugs of the original bug-for-bug (see `ORIGINAL_CODE_FINDINGS.md`, F1/F3) —
these are behaviors of the ORIGINAL code, faithfully reproduced, not defects introduced by
the port.

## Python dependencies
NumPy, SciPy, PyYAML, matplotlib (and, for the dev/report extras, pytest, ruff, jupyter)
— each under its own permissive license (BSD/MIT/PSF-style). See each package's
distribution for the exact text. Pinned versions in `requirements.txt` / `pyproject.toml`.
