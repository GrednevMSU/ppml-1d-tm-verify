#!/usr/bin/env bash
# One command: regenerate references (auto-detect MATLAB/Octave) -> run the Python port
# -> build report.html -> open it. PPML 1-D TM differential test harness (Phase 6).
set -euo pipefail
cd "$(dirname "$0")"
HERE="$PWD"

echo "== PPML 1-D TM differential test harness =="

# ---- 1. Python environment -------------------------------------------------
PY="$HERE/.venv/bin/python"
if [ ! -x "$PY" ]; then
  echo "[1/5] creating Python venv..."
  python3 -m venv .venv
fi
echo "[1/5] installing Python deps..."
"$HERE/.venv/bin/pip" install -q --disable-pip-version-check -r requirements.txt

# ---- 2. sanity: frozen corpus must be present ------------------------------
if [ ! -f fixtures/manifest.json ]; then
  echo "ERROR: fixtures/manifest.json missing. The frozen corpus must be committed." >&2
  exit 3
fi

# ---- 3. detect engines -----------------------------------------------------
MATLAB=""; OCTAVE=""
if command -v matlab >/dev/null 2>&1; then
  MATLAB="matlab"
else
  for m in /Applications/MATLAB_*.app/bin/matlab "/Applications/MATLAB/"*/bin/matlab; do
    [ -x "$m" ] && MATLAB="$m" && break
  done
fi
command -v octave >/dev/null 2>&1 && OCTAVE="octave"

if [ -z "$MATLAB" ] && [ -z "$OCTAVE" ]; then
  echo "ERROR: need MATLAB or Octave to generate the reference outputs." >&2
  echo "       Install GNU Octave (free): https://octave.org  then re-run." >&2
  exit 2
fi
echo "[2/5] engines: MATLAB='${MATLAB:-none}'  Octave='${OCTAVE:-none}'"

# ---- 4. regenerate reference outputs over the FROZEN fixtures ---------------
if [ -n "$MATLAB" ]; then
  echo "[3/5] MATLAB reference..."
  "$MATLAB" -sd "$HERE" -batch "run_reference" || { echo "MATLAB run_reference failed" >&2; exit 4; }
fi
if [ -n "$OCTAVE" ]; then
  echo "[3/5] Octave reference..."
  "$OCTAVE" -q run_reference.m || { echo "Octave run_reference failed" >&2; exit 4; }
fi

# ---- 5. reference of record: MATLAB if present, else Octave -----------------
REF="octave"; [ -n "$MATLAB" ] && REF="matlab"
echo "[4/5] comparing Python port vs reference ($REF)..."
set +e
"$PY" verify/compare.py --reference "$REF"
RC=$?
set -e

# ---- 6. open the report ----------------------------------------------------
REPORT="$HERE/report.html"
echo "[5/5] report: $REPORT"
case "$(uname)" in
  Darwin) open "$REPORT" 2>/dev/null || true ;;
  Linux)  xdg-open "$REPORT" 2>/dev/null || true ;;
esac

if [ "$RC" -eq 0 ]; then echo "RESULT: PASSED"; else echo "RESULT: FAILED (exit $RC)"; fi
exit "$RC"
