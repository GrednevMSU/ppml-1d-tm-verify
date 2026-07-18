@echo off
REM One command (Windows): regenerate references (auto-detect MATLAB/Octave) -> run the
REM Python port -> build report.html -> open it. PPML 1-D TM differential test harness.
setlocal enabledelayedexpansion
cd /d "%~dp0"
echo == PPML 1-D TM differential test harness ==

REM ---- 1. Python venv + deps ----
if not exist ".venv\Scripts\python.exe" (
  echo [1/5] creating Python venv...
  python -m venv .venv
)
echo [1/5] installing Python deps...
".venv\Scripts\python.exe" -m pip install -q --disable-pip-version-check -r requirements.txt

REM ---- 2. frozen corpus present? ----
if not exist "fixtures\manifest.json" (
  echo ERROR: fixtures\manifest.json missing. The frozen corpus must be committed.
  exit /b 3
)

REM ---- 3. detect engines ----
set "MATLAB="
set "OCTAVE="
where matlab >nul 2>&1 && set "MATLAB=matlab"
where octave >nul 2>&1 && set "OCTAVE=octave"
if "%MATLAB%"=="" if "%OCTAVE%"=="" (
  echo ERROR: need MATLAB or Octave to generate reference outputs.
  echo        Install GNU Octave ^(free^): https://octave.org  then re-run.
  exit /b 2
)
echo [2/5] engines: MATLAB='%MATLAB%'  Octave='%OCTAVE%'

REM ---- 4. references over FROZEN fixtures ----
if not "%MATLAB%"=="" (
  echo [3/5] MATLAB reference...
  matlab -sd "%CD%" -batch "run_reference" || (echo MATLAB run_reference failed & exit /b 4)
)
if not "%OCTAVE%"=="" (
  echo [3/5] Octave reference...
  octave -q run_reference.m || (echo Octave run_reference failed & exit /b 4)
)

REM ---- 5. reference of record ----
set "REF=octave"
if not "%MATLAB%"=="" set "REF=matlab"
echo [4/5] comparing Python port vs reference (%REF%)...
".venv\Scripts\python.exe" verify\compare.py --reference %REF%
set RC=%ERRORLEVEL%

REM ---- 6. open report ----
echo [5/5] report: %CD%\report.html
start "" "%CD%\report.html"

if "%RC%"=="0" (echo RESULT: PASSED) else (echo RESULT: FAILED ^(exit %RC%^))
exit /b %RC%
