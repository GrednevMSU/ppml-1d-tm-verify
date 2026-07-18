"""PPML v3.0 — 1-D TM group, Python port (NumPy/SciPy only).

Public API mirrors matlab_src/1d_tm/. See NAME_MAP.md for renames and
ORIGINAL_CODE_FINDINGS.md for replicated original behaviors (F1/F2/F3).
"""
from .general import smpropag_bw_cond, smpropag_fw_cond, sqrt_whittaker
from .tm_1d import PPMLEnergyError, RTA_1d_tm, SM_1d_tm, field_1d_tm

__all__ = [
    "RTA_1d_tm", "SM_1d_tm", "field_1d_tm", "PPMLEnergyError",
    "sqrt_whittaker", "smpropag_fw_cond", "smpropag_bw_cond",
]
