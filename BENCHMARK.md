# Benchmark

**Speed is not a goal of this project.** The deliverable is numerical *equivalence* to the
original MATLAB, not a faster code. This page exists so that any performance statement is
*measured*, not asserted — and so you can reproduce it.

Reproduce:
```bash
python bench/benchmark.py                                   # Python port
matlab -batch "addpath('bench'); bench_matlab"              # original, in your MATLAB
```

## Measured (single machine — indicative only)

- Structure: 3-layer GaAs / Au-slit / GaAs grating (UniversalLineshapes), θ = 0.1°.
- Metric: wall-clock per `RTA_1d_tm` call, median of 200 (N ≤ 41) / 80 (N = 61) reps,
  after a warm-up call.
- Machine: Apple Silicon (macOS, arm64), Python 3.14 + NumPy 2.5.1 (Accelerate BLAS);
  MATLAB R2025a. **Your numbers will differ** — hardware, BLAS, MATLAB JIT all matter.

| N (=2·halfnpw+1) | Python median [ms] | MATLAB median [ms] | ratio (MATLAB/Python) |
|---:|---:|---:|---:|
| 11 |  0.72 |  1.83 | 2.5× |
| 21 |  1.41 |  3.21 | 2.3× |
| 41 |  5.50 | 11.45 | 2.1× |
| 61 | 12.19 | 31.10 | 2.6× |

On this machine the Python port runs about **2–2.5× faster** than the original MATLAB for
this structure. Take that with salt:

- One structure, one machine, one BLAS. Not a portfolio-wide claim.
- MATLAB's first call carries JIT/parse overhead (its *mean* at N=11 was ~8.8 ms vs a
  1.8 ms median — startup, not steady state); the medians above are the fair comparison.
- Neither side is performance-tuned. The port mirrors the MATLAB structure faithfully
  (correctness over style), and the MATLAB is the original, unmodified.
- Complexity is the same on both sides: the per-call cost is dominated by an `eig` on an
  `npw × npw` matrix per layer, so both scale ~`O(L · npw³)`.

**Bottom line:** the port is not slower, and is somewhat faster here — but if you need a
speed guarantee for a specific workload, measure it on your own hardware with the two
scripts above.
