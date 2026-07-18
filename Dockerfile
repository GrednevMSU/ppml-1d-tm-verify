# Reproduce report.html WITHOUT MATLAB/Octave, from the committed golden reference.
# The image installs only Python + deps; reference_outputs/matlab/ is baked in, so
# `compare.py --reference matlab` regenerates the report deterministically.
#
#   docker build -t ppml-verify .
#   docker run --rm -v "$PWD/out:/harness/out" ppml-verify
#   # -> out/report.html
FROM python:3.12-slim

WORKDIR /harness
COPY requirements.txt pyproject.toml VERSION ./
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

# Build the report from the committed MATLAB golden; copy it to a mounted volume.
CMD ["sh", "-c", "python verify/compare.py --reference matlab && mkdir -p out && cp report.html verification_report.md out/ && echo 'report written to out/'"]
