#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
STUB_DIR="${REPO_DIR}/tests/stub"
REF="${STUB_DIR}/reference.fa"
KNOWN_SITES="${STUB_DIR}/known_sites.vcf.gz"
GERMLINE="${STUB_DIR}/known_sites.vcf.gz"
SAMPLE_SHEET="${STUB_DIR}/samples.csv"
RUN_ROOT="${STUB_DIR}/results"
WORK_ROOT="${STUB_DIR}/work"

if command -v nextflow >/dev/null 2>&1; then
  NEXTFLOW_BIN="$(command -v nextflow)"
elif [[ -x "/Users/adigenova/bin/nextflow" ]]; then
  NEXTFLOW_BIN="/Users/adigenova/bin/nextflow"
elif [[ -x "/Users/adigenova/nextflow_pmd/nextflow" ]]; then
  NEXTFLOW_BIN="/Users/adigenova/nextflow_pmd/nextflow"
else
  echo "ERROR: nextflow not found in PATH and fallback binary is missing." >&2
  exit 1
fi

run_case() {
  local mode="$1"
  shift
  local outdir="${RUN_ROOT}/${mode}"
  local workdir="${WORK_ROOT}/${mode}"

  echo "=== Running stub mode: ${mode} ==="
  "${NEXTFLOW_BIN}" run "${REPO_DIR}/mutect2_pipeline.nf" \
    -profile local,stub \
    -stub-run \
    -ansi-log false \
    -w "${workdir}" \
    --reference "${REF}" \
    --known_sites "${KNOWN_SITES}" \
    --germline_resource "${GERMLINE}" \
    --outdir "${outdir}" \
    "$@"
}

# 1) Sample sheet, no PoN (paired + tumor-only paths)
run_case "sample_sheet" \
  --sample_sheet "${SAMPLE_SHEET}"

# 2) Build PoN from normal CRAM(s)
run_case "sample_sheet_build_pon" \
  --sample_sheet "${SAMPLE_SHEET}" \
  --pon_crams "${STUB_DIR}/data/NOR1.cram"

# 3) Reuse prebuilt PoN in calling mode
run_case "sample_sheet_with_pon" \
  --sample_sheet "${SAMPLE_SHEET}" \
  --panel_of_normals "${RUN_ROOT}/sample_sheet_build_pon/pon/panel_of_normals.vcf.gz"

# 4) Reads mode (no sample sheet)
run_case "reads_mode" \
  --reads "${STUB_DIR}/data/*_{R1,R2}.fastq.gz" \
  --tumor_sample TUM3

# 5) CRAM mode (no sample sheet)
run_case "cram_mode" \
  --crams "${STUB_DIR}/data/TUM2.cram" \
  --tumor_sample TUM2

echo "All stub modes completed successfully."
