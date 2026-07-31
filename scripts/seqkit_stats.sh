#!/bin/bash
#SBATCH --job-name=seqkit_stats
#SBATCH --output=seqkit_stats_%j.out
#SBATCH --error=seqkit_stats_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=day

# =============================================================================
# seqkit_stats.sh
# Run seqkit stats on one or more FASTQ or FASTA files.
#
# Usage:
#   sbatch seqkit_stats.sh <input_path> <output_dir>
#   bash   seqkit_stats.sh <input_path> <output_dir>
#
# Arguments:
#   $1  input_path   Path to a single FASTQ/FASTA file, OR a directory
#                    containing FASTQ and/or FASTA files (searched recursively).
#                    Supported extensions: .fastq, .fastq.gz, .fq, .fq.gz,
#                                         .fasta, .fasta.gz, .fa, .fa.gz
#   $2  output_dir   Directory to write results into (created if absent).
#
# Outputs (written to <output_dir>/):
#   seqkit_stats.tsv          Full seqkit stats table (TSV, all metrics).
#   seqkit_stats_summary.txt  Human-readable run summary.
# =============================================================================

set -euo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate orca-seq

# ---------------------------------------------------------------------------
# 1. Argument validation
# ---------------------------------------------------------------------------
if [[ $# -lt 2 ]]; then
    echo "[ERROR] Usage: $(basename "$0") <input_path> <output_dir>" >&2
    exit 1
fi

INPUT_PATH="${1}"
OUTPUT_DIR="${2}"

if [[ ! -e "${INPUT_PATH}" ]]; then
    echo "[ERROR] Input path does not exist: ${INPUT_PATH}" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# 2. Build file list
# ---------------------------------------------------------------------------
declare -a FASTQ_FILES

if [[ -f "${INPUT_PATH}" ]]; then
    # Single file supplied — validate extension
    if [[ "${INPUT_PATH}" =~ \.(fastq|fastq\.gz|fq|fq\.gz|fasta|fasta\.gz|fa|fa\.gz)$ ]]; then
        FASTQ_FILES=("${INPUT_PATH}")
    else
        echo "[ERROR] File does not appear to be a FASTQ or FASTA: ${INPUT_PATH}" >&2
        exit 1
    fi
elif [[ -d "${INPUT_PATH}" ]]; then
    # Directory supplied — find all FASTQ/FASTA files recursively
    mapfile -t FASTQ_FILES < <(
        find "${INPUT_PATH}" -type f \
            \( -name "*.fastq"    -o -name "*.fastq.gz" \
               -o -name "*.fq"    -o -name "*.fq.gz"    \
               -o -name "*.fasta" -o -name "*.fasta.gz" \
               -o -name "*.fa"    -o -name "*.fa.gz"    \) \
        | sort
    )
else
    echo "[ERROR] Input path is neither a file nor a directory: ${INPUT_PATH}" >&2
    exit 1
fi

if [[ ${#FASTQ_FILES[@]} -eq 0 ]]; then
    echo "[ERROR] No FASTQ or FASTA files found under: ${INPUT_PATH}" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# 3. Prepare output directory (and SLURM log directory)
# ---------------------------------------------------------------------------
mkdir -p "${OUTPUT_DIR}"

OUTPUT_TSV="${OUTPUT_DIR}/seqkit_stats.tsv"
OUTPUT_SUMMARY="${OUTPUT_DIR}/seqkit_stats_summary.txt"

# ---------------------------------------------------------------------------
# 4. Environment loading & validation
# ---------------------------------------------------------------------------
# Verify seqkit is available
if ! command -v seqkit &>/dev/null; then
    echo "[ERROR] seqkit not found on PATH. Load the appropriate module or activate conda env." >&2
    exit 1
fi

SEQKIT_VERSION=$(seqkit version 2>&1 | head -1)
THREADS="${SLURM_CPUS_PER_TASK:-8}"

# ---------------------------------------------------------------------------
# 5. Log run info
# ---------------------------------------------------------------------------
echo "============================================================"
echo "  seqkit stats run"
echo "============================================================"
echo "  Script       : $0"
echo "  Input path   : ${INPUT_PATH}"
echo "  Output dir   : ${OUTPUT_DIR}"
echo "  Files found  : ${#FASTQ_FILES[@]}"
echo "  Threads      : ${THREADS}"
echo "  seqkit       : ${SEQKIT_VERSION}"
echo "  Start time   : $(date '+%Y-%m-%d %H:%M:%S')"
echo "  SLURM job ID : ${SLURM_JOB_ID:-N/A (interactive)}"
echo "============================================================"
echo ""
printf '  %s\n' "${FASTQ_FILES[@]}"
echo ""

# ---------------------------------------------------------------------------
# 6. Run seqkit stats
# ---------------------------------------------------------------------------
# Flags:
#   -a / --all        Include all metrics: GC%, N50, Q20/Q30 for FASTQ, etc.
#   -T / --tabular    TSV output (machine-friendly)
#   -j THREADS        Parallel threads
#   -e / --skip-err   Skip files with read errors rather than aborting
#   -b / --basename   Store only the filename in the 'file' column

echo "[INFO] Running seqkit stats..."

seqkit stats \
    --all \
    --tabular \
    --threads "${THREADS}" \
    --skip-err \
    "${FASTQ_FILES[@]}" \
    > "${OUTPUT_TSV}"

echo "[INFO] Stats written to: ${OUTPUT_TSV}"

# ---------------------------------------------------------------------------
# 7. Write human-readable summary
# ---------------------------------------------------------------------------
{
    echo "seqkit stats summary"
    echo "Run date   : $(date '+%Y-%m-%d %H:%M:%S')"
    echo "Input path : ${INPUT_PATH}"
    echo "seqkit     : ${SEQKIT_VERSION}"
    echo "Files      : ${#FASTQ_FILES[@]}"
    echo ""
    # Pretty-print the TSV with column alignment
    column -t -s $'\t' "${OUTPUT_TSV}"
} > "${OUTPUT_SUMMARY}"

echo "[INFO] Summary written to: ${OUTPUT_SUMMARY}"

# ---------------------------------------------------------------------------
# 8. Quick console preview
# ---------------------------------------------------------------------------
echo ""
echo "--- Preview (first 5 data rows) ---"
head -6 "${OUTPUT_TSV}" | column -t -s $'\t'
echo ""

echo "============================================================"
echo "  Finished : $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================================"