#!/bin/bash
#SBATCH --job-name=demux_and_clean
#SBATCH --mem=8G
#SBATCH --cpus-per-task=16
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.log


set -euo pipefail  # Exit on error, undefined variables, and pipe failures

# Check if input file argument is provided
if [ $# -eq 0 ]; then
    echo "Error: No input file provided"
    echo "Usage: sbatch $0 <pychopped_input_file>"
    echo "Example: sbatch $0 /path/to/pychopped/pychopped_sample1.gz"
    exit 1
fi

# Define variables
infile="$1"
threads=16
e_rate=0.1

# Extract directory and filename components
indir=$(dirname "$infile")
filename=$(basename "$infile")

# Extract dataset name by removing "pychopped_" prefix and extensions
dataset="${filename#pychopped_}"
dataset="${dataset%.fastq.gz}"
dataset="${dataset%.fastq}"
dataset="${dataset%.fq.gz}"
dataset="${dataset%.fq}"
dataset="${dataset%.gz}"
dataset="${dataset%_pass}"

# Setup directories
parent_dir=$(dirname "$indir")
outdir="$parent_dir/02_demuxed"
mkdir -p "$outdir"/{SP5,SP27}

# Get the directory of this script
SCRIPT_DIR="/hpc/scratch/DP/OrCA-seq-processing/scripts"

# Define adapter and primer paths
adapters_SP5="${SCRIPT_DIR}/../adapters_primers/M13_amplicon_indices_forward.fa"
adapters_SP27="${SCRIPT_DIR}/../adapters_primers/M13_amplicon_indices_reverse_rc.fa"
echo "========================================="
echo "Processing: $infile"
echo "Dataset name: $dataset"
echo "Output directory: $outdir"
echo "========================================="

# Verify required files exist
for f in "$infile" "$adapters_SP5" "$adapters_SP27"; do
    if [ ! -f "$f" ]; then
        echo "Error: Required file not found: $f"
        exit 1
    fi
done

# Activate conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate orca-seq

# First round: demultiplex with SP5 adapters
echo "Round 1: Demultiplexing with SP5 adapters..."
cutadapt \
    --action=trim \
    -e "$e_rate" \
    -j "$threads" \
    --rc \
    -g file:"$adapters_SP5" \
    -o "$outdir"/SP5/{name}_"$dataset".fastq.gz \
    "$infile" \
    --json="$outdir"/SP5/cutadapt_SP5_"$dataset".json

# Get list of SP5 identifiers from output files
suffix="_${dataset}.fastq.gz"
mapfile -t identifiers < <(
    find "$outdir/SP5" -name "*$suffix" -type f -exec basename {} \; | 
    sed "s/$suffix$//" | 
    grep -v "unknown"
)

if [ ${#identifiers[@]} -eq 0 ]; then
    echo "Error: No valid identifiers found from SP5 demultiplexing"
    exit 1
fi

echo "Found ${#identifiers[@]} identifiers from SP5 demultiplexing"

# Second round: demultiplex each SP5 file with SP27 adapters
echo "Round 2: Demultiplexing with SP27 adapters..."
for identifier in "${identifiers[@]}"; do
    echo "  Processing: $identifier"
    
    cutadapt \
        --action=trim \
        -e "$e_rate" \
        -j "$threads" \
        --rc \
        -a file:"$adapters_SP27" \
        -o "$outdir"/SP27/{name}_"$identifier"_"$dataset".fastq.gz \
        "$outdir"/SP5/"$identifier"_"$dataset".fastq.gz \
        --json="$outdir"/SP27/"$identifier"_"$dataset".json
done

echo "Demultiplexing complete!"

# Cleaning files that are 'unknown' demuxes and non-existent adapter combinations
# Remove unknown files (i.e., the adapter sequences were unidentifiable, so the read was not able to be assigned)
echo "Removing files with 'unknown' in name..."
removed_unknown=$(find "$outdir" -type f -name "*unknown*" -delete -print | wc -l)
echo "Removed $removed_unknown files with 'unknown'"

# Remove invalid SP27 combinations (the structure of a 96 well plate allows all 12 SP5 adapter and SP27 adapters 1 to 8 only)
removed_invalid=$(find "$outdir" -type f \( \
    -name "*SP27_009*" -o \
    -name "*SP27_010*" -o \
    -name "*SP27_011*" -o \
    -name "*SP27_012*" \) -print -delete | wc -l)
echo "Removed $removed_invalid files with invalid SP27 indices"

echo "Pipeline complete!"
echo "Results in: $outdir"


