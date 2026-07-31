#!/bin/bash
#SBATCH --job-name=pychopper
#SBATCH --mem=8G
#SBATCH --cpus-per-task=16
#SBATCH --output=%x_%j.log

# Check if input file argument is provided
if [ $# -lt 1 ]; then
    echo "Error: No input file provided"
    echo "Usage: sbatch $0 <input_fastq_file> [output_dir]"
    exit 1
fi

# variables
Qscore=10  # Set the quality score threshold
THREADS=16
infile="$1"  # Input FASTQ file from command line argument

# Extract directory and filename components
indir=$(dirname "$infile")
filename=$(basename "$infile")

# Remove .gz extension if present, then remove .fastq/.fq extension
basename_no_ext="${filename%.gz}"
basename_no_ext="${basename_no_ext%.fastq}"
basename_no_ext="${basename_no_ext%.fq}"

# Set output directory: use user-specified or default to input dir
if [ -n "$2" ]; then
    outdir="$2"
else
    outdir="$indir/pychopped"
fi

outfile="$outdir/pychopped_${basename_no_ext}.fastq.gz"

# Get the directory of this script
SCRIPT_DIR="/hpc/scratch/DP/OrCA-seq-processing/scripts"

# Paths for primer sequences and config
primer_seqs="${SCRIPT_DIR}/../adapters_primers/M13_seqs_for_pychopper.fa"
config="${SCRIPT_DIR}/../adapters_primers/M13_config_for_pychopper.txt"

# Make directories if not existing already
mkdir -p "$outdir"

# for my specific cluster, I need source activate my environments. Most clusters use conda activate (or mamba activate)
source $(conda info --base)/etc/profile.d/conda.sh

conda activate pychopper

echo "Processing: $infile"
echo "Output directory: $outdir"
echo "Output file: $outfile"

pychopper \
 -b "$primer_seqs" \
 -c "$config" \
 -k LSK114 \
 -Q "$Qscore" \
 -w "$outdir"/${basename_no_ext}_rescued.fastq \
 -u "$outdir"/${basename_no_ext}_unclass.fastq \
 -l "$outdir"/${basename_no_ext}_short.fastq \
 -S "$outdir"/${basename_no_ext}_stats.out \
 -p \
 -t "${THREADS}" \
 -m edlib \
 "$infile" | gzip > "$outfile"