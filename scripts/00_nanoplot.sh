#!/bin/bash
#SBATCH --job-name=nanoplot
#SBATCH --mem=4G
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.log

# Check if input file argument is provided
if [ $# -lt 1 ]; then
    echo "Error: No input file provided"
    echo "Usage: sbatch $0 <input_fastq_file> [output_dir]"
    exit 1
fi

# Variable assignment
infile="$1"

# Check if input file exists and is readable
if [ ! -f "$infile" ]; then
    echo "Error: Input file '$infile' not found"
    exit 1
fi

if [ ! -r "$infile" ]; then
    echo "Error: Input file '$infile' is not readable"
    exit 1
fi

# Get directory and filename
indir=$(dirname "$infile")
filename=$(basename "$infile" .gz)

# Set output directory: use user-specified or default to input dir
if [ -n "$2" ]; then
    outdir="$2"
else
    outdir="${indir}/${filename}_nanoplot"
fi

echo "Input file: $infile"
echo "Output directory: $outdir"

# Create output directory if it doesn't exist
mkdir -p "$outdir" || { echo "Error: Failed to create output directory '$outdir'"; exit 1; }

# Load required modules
conda activate orca-seq

# Run NanoPlot with specified parameters
echo "Starting NanoPlot analysis..."
NanoPlot \
    -t 24 \
    --huge \
    --info_in_report \
    --tsv_stats \
    --format pdf \
    --N50 \
    --fastq "$infile" \
    --outdir "$outdir"

# Check if NanoPlot completed successfully
if [ $? -eq 0 ]; then
    echo "NanoPlot completed successfully"
    echo "Results saved in: $outdir"
else
    echo "Error: NanoPlot failed"
    exit 1
fi