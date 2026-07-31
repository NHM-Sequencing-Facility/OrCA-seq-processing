#!/bin/bash
#SBATCH --job-name=AS_array
#SBATCH --mem=2G
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --output=%x_%j.log
#SBATCH --array=1-96

set -euo pipefail  # Exit on error, undefined variables, and pipe failures

# Usage function
usage() {
    cat << EOF
Usage: sbatch $0 -input <directory> [options]

Required:
  -input <dir>     Input directory (demuxed SP27 directory)

Optional:
  -min <int>       Minimum amplicon size
  -max <int>       Maximum amplicon size
  -prefix <name>   Output folder prefix/amplicon name (default: 'amplicon')
  -help            Show this help message

Examples:
  # Minimal - no size filtering
  sbatch $0 -input /path/to/demuxed/SP27

  # With output prefix only
  sbatch $0 -input /path/to/demuxed/SP27 -prefix rRNAs

  # With minimum size only
  sbatch $0 -input /path/to/demuxed/SP27 -min 500 -prefix rRNAs

  # With maximum size only
  sbatch $0 -input /path/to/demuxed/SP27 -max 2000 -prefix rRNAs

  # With both size filters
  sbatch $0 -input /path/to/demuxed/SP27 -min 500 -max 2000 -prefix rRNAs
EOF
    exit 1
}

# Initialize variables with defaults
input_dir=""
minimum_size=""
maximum_size=""
outfolder_prefix="amplicon"

# Manual argument parsing to support multi-character flags
# Claude Sonnet v4.5 was used to help wrap argument parsing in shell for this script
while [ $# -gt 0 ]; do
    case "$1" in
        -input)
            if [ -z "${2:-}" ]; then
                echo "Error: -input requires a directory argument"
                usage
            fi
            input_dir="$2"
            shift 2
            ;;
        -min)
            if [ -z "${2:-}" ]; then
                echo "Error: -min requires an integer argument"
                usage
            fi
            minimum_size="$2"
            shift 2
            ;;
        -max)
            if [ -z "${2:-}" ]; then
                echo "Error: -max requires an integer argument"
                usage
            fi
            maximum_size="$2"
            shift 2
            ;;
        -prefix)
            if [ -z "${2:-}" ]; then
                echo "Error: -prefix requires a name argument"
                usage
            fi
            outfolder_prefix="$2"
            shift 2
            ;;
        -help)
            usage
            ;;
        *)
            echo "Error: Unknown option '$1'"
            usage
            ;;
    esac
done

# Validate required arguments
if [ -z "$input_dir" ]; then
    echo "Error: Input directory (-input) is required"
    usage
fi

if [ ! -d "$input_dir" ]; then
    echo "Error: Input directory does not exist: $input_dir"
    exit 1
fi

# Setup directories
basedir=$(dirname $(dirname "$input_dir"))  # Go up two levels from SP27
output_base_dir="${basedir}/03_amplicon_sorted_user"
mkdir -p "${output_base_dir}"

echo "Input directory: $input_dir"
echo "Base directory: $basedir"
echo "Output base directory: $output_base_dir"
echo "Minimum size: ${minimum_size:-'not specified'}"
echo "Maximum size: ${maximum_size:-'not specified'}"
echo "Output folder prefix: $outfolder_prefix"

# Get the file for current array task using mapfile for safety
mapfile -t files < <(find "$input_dir" -maxdepth 1 -name "*.fastq.gz" -type f | sort)

# Check if we have files
if [ ${#files[@]} -eq 0 ]; then
    echo "Error: No fastq.gz files found in $input_dir"
    exit 1
fi

# Check if array task ID is within range
if [ "$SLURM_ARRAY_TASK_ID" -gt "${#files[@]}" ]; then
    echo "Error: Array task ID $SLURM_ARRAY_TASK_ID exceeds number of files (${#files[@]})"
    exit 1
fi

# Get file for this task (arrays are 1-indexed, bash arrays are 0-indexed)
file="${files[$((SLURM_ARRAY_TASK_ID - 1))]}"
echo "Processing file: $file"

# Verify file exists
if [ ! -f "$file" ]; then
    echo "Error: File not found: $file"
    exit 1
fi

# Extract identifier by removing all possible extensions
filename=$(basename "$file")
identifier_nopass="${filename%.fastq.gz}"
identifier_nopass="${identifier_nopass%_pass}"

# Extract dataset name by removing SP27_XXX_SP5_XXX_ pattern
dataset=$(echo "$identifier_nopass" | sed 's/^.*SP27_[0-9]\+_SP5_[0-9]\+_//')

echo "Identifier: $identifier_nopass"
echo "Dataset: $dataset"

# Create output folder
outfolder_rRNAs="${output_base_dir}/${identifier_nopass}/${outfolder_prefix}"
mkdir -p "${outfolder_rRNAs}"

# Build amplicon_sorter command with optional size parameters
source $(conda info --base)/etc/profile.d/conda.sh
conda activate orca-seq

# Construct the command with conditional size parameters
as_cmd="python3 /hpc/scratch/DP/ont_pipeline_testing/amplicon_sorter.py \
 -i ${file} \
 -o ${outfolder_rRNAs} \
 -ar \
 -np ${SLURM_CPUS_PER_TASK}"

# Add minimum size if provided
if [ -n "$minimum_size" ]; then
    as_cmd="$as_cmd -min ${minimum_size}"
fi

# Add maximum size if provided
if [ -n "$maximum_size" ]; then
    as_cmd="$as_cmd -max ${maximum_size}"
fi

echo "Running: $as_cmd"
eval "$as_cmd"

# Verify expected output exists
consensus_file="${outfolder_rRNAs}/consensusfile.fasta"
if [ ! -f "${consensus_file}" ]; then
    echo "Error: Output file '${consensus_file}' not created!"
    echo "Check ${outfolder_rRNAs}/results.txt for details"
    exit 2
fi

echo "Consensus file created successfully"

# Process consensus file headers
conda deactivate && conda activate orca-seq

echo "Reformatting consensus headers..."

# Step 1: Replace parentheses with _readcount_
seqkit replace -p '\((\d+)\)$' -r '_readcount_$1' "${consensus_file}" > "${consensus_file}.tmp1"

# Step 2: Replace _X_Y patterns with sequential _groupN
# Claude Sonnet v4.5 was used to help debug the awk lines for this script
awk '
BEGIN {counter = 1}
/^>/ {
    # Match _X_Y_readcount pattern and replace with _groupN_readcount
    if (match($0, /_[0-9]+_[0-9]+_readcount/)) {
        sub(/_[0-9]+_[0-9]+_readcount/, "_group" counter "_readcount")
        counter++
    }
}
{print}
' "${consensus_file}.tmp1" > "${outfolder_rRNAs}/${identifier_nopass}_consensus_${outfolder_prefix}.fasta"

# Clean up temporary files
rm "${consensus_file}.tmp1" "${consensus_file}"

echo "Processing complete!"
echo "Output: ${outfolder_rRNAs}/${identifier_nopass}_consensus_${outfolder_prefix}.fasta"
