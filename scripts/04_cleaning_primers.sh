#!/bin/bash
#SBATCH --job-name=trim_primers
#SBATCH --mem=2G
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --output=trim_%A_%a.log
#SBATCH --array=1-96

# Usage information
usage() {
    cat << EOF
Usage: sbatch $0 <amplicon_sorted_dir> <amplicon_type> --r1-primers <file> [options]

Required arguments:
  <amplicon_sorted_dir>    Base directory containing sample subdirectories (e.g., Lakes_day1/amplicon_sorted)
  <amplicon_type>          Amplicon type to process (e.g., rRNAs, COIs)
  --r1-primers <file>      FASTA file with Round 1 primers (linked pairs, degenerate bases permitted)

Optional arguments:
  --r2-primers <file>      FASTA file with Round 2 primers (unlinked, degenerate bases permitted)
  --run-round2             Run Round 2 trimming on untrimmed sequences (default: skip Round 2)

Primer FASTA format:
  - Headers must contain "Forward" or "Reverse" to indicate primer orientation
  - For Round 1 (linked): Headers must contain pair identifiers (e.g., "_A", "_B")
    Format: >PrimerName|Description_Forward_A
            >PrimerName|Description_Reverse_A
  - For Round 2 (unlinked): Headers just need "Forward" or "Reverse"
    Format: >PrimerName|Description_Forward
            >PrimerName|Description_Reverse

Examples:
  Round 1 primers (COI_primers_r1.fasta):
    >[name]|Forward_A
    TGAYATTGG
    >[name]|Reverse_A
    TANARYNCA
    >[name]|Forward_B
    WTAYCCNCC
    >[name]|Reverse_B
    TANARAAYCA

Usage:
  # Default (skip Round 2, discard untrimmed sequences):
  sbatch $0 Lakes_day1/amplicon_sorted COIs --r1-primers primers_r1.fasta --cluster-round1
  
  # With Round 2 trimming:
  sbatch $0 Lakes_day1/amplicon_sorted COIs --r1-primers primers_r1.fasta --r2-primers primers_r2.fasta --run-round2

Note: The primers MUST BE IN THE CORRECT ORIENTATION (5' to 3'), so check that you've correctly oriented them (might need to reverse complement the 3' primer hehe).
EOF
    exit 1
}

# Initialize variables
amplicon_sorted_dir=""
amplicon_type=""
r1_primers_file=""
r2_primers_file=""
run_round2=false
cluster_round1=false

# Parse arguments
if [ $# -lt 2 ]; then
    usage
fi

amplicon_sorted_dir="$1"
amplicon_type="$2"
shift 2

# Claude Sonnet v4.5 was used to help wrap argument parsing in shell for this script
while [[ $# -gt 0 ]]; do
    case $1 in
        --r1-primers)
            r1_primers_file="$2"
            shift 2
            ;;
        --r2-primers)
            r2_primers_file="$2"
            shift 2
            ;;
        --run-round2)
            run_round2=true
            shift
            ;;
        --cluster-round1)
            cluster_round1=true
            shift
            ;;
        --help|-h)
            usage
            ;;
        *)
            echo "Error: Unknown argument $1"
            usage
            ;;
    esac
done

# Validate required arguments
if [ -z "$r1_primers_file" ]; then
    echo "Error: Missing required argument --r1-primers"
    usage
fi

if [ ! -f "$r1_primers_file" ]; then
    echo "Error: Round 1 primers file not found: $r1_primers_file"
    exit 1
fi

if [ -n "$r2_primers_file" ] && [ ! -f "$r2_primers_file" ]; then
    echo "Error: Round 2 primers file not found: $r2_primers_file"
    exit 1
fi

if [ ! -d "$amplicon_sorted_dir" ]; then
    echo "Error: Directory $amplicon_sorted_dir does not exist"
    exit 1
fi

# Build array of sample directories
mapfile -t sample_dirs < <(find "$amplicon_sorted_dir" -mindepth 2 -maxdepth 2 -type d -name "$amplicon_type" | sort)

if [ ${#sample_dirs[@]} -eq 0 ]; then
    echo "Error: No $amplicon_type directories found in $amplicon_sorted_dir"
    exit 1
fi

# Get the sample directory for this array task
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    echo "Error: This script must be run as an SLURM array job"
    echo "Example: sbatch --array=1-${#sample_dirs[@]} $0 $amplicon_sorted_dir $amplicon_type --r1-primers $r1_primers_file"
    exit 1
fi

sample_amplicon_dir="${sample_dirs[$((SLURM_ARRAY_TASK_ID - 1))]}"

if [ -z "$sample_amplicon_dir" ]; then
    echo "Error: Invalid array task ID $SLURM_ARRAY_TASK_ID (max: ${#sample_dirs[@]})"
    exit 1
fi


echo "Array Task ID: $SLURM_ARRAY_TASK_ID / ${#sample_dirs[@]}"
echo "Processing: $sample_amplicon_dir"


# Get consensus file
consensus_file=$(find "$sample_amplicon_dir" -type f -name "*_consensus_${amplicon_type}.fasta" | head -n 1)

if [ ! -f "$consensus_file" ]; then
    echo "Error: No consensus file found in $sample_amplicon_dir"
    echo "Expected pattern: *_consensus_${amplicon_type}.fasta"
    exit 1
fi

echo "Input file: $consensus_file"

# Extract sample identifier from directory structure
identifier=$(basename "$(dirname "$sample_amplicon_dir")")
basedir=$(dirname $(dirname "$(dirname "$sample_amplicon_dir")"))
echo "Sample identifier: $identifier"
echo "Amplicon type: $amplicon_type"

# Create output directory
output_dir="${basedir}/04_primerless/${identifier}/${amplicon_type}"
mkdir -p "$output_dir"
echo "Output directory: $output_dir"


# Define output files
primerless_fasta_round1="${output_dir}/round1_amplicon_${identifier}.fasta"
cleanest_fasta_round1="${output_dir}/cleaned_amplicon_${identifier}.fasta"
untrimmed_fasta_round1="${output_dir}/untrimmed_round1_${identifier}.fasta"
primerless_fasta_round2="${output_dir}/wrong_primers_${identifier}.fasta"

# Define temporary files
temp_primers="${output_dir}/temp_primers_${identifier}.fa"
temp_primer_locations="${output_dir}/primer_locations_${identifier}.txt"
temp_sequences_with_primers="${output_dir}/sequences_with_primers_${identifier}.txt"


# Parse Round 1 primers and build pairs

echo "--- Parsing Round 1 primers from $r1_primers_file ---"

declare -A r1_forward_primers
declare -A r1_reverse_primers
declare -a r1_pair_ids

current_header=""
current_seq=""

while IFS= read -r line || [ -n "$line" ]; do
    # Skip empty lines
    if [ -z "$line" ]; then
        continue
    fi
    
    if [[ $line == ">"* ]]; then
        # Process previous sequence if exists
        if [ -n "$current_header" ] && [ -n "$current_seq" ]; then
            # Clean the sequence of any whitespace
            current_seq=$(echo "$current_seq" | tr -d '[:space:]')
            
            # Extract ALL pair IDs (A, B, C, etc.) from the header
            pair_ids_in_header=($(echo "$current_header" | grep -oP '_[A-Z](?=\s|$|_)' | sed 's/_//g'))
            
            if [[ $current_header == *"Forward"* ]]; then
                # Store forward primer for each pair ID found
                for pair_id in "${pair_ids_in_header[@]}"; do
                    r1_forward_primers[$pair_id]="$current_seq"
                    echo "  Round 1 Forward (pair $pair_id): $current_seq"
                    
                    # Track unique pair IDs
                    if [[ ! " ${r1_pair_ids[@]} " =~ " ${pair_id} " ]]; then
                        r1_pair_ids+=("$pair_id")
                    fi
                done
            elif [[ $current_header == *"Reverse"* ]]; then
                # Store reverse primer for each pair ID found
                for pair_id in "${pair_ids_in_header[@]}"; do
                    r1_reverse_primers[$pair_id]="$current_seq"
                    echo "  Round 1 Reverse (pair $pair_id): $current_seq"
                    
                    # Track unique pair IDs
                    if [[ ! " ${r1_pair_ids[@]} " =~ " ${pair_id} " ]]; then
                        r1_pair_ids+=("$pair_id")
                    fi
                done
            fi
        fi
        
        current_header="$line"
        current_seq=""
    else
        # Accumulate sequence, removing any whitespace
        current_seq="${current_seq}$(echo "$line" | tr -d '[:space:]')"
    fi
done < "$r1_primers_file"

# Process last sequence - this is critical!
if [ -n "$current_header" ] && [ -n "$current_seq" ]; then
    # Clean the sequence
    current_seq=$(echo "$current_seq" | tr -d '[:space:]')
    
    pair_ids_in_header=($(echo "$current_header" | grep -oP '_[A-Z](?=\s|$|_)' | sed 's/_//g'))
    
    if [[ $current_header == *"Forward"* ]]; then
        for pair_id in "${pair_ids_in_header[@]}"; do
            r1_forward_primers[$pair_id]="$current_seq"
            echo "  Round 1 Forward (pair $pair_id): $current_seq"
            
            if [[ ! " ${r1_pair_ids[@]} " =~ " ${pair_id} " ]]; then
                r1_pair_ids+=("$pair_id")
            fi
        done
    elif [[ $current_header == *"Reverse"* ]]; then
        for pair_id in "${pair_ids_in_header[@]}"; do
            r1_reverse_primers[$pair_id]="$current_seq"
            echo "  Round 1 Reverse (pair $pair_id): $current_seq"
            
            if [[ ! " ${r1_pair_ids[@]} " =~ " ${pair_id} " ]]; then
                r1_pair_ids+=("$pair_id")
            fi
        done
    fi
fi


# Parse Round 2 primers if provided

declare -A r2_forward_primers
declare -A r2_reverse_primers
declare -a r2_pair_ids

if [ "$run_round2" = true ] && [ -n "$r2_primers_file" ]; then
    echo "--- Parsing Round 2 primers from $r2_primers_file ---"
    
    current_header=""
    current_seq=""
    
    while IFS= read -r line || [ -n "$line" ]; do
        # Skip empty lines
        if [ -z "$line" ]; then
            continue
        fi
        
        if [[ $line == ">"* ]]; then
            # Process previous sequence if exists
            if [ -n "$current_header" ] && [ -n "$current_seq" ]; then
                # Clean the sequence of any whitespace
                current_seq=$(echo "$current_seq" | tr -d '[:space:]')
                
                # Extract ALL pair IDs (A, B, C, etc.) from the header
                pair_ids_in_header=($(echo "$current_header" | grep -oP '_[A-Z](?=\s|$|_)' | sed 's/_//g'))
                
                if [[ $current_header == *"Forward"* ]]; then
                    # Store forward primer for each pair ID found
                    for pair_id in "${pair_ids_in_header[@]}"; do
                        r2_forward_primers[$pair_id]="$current_seq"
                        echo "  Round 2 Forward (pair $pair_id): $current_seq"
                        
                        # Track unique pair IDs
                        if [[ ! " ${r2_pair_ids[@]} " =~ " ${pair_id} " ]]; then
                            r2_pair_ids+=("$pair_id")
                        fi
                    done
                elif [[ $current_header == *"Reverse"* ]]; then
                    # Store reverse primer for each pair ID found
                    for pair_id in "${pair_ids_in_header[@]}"; do
                        r2_reverse_primers[$pair_id]="$current_seq"
                        echo "  Round 2 Reverse (pair $pair_id): $current_seq"
                        
                        # Track unique pair IDs
                        if [[ ! " ${r2_pair_ids[@]} " =~ " ${pair_id} " ]]; then
                            r2_pair_ids+=("$pair_id")
                        fi
                    done
                fi
            fi
            
            current_header="$line"
            current_seq=""
        else
            # Accumulate sequence, removing any whitespace
            current_seq="${current_seq}$(echo "$line" | tr -d '[:space:]')"
        fi
    done < "$r2_primers_file"

    # Process last sequence - this is critical!
    if [ -n "$current_header" ] && [ -n "$current_seq" ]; then
        # Clean the sequence
        current_seq=$(echo "$current_seq" | tr -d '[:space:]')
        
        pair_ids_in_header=($(echo "$current_header" | grep -oP '_[A-Z](?=\s|$|_)' | sed 's/_//g'))
        
        if [[ $current_header == *"Forward"* ]]; then
            for pair_id in "${pair_ids_in_header[@]}"; do
                r2_forward_primers[$pair_id]="$current_seq"
                echo "  Round 2 Forward (pair $pair_id): $current_seq"
                
                if [[ ! " ${r2_pair_ids[@]} " =~ " ${pair_id} " ]]; then
                    r2_pair_ids+=("$pair_id")
                fi
            done
        elif [[ $current_header == *"Reverse"* ]]; then
            for pair_id in "${pair_ids_in_header[@]}"; do
                r2_reverse_primers[$pair_id]="$current_seq"
                echo "  Round 2 Reverse (pair $pair_id): $current_seq"
                
                if [[ ! " ${r2_pair_ids[@]} " =~ " ${pair_id} " ]]; then
                    r2_pair_ids+=("$pair_id")
                fi
            done
        fi
    fi
fi


# Activate cutadapt environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate orca-seq


# ROUND 1: Linked primer trimming

echo "--- Round 1: Trimming with linked primers ---"

# Build cutadapt command with linked primer pairs
cutadapt_cmd="cutadapt -j $SLURM_CPUS_PER_TASK"

for pair_id in "${r1_pair_ids[@]}"; do
    fwd="${r1_forward_primers[$pair_id]}"
    rev="${r1_reverse_primers[$pair_id]}"
    if [ -n "$fwd" ] && [ -n "$rev" ]; then
        cutadapt_cmd="$cutadapt_cmd -g ${fwd}...${rev}"
        echo "Added linked pair $pair_id: ${fwd}...${rev}"
    else
        echo "Warning: Incomplete pair $pair_id (missing forward or reverse)"
    fi
done

# Add output options
cutadapt_cmd="$cutadapt_cmd --untrimmed-output=$untrimmed_fasta_round1 -o $primerless_fasta_round1 $consensus_file"

# Execute cutadapt
eval $cutadapt_cmd

echo "Round 1 completed."
echo "  Trimmed sequences: $primerless_fasta_round1"
echo "  Untrimmed sequences: $untrimmed_fasta_round1"


# FAILSAFE: Check for residual primers using seqkit locate

if [ -f "$primerless_fasta_round1" ] && [ -s "$primerless_fasta_round1" ]; then
    echo "--- Failsafe: Checking for residual primers with seqkit locate ---"
    
    # Concatenate Round 1 primers
    cat "$r1_primers_file" > "$temp_primers"
    
    # Add Round 2 primers if provided
    if [ -n "$r2_primers_file" ]; then
        cat "$r2_primers_file" >> "$temp_primers"
    fi
    
    primer_count=$(grep -c "^>" "$temp_primers")
    echo "Checking for $primer_count primers..."
    
    # Extract first and last 100 bases from each sequence
    temp_ends="${output_dir}/${identifier}_sequence_ends.fasta"
    seqkit subseq -r 1:100 "$primerless_fasta_round1" > "$temp_ends"
    head "$temp_ends"
    seqkit subseq -r -100:-1 "$primerless_fasta_round1" >> "$temp_ends"
    tail "$temp_ends"
    
    # Use seqkit locate to find sequences with primers in the terminal regions
    seqkit locate -d --pattern-file "$temp_primers" "$temp_ends" > "$temp_primer_locations" 2>/dev/null
    
    # Extract unique sequence IDs that have primer hits (skip header line)
    sequences_with_primers=$(tail -n +2 "$temp_primer_locations" | cut -f1 | uniq)
    
    if [ -n "$sequences_with_primers" ]; then
        seq_count=$(echo "$sequences_with_primers" | wc -l)
        echo "WARNING: Found $seq_count sequence(s) with residual primers!"
        echo "Sample: $identifier"
        echo "Sequences with primers:"
        echo "$sequences_with_primers"
        echo "These sequences still contain primer sequences after Round 1 trimming."
        
        # Extract unique sequence IDs that have primer hits (skip header line)
        cut -f1 "$temp_primer_locations" | sort -u > "$temp_sequences_with_primers"
        
        # Use that file to exclude sequences
        seqkit grep -v -f "$temp_sequences_with_primers" "$primerless_fasta_round1" > "$cleanest_fasta_round1"
        
        removed_count=$(wc -l < "$temp_sequences_with_primers" | tr -d '[:space:]')
        echo "Failsafe FAILED. Residual primers detected, removed ${removed_count} sequences."
        
        # Clean up temp files
        rm -f "$temp_primers" "$temp_primer_locations" "$temp_sequences_with_primers" "$temp_ends"
        
        exit 0
    else
        echo "Failsafe passed. No residual primers detected in Round 1 output."
        # Copy to cleanest output
        cp "$primerless_fasta_round1" "$cleanest_fasta_round1"
    fi
    
    # Clean up temp files after failsafe
    rm -f "$temp_primers" "$temp_primer_locations" "$temp_sequences_with_primers" "$temp_ends"
    
else
    echo "WARNING: No sequences produced in Round 1 for sample: $identifier"
    exit 0
fi


# ROUND 2: Unlinked primer trimming (Only if requested)
if [ "$run_round2" = true ] && [ -s "$untrimmed_fasta_round1" ]; then
    echo "--- Round 2: Trimming untrimmed sequences with unlinked primers ---"
    
    # Build cutadapt command with unlinked primers
    cutadapt_r2_cmd="cutadapt -j $SLURM_CPUS_PER_TASK"
    
    # Add Round 1 primers as unlinked (-g for forward, -a for reverse)
    for pair_id in "${r1_pair_ids[@]}"; do
        fwd="${r1_forward_primers[$pair_id]}"
        rev="${r1_reverse_primers[$pair_id]}"
        
        if [ -n "$fwd" ]; then
            cutadapt_r2_cmd="$cutadapt_r2_cmd -g $fwd"
            echo "Round 2 using forward primer (pair $pair_id): $fwd"
        fi
        
        if [ -n "$rev" ]; then
            cutadapt_r2_cmd="$cutadapt_r2_cmd -a $rev"
            echo "Round 2 using reverse primer (pair $pair_id): $rev"
        fi
    done
    
    # Add Round 2 primers if provided
    if [ ${#r2_forward_primers[@]} -gt 0 ]; then
        echo "Round 2 additional forward primers:"
        for fwd_primer in "${r2_forward_primers[@]}"; do
            cutadapt_r2_cmd="$cutadapt_r2_cmd -g $fwd_primer"
            echo "  $fwd_primer"
        done
    fi
    
    if [ ${#r2_reverse_primers[@]} -gt 0 ]; then
        echo "Round 2 additional reverse primers:"
        for rev_primer in "${r2_reverse_primers[@]}"; do
            cutadapt_r2_cmd="$cutadapt_r2_cmd -a $rev_primer"
            echo "  $rev_primer"
        done
    fi
    
    # Add output and input files
    cutadapt_r2_cmd="$cutadapt_r2_cmd -o $primerless_fasta_round2 $untrimmed_fasta_round1"

    # Execute cutadapt
    eval $cutadapt_r2_cmd
    
    echo "Round 2 completed."
    echo "  Trimmed sequences: $primerless_fasta_round2"
    
elif [ "$run_round2" = true ]; then
    echo "--- Round 2 requested but no untrimmed sequences from Round 1 ---"
    
else
    echo "--- Round 2 skipped (untrimmed sequences discarded) ---"
    # Remove untrimmed file if Round 2 is not requested
    if [ -f "$untrimmed_fasta_round1" ]; then
        rm -f "$untrimmed_fasta_round1"
        echo "  Removed untrimmed sequences file"
    fi
fi


# Summary

echo "Processing completed for $identifier"
echo "Final outputs:"
echo "  Round 1 (clean): $cleanest_fasta_round1"
if [ "$run_round2" = true ] && [ -f "$primerless_fasta_round2" ]; then
    echo "  Round 2 (wrong primers): $primerless_fasta_round2"
fi
