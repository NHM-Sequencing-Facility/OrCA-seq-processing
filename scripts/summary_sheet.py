#!/usr/bin/env python3

import re
import pandas as pd
import glob
import os
import openpyxl 
from openpyxl.styles import Font, Alignment, PatternFill

# This will need to change depending on the master script, can implement argparser
sample_sheet= "/home/sophr2/OrCA-seq-processing/OC SSI Indexing_28S.csv"
cutadapt = "/hpc/scratch/sophr2/DA_2026-04-01_OC-SSI-Run2/no_sample_id/20260401_1152_X2_FBF82196_4b65a8cc/fastq_pass/barcode02/28S_log/demux_samples_summary.tsv"
amplicon_sorter = "/hpc/scratch/sophr2/DA_2026-04-01_OC-SSI-Run2/no_sample_id/20260401_1152_X2_FBF82196_4b65a8cc/fastq_pass/barcode02/28S_log/amplicon_sorter_summary.tsv"
cleaning_primers = "/hpc/scratch/sophr2/DA_2026-04-01_OC-SSI-Run2/no_sample_id/20260401_1152_X2_FBF82196_4b65a8cc/fastq_pass/barcode02/28S_log/trim_summary.tsv"
seqkit_stats = "/hpc/scratch/sophr2/DA_2026-04-01_OC-SSI-Run2/no_sample_id/20260401_1152_X2_FBF82196_4b65a8cc/fastq_pass/barcode02/rRNA_genes/seqkit_stats.tsv"
barrnap = "/hpc/scratch/sophr2/DA_2026-04-01_OC-SSI-Run2/no_sample_id/20260401_1152_X2_FBF82196_4b65a8cc/fastq_pass/barcode02/rRNA_genes/"
blast = "/hpc/scratch/sophr2/DA_2026-04-01_OC-SSI-Run2/no_sample_id/20260401_1152_X2_FBF82196_4b65a8cc/fastq_pass/barcode02/rRNA_genes/phylum_fastas/"

# Define the amplicon
amplicon= "28S"

# parse the index combinations from the identifier column in amplicon sorter and cleaning primers
def parse_id(identifier):
    """
    Extract SP27 and SP5 index combinations from an identifier string.
    e.g. 'SP27_001_SP5_001_barcode01_reads-coi' -> ('SP27_001', 'SP5_001')
    """
    sp27 = re.search(r"SP27_\d+", identifier)
    sp5  = re.search(r"SP5_\d+", identifier)
    return sp27.group() if sp27 else None, sp5.group() if sp5 else None


# Load the sample sheet as the master reference in which everything else will join onto
df_samples = pd.read_csv(sample_sheet)

# Load in cutadapt results
df_cutadapt = pd.read_csv(cutadapt, sep="\t")

df_cutadapt["sample F"] = df_cutadapt["sp5_index"]
df_cutadapt["sample R"] = df_cutadapt["sp27_index"]


# Load in amplicon sorter results
df_amplicon_sorter = pd.read_csv(amplicon_sorter, sep="\t")

# parse_id through amplicon sorter 
df_amplicon_sorter[["sample R", "sample F"]] = df_amplicon_sorter["identifier"].apply(
    lambda x: pd.Series(parse_id(x))
)

# Load in cleaning_primers results
df_cleaning_primers = pd.read_csv(cleaning_primers, sep="\t")

# parse_id through cleaning primers
df_cleaning_primers[["sample R", "sample F"]] = df_cleaning_primers["identifier"].apply(
    lambda x: pd.Series(parse_id(x))
)

# Load in seqkit_stats results
df_seqkit_stats = pd.read_csv(seqkit_stats, sep="\t")

# Filter to relevant amplicon file, and parse sample index
df_seqkit_stats = df_seqkit_stats[df_seqkit_stats["file"].str.contains(f"_{amplicon}.fa")]
df_seqkit_stats[["sample R", "sample F"]] = df_seqkit_stats["file"].apply(
    lambda x: pd.Series(parse_id(x))
)
df_seqkit_stats = df_seqkit_stats.drop_duplicates(subset=["sample F", "sample R"], keep="first")

# Parse the barrnap fields if amplicon is 28s or 18S
if amplicon in ["28S", "18S"]:
    # Find all .gff3 files recursively under the barrnap directory
    gff3_files = glob.glob(os.path.join(barrnap, "**", "*.gff3"), recursive=True)   
    rows = []
    for f in gff3_files:
        with open(f) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 9:
                    continue
                rows.append({
                    "filename": os.path.basename(f),
                    "seqname": parts[0],
                    "barrnap_evalue": parts[5],
                    "attributes": parts[8]
                })
    
    df_barrnap = pd.DataFrame(rows)
    df_barrnap[["sample R", "sample F"]] = df_barrnap["filename"].apply(
        lambda x: pd.Series(parse_id(x))
    )

    # Extract readcount from sequence name and keep highest readcount per sample
    df_barrnap["readcount"] = (
        df_barrnap["seqname"]
        .str.extract(r"readcount_(\d+)")[0]
        .astype(float)
    )
    df_barrnap = (
        df_barrnap
        .sort_values("readcount", ascending=False)
        .groupby(["sample F", "sample R"], as_index=False)
        .first()
        .drop(columns=["readcount", "seqname"])
    )

# Load BLAST results -accepts a single file or a directory of TSV files
if os.path.isfile(blast):
    df_blast = pd.read_csv(blast, sep="\t")
elif os.path.isdir(blast):
    blast_files = glob.glob(os.path.join(blast, "*.tsv"))
    df_blast = pd.concat(
        [pd.read_csv(f, sep="\t") for f in blast_files],
        ignore_index=True
    )
else:
    raise ValueError(f"blast results not found: {blast}")

# parse index combinations from qseqid
df_blast[["sample R", "sample F"]] = df_blast["qseqid"].apply(
    lambda x: pd.Series(parse_id(x))
)

# Extract the readcount from the sequence header
df_blast["readcount"] = df_blast["qseqid"].apply(
    lambda x: int(re.search(r"readcount_(\d+)", x).group(1))
)

# Sort by bitscore to extract the top BLAST hit for each HRCB
df_blast = df_blast.sort_values("bitscore", ascending=False)
df_blast = df_blast.drop_duplicates(subset="qseqid", keep="first")
idx = df_blast.groupby(["sample F", "sample R"])["readcount"].idxmax()
df_hrcb = df_blast.loc[idx]

df_amplicon_sorter = df_amplicon_sorter.rename(columns={"status": "amplicon_status"})
df_cleaning_primers = df_cleaning_primers.rename(columns={"status": "primers_status"})

# Merge all the results into a table
df_output = (
    df_samples
    .merge(df_cutadapt[["sample F", "sample R", "reads", "rc_reads", "sp5_reads"]], on=["sample F", "sample R"], how="left")
    .merge(df_amplicon_sorter[["sample F", "sample R", "sorted_status", "reads_in_window", "used_reads", "n_gene_groups", "n_consensus", "top_consensus_reads", "top_consensus_pct", "top2_pct"]], on=["sample F", "sample R"], how="left")
    .merge(df_cleaning_primers[["sample F", "sample R", "trim_status", "primer_fwd", "primer_rev"]], on=["sample F", "sample R"], how="left")
    .merge(df_seqkit_stats[["sample F", "sample R", "sum_len", "min_len", "avg_len", "max_len"]], on=["sample F", "sample R"], how="left")
    .merge(df_hrcb[["sample F", "sample R", "qseqid", "sseqid", "pident", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "stitle"]], on=["sample F", "sample R"], how="left")
)

# write out barrnap
if amplicon in ["28S", "18S"]:
    df_output = df_output.merge(
        df_barrnap[["sample F", "sample R", "barrnap_evalue", "attributes"]],
        on=["sample F", "sample R"],
        how="left"
    )


# Save the spreadsheet
output_path = f"pipeline_summary_{amplicon}.xlsx"
df_output.to_excel(output_path, index=False)

# Define column groups for the merged header row for user visibility
groups = [
    ("SAMPLE INFO",      ["Plate Location", "Sample Number", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Forward Indexing Primer", "Reverse Indexing Primer", "sample F", "sample R"]),
    ("CUTADAPT",         ["reads", "rc_reads", "sp5_reads"]),
    ("AMPLICON SORTER",  ["sorter_status", "reads_in_window", "used_reads", "n_gene_groups", "n_consensus", "top_consensus_reads", "top_consensus_pct", "top2_pct"]),
    ("CLEANING PRIMERS", ["trim_status", "primer_fwd", "primer_rev"]),
    ("SEQKIT STATS",     ["sum_len", "min_len", "avg_len", "max_len"]),
    ("BLAST",            ["qseqid", "sseqid", "pident", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "stitle"]),
]

if amplicon in ["28S", "18S"]:
    groups.append(("BARRNAP", ["barrnap_evalue", "attributes"]))

# Add merged group headers above the column names
def add_group_headers(filepath, df, groups):
    wb = openpyxl.load_workbook(filepath)
    ws = wb.active
    ws.insert_rows(1)
    col_names = list(df.columns)
    for group_name, cols in groups:
        positions = [col_names.index(c) + 1 for c in cols if c in col_names]
        if not positions:
            continue
        first_col = min(positions)
        last_col  = max(positions)
        cell = ws.cell(row=1, column=first_col, value=group_name)
        cell.font      = Font(bold=True)
        cell.alignment = Alignment(horizontal="left")
        if first_col != last_col:
            ws.merge_cells(start_row=1, start_column=first_col,
                           end_row=1,   end_column=last_col)
    wb.save(filepath)

add_group_headers(output_path, df_output, groups)
print(f"Summary spreadsheet written to {output_path}")

