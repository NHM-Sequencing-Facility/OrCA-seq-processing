#!/usr/bin/env python3
"""
index_mapping.py

Reads both forwards and reverse index sequences from sample sheet (in the form of an excel file),
searches them against the FASTA file, fetches the SP5 and SP27 combination from the FASTA file and pastes
them into the columns 'Sample F' and 'Sample R'.

Sequences need to be reformatted in order to be searched against the FASTA file where
sample sheet sequences are:
<constant region> <index sequence>*<nucleotide>

Usage:
python index_mapping.py --fasta <indexes>.fasta --sample <sample_sheet>.xlsx --sheet <SheetName>

Output: <sample sheet>_mapped.xlsx (written in the current directory)
"""

import argparse
import os
import sys
import openpyxl

# Helper functions

def parse_fasta(path):
    """Return dict {header_id: sequence} from a FASTA file."""
    seqs = {}
    current_id = None
    current_seq = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    seqs[current_id] = "".join(current_seq)
                current_id = line[1:].split()[0]   # first word after '>'
                current_seq = []
            else:
                current_seq.append(line)
    if current_id is not None:
        seqs[current_id] = "".join(current_seq)
    return seqs

def strip_formatting(seq):
    """Remove spaces and asterisks from a formatted sequence."""
    if seq is None:
        return ""
    return str(seq).replace(" ", "").replace("*", "").upper().strip()

def build_reverse_lookup(fasta_seqs):
    """
    Return dict {raw_sequence: header_id}.
    Warns if two headers share the same sequence.
    """
    lookup = {}
    for id_, seq in fasta_seqs.items():
        seq_clean = seq.upper().strip()
        if seq_clean in lookup:
            print("WARNING: duplicate sequence for '{}' and '{}' - "
                  "keeping first.".format(lookup[seq_clean], id_))
        else:
            lookup[seq_clean] = id_
    return lookup

def get_or_create_col(ws, headers, col_name):
    """
    Return 1-based column index for col_name.
    If it doesn't exist, append it to row 1 and return the new index.
    """
    if col_name in headers:
        return headers.index(col_name) + 1
    # append
    new_col = len(headers) + 1
    ws.cell(row=1, column=new_col).value = col_name
    headers.append(col_name)
    return new_col


# Main

def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--fasta",   required=True,
                    help="FASTA file containing SP5/SP27 index sequences")
    ap.add_argument("--sample",   required=True,
                    help="Excel sample sheet (.xlsx)")
    ap.add_argument("--sheet",   default=None,
                    help="Sheet name to edit (default: first sheet)")
    ap.add_argument("--fwd-col", default="Forward Indexing Primer",
                    help="Column header for forward primer sequences")
    ap.add_argument("--rev-col", default="Reverse Indexing Primer",
                    help="Column header for reverse primer sequences")
    args = ap.parse_args()


## Load the FASTA file and build the reverse lookups

    raw_seqs = parse_fasta(args.fasta)
    print("Loaded {} sequences from FASTA.".format(len(raw_seqs)))
    reverse_lookup = build_reverse_lookup(raw_seqs)

## Load workbook

    wb = openpyxl.load_workbook(args.sample)
    ws = wb[args.sheet] if args.sheet else wb.active
    print("Editing sheet: '{}'".format(ws.title))

## Create column for index IDs

    headers = [cell.value for cell in ws[1]]
    while headers and headers[-1] is None:
        headers.pop()

    def require_col(name):
        if name not in headers:
            sys.exit("ERROR: column '{}' not found.\nColumns: {}".format(name, headers))
        return headers.index(name) + 1

    col_fwd_seq  = require_col(args.fwd_col)
    col_rev_seq  = require_col(args.rev_col)

    # Create output columns if they don't exist yet
    col_sample_f = get_or_create_col(ws, headers, "sample F")
    col_sample_r = get_or_create_col(ws, headers, "sample R")


## Fill IDs into Columns

    filled_f  = 0
    filled_r  = 0
    missing   = []

    for row_idx in range(2, ws.max_row + 1):
        fwd_raw = ws.cell(row=row_idx, column=col_fwd_seq).value
        rev_raw = ws.cell(row=row_idx, column=col_rev_seq).value

        if not fwd_raw and not rev_raw:
            continue   # blank row

        fwd_clean = strip_formatting(fwd_raw)
        rev_clean = strip_formatting(rev_raw)

        fwd_id = reverse_lookup.get(fwd_clean)
        rev_id = reverse_lookup.get(rev_clean) or reverse_lookup.get(rev_clean[:-1])

        if fwd_id:
            ws.cell(row=row_idx, column=col_sample_f).value = fwd_id
            filled_f += 1
        elif fwd_clean:
            missing.append("row {}: forward sequence not found in FASTA "
                           "(first 20 nt: {})".format(row_idx, fwd_clean[:20]))

        if rev_id:
            ws.cell(row=row_idx, column=col_sample_r).value = rev_id
            filled_r += 1
        elif rev_clean:
            missing.append("row {}: reverse sequence not found in FASTA "
                           "(first 20 nt: {})".format(row_idx, rev_clean[:20]))

## Save sample sheet

    base, ext = os.path.splitext(args.sample)
    sheet_label = args.sheet.replace(" ", "_") if args.sheet else ws.title.replace(" ", "_")
    out_path = base + "_" + sheet_label + "_mapped" + ext
    wb.save(out_path)

    print("Done. sample F filled: {}  |  sample R filled: {}".format(filled_f, filled_r))
    print("Output saved to: {}".format(out_path))

    if missing:
        print("\nWarnings ({} sequences not matched):".format(len(missing)))
        for m in missing:
            print("  " + m)

if __name__ == "__main__":
    main()