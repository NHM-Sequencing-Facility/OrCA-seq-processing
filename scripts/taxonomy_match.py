#!/usr/bin/env python3
"""
taxonomy_match.py

Compares the expected taxonomy of each sample (the Phylum, Class, Order,
Family, Genus and Species columns from the sample sheet) with the taxonomy of
its top BLAST hit (the staxid column), using a local NCBI taxdump.

Two columns are added to the summary spreadsheet made by summary_sheet.py:
  tax_match   TRUE  - one of the expected taxa is in the hit's NCBI lineage
              FALSE - none of the expected taxa are in the hit's NCBI lineage
              NA    - no usable staxid, or no expected taxonomy for the sample
  rank_match  The lowest expected rank that matched (Species ... Phylum).
              Blank when tax_match is FALSE, NA when tax_match is NA.

Matching:
  The hit's lineage is the scientific name of its taxid and of every ancestor
  up to the root. Starting at the lowest filled expected rank, each expected 
  name is looked for in that lineage, moving up a rank at a time until one is
  found. Blank expected ranks are skipped.

The input file is overwritten. If it already has tax_match / rank_match
columns (from an earlier run), they are replaced.

Usage:
  python3 taxonomy_match.py --input pipeline_summary_COI.xlsx --taxdump /path/to/taxdump

  --taxdump is a directory containing names.dmp, nodes.dmp, merged.dmp and
  delnodes.dmp, from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
"""

import argparse
import os
import sys
import tempfile
from collections import Counter
from copy import copy

import openpyxl
from openpyxl.styles import Font, Alignment

# Expected rank columns, lowest first (the order they are checked in)
EXPECTED_RANKS = ["Species", "Genus", "Family", "Order", "Class", "Phylum"]

TAXID_COL   = "staxid"
TITLE_COL   = "stitle"
SAMPLE_COL  = "Sample Number"
BLAST_GROUP = "BLAST"
NEW_GROUP   = "TAXONOMIC RECONCILIATION"
NEW_COLS    = ["tax_match", "rank_match"]

GROUP_ROW      = 1   # merged group headers (SAMPLE INFO, BLAST, ...)
HEADER_ROW     = 2   # column names
FIRST_DATA_ROW = 3

TAXDUMP_FILES = ["nodes.dmp", "names.dmp", "merged.dmp", "delnodes.dmp"]

# Expected-taxonomy cell values treated as blank
BLANK_VALUES = {"na", "n/a", "nan"}

# NA reasons, in the order they are reported
NA_NO_TAXID    = "no staxid (no BLAST hit)"
NA_BAD_TAXID   = "staxid is 0 / N/A / not a taxid (database without taxonomy?)"
NA_DELETED     = "taxid deleted from NCBI taxonomy (delnodes.dmp)"
NA_NOT_IN_DUMP = "taxid not in taxdump (taxdump older than the BLAST database?)"
NA_NO_EXPECTED = "no expected taxonomy for the sample"
NA_REASONS     = [NA_NO_TAXID, NA_BAD_TAXID, NA_DELETED, NA_NOT_IN_DUMP, NA_NO_EXPECTED]


# -----------------------------------------------------------------------------
# Small helpers
# -----------------------------------------------------------------------------
def normalise(name):
    """Lower-cased name with whitespace collapsed, or None if blank."""
    if name is None:
        return None
    s = " ".join(str(name).split())
    if not s or s.casefold() in BLANK_VALUES:
        return None
    return s.casefold()


def parse_taxid(value):
    """
    Returns (taxid, None) for a valid taxid, or (None, reason) otherwise.
    """
    if value is None or (isinstance(value, str) and not value.strip()):
        return None, NA_NO_TAXID
    if isinstance(value, bool):
        return None, NA_BAD_TAXID
    if isinstance(value, int):
        taxid = value
    else:
        try:
            f = float(value)
        except (TypeError, ValueError):
            return None, NA_BAD_TAXID
        if f != f:  # NaN
            return None, NA_NO_TAXID
        if not f.is_integer():
            return None, NA_BAD_TAXID
        taxid = int(f)
    if taxid <= 0:
        return None, NA_BAD_TAXID
    return taxid, None


def dmp_fields(line, maxsplit=-1):
    """Splits a line of an NCBI .dmp file ("a\\t|\\tb\\t|\\n") into fields."""
    line = line.rstrip("\n")
    if line.endswith("\t|"):
        line = line[:-2]
    return line.split("\t|\t", maxsplit)


# -----------------------------------------------------------------------------
# NCBI taxonomy
# -----------------------------------------------------------------------------
def load_parents(path):
    """nodes.dmp -> {taxid: parent_taxid}"""
    parents = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            f = dmp_fields(line, 2)
            parents[int(f[0])] = int(f[1])
    return parents


def load_merged(path):
    """merged.dmp -> {old_taxid: new_taxid}"""
    merged = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            f = dmp_fields(line)
            merged[int(f[0])] = int(f[1])
    return merged


def load_deleted(path):
    """delnodes.dmp -> {taxid, ...}"""
    with open(path, encoding="utf-8", errors="replace") as fh:
        return {int(dmp_fields(line)[0]) for line in fh if line.strip()}


def load_scientific_names(path, wanted_taxids, wanted_names):
    """
    Reads names.dmp once, keeping only scientific names.
    Returns:
      names - {taxid: scientific name} for taxids in wanted_taxids
      found - the normalised names from wanted_names that exist as a
              scientific name anywhere in the NCBI taxonomy
    """
    names, found = {}, set()
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            f = dmp_fields(line)
            if len(f) < 4 or f[3] != "scientific name":
                continue
            taxid = int(f[0])
            if taxid in wanted_taxids:
                names[taxid] = f[1]
            n = normalise(f[1])
            if n in wanted_names:
                found.add(n)
    return names, found


def lineage_taxids(taxid, parents):
    """taxid followed by all of its ancestors, up to and including the root."""
    out = [taxid]
    for _ in range(1000):  # guards against a malformed nodes.dmp
        parent = parents.get(taxid)
        if parent is None or parent == taxid:
            break
        out.append(parent)
        taxid = parent
    return out


# -----------------------------------------------------------------------------
# Spreadsheet layout
# -----------------------------------------------------------------------------
def header_positions(ws):
    """{column name: column index} from the column-name row."""
    positions = {}
    for col in range(1, ws.max_column + 1):
        value = ws.cell(HEADER_ROW, col).value
        if value is not None:
            positions.setdefault(str(value).strip(), col)
    return positions


def take_group_headers(ws):
    """
    Records the group headers in the first row as [(group name, [column names])]
    and removes them (unmerges and clears the row), so columns can be deleted
    and inserted safely. They are drawn again by draw_group_headers().
    """
    names = {col: name for name, col in header_positions(ws).items()}
    spans, covered = [], set()

    row_ranges = [r for r in ws.merged_cells.ranges
                  if r.min_row == GROUP_ROW and r.max_row == GROUP_ROW]
    for r in row_ranges:
        spans.append((r.min_col, r.max_col, ws.cell(GROUP_ROW, r.min_col).value))
        covered.update(range(r.min_col, r.max_col + 1))
    for col in range(1, ws.max_column + 1):
        value = ws.cell(GROUP_ROW, col).value
        if col not in covered and value is not None:
            spans.append((col, col, value))

    groups = []
    for first, last, name in sorted(spans, key=lambda s: s[0]):
        if name is None:
            continue
        cols = [names[c] for c in range(first, last + 1) if c in names]
        groups.append((str(name).strip(), cols))

    for r in row_ranges:
        ws.unmerge_cells(str(r))
    for col in range(1, ws.max_column + 1):
        cell = ws.cell(GROUP_ROW, col)
        cell.value = None
        cell.font = Font()
        cell.alignment = Alignment()
    return groups


def draw_group_headers(ws, groups):
    """Writes and merges the group headers, as summary_sheet.py does."""
    positions = header_positions(ws)
    for name, cols in groups:
        cols_at = sorted(positions[c] for c in cols if c in positions)
        if not cols_at:
            continue
        first, last = cols_at[0], cols_at[-1]
        if last - first + 1 != len(cols_at):
            print(f"  Warning: the columns of group '{name}' are no longer next to "
                  f"each other; the header is not merged", file=sys.stderr)
            last = first
        cell = ws.cell(GROUP_ROW, first, value=name)
        cell.font = Font(bold=True)
        cell.alignment = Alignment(horizontal="left")
        if first != last:
            ws.merge_cells(start_row=GROUP_ROW, start_column=first,
                           end_row=GROUP_ROW, end_column=last)


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Match the expected taxonomy of each sample against its top "
                    "BLAST hit (staxid) and add tax_match / rank_match columns "
                    "to the summary spreadsheet (overwritten in place).")
    parser.add_argument("--input", required=True,
                        help="Summary spreadsheet (.xlsx) from summary_sheet.py")
    parser.add_argument("--taxdump", required=True,
                        help="Directory with names.dmp, nodes.dmp, merged.dmp and delnodes.dmp")
    args = parser.parse_args()

    # --- Check inputs -------------------------------------------------------
    if not os.path.isfile(args.input):
        sys.exit(f"Error: input spreadsheet not found: {args.input}")
    missing = [f for f in TAXDUMP_FILES
               if not os.path.isfile(os.path.join(args.taxdump, f))]
    if missing:
        sys.exit(f"Error: missing from --taxdump {args.taxdump}: {', '.join(missing)}\n"
                 f"Download and unpack https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz")

    print(f"Taxonomy matching: {args.input}")
    wb = openpyxl.load_workbook(args.input)
    ws = wb.active

    positions = header_positions(ws)
    required = EXPECTED_RANKS + [TAXID_COL]
    absent = [c for c in required if c not in positions]
    if absent:
        msg = f"Error: column(s) not found in row {HEADER_ROW} of {args.input}: {', '.join(absent)}"
        if TAXID_COL in absent:
            msg += ("\nThe summary sheet needs a staxid column: re-run 06_BLAST.sh (with "
                    "staxid in BLAST_COLS) and summary_sheet.py.")
        sys.exit(msg)

    # --- Remove group headers and any columns from an earlier run -----------
    groups = take_group_headers(ws)
    old_cols = sorted((positions[c] for c in NEW_COLS if c in positions), reverse=True)
    for col in old_cols:
        ws.delete_cols(col)
    if old_cols:
        print("  Replacing tax_match / rank_match from an earlier run")
    groups = [(name, [c for c in cols if c not in NEW_COLS])
              for name, cols in groups if name.upper() != NEW_GROUP]
    positions = header_positions(ws)

    # --- Read the rows ------------------------------------------------------
    rows = []  # (row number, sample label, expected [(rank, name)], raw staxid)
    spelling = {}  # normalised name -> name as written in the sheet (for the log)
    for r in range(FIRST_DATA_ROW, ws.max_row + 1):
        if all(ws.cell(r, c).value is None for c in range(1, ws.max_column + 1)):
            continue
        sample = ws.cell(r, positions[SAMPLE_COL]).value if SAMPLE_COL in positions else None
        sample = str(sample) if sample is not None else f"row {r}"
        expected = []
        for rank in EXPECTED_RANKS:  # lowest first
            value = ws.cell(r, positions[rank]).value
            name = normalise(value)
            if name:
                expected.append((rank, name))
                spelling.setdefault(name, " ".join(str(value).split()))
        rows.append((r, sample, expected, ws.cell(r, positions[TAXID_COL]).value))

    # --- Load the taxonomy and build each hit's lineage ---------------------
    print(f"  Loading NCBI taxonomy from {args.taxdump}")
    merged  = load_merged(os.path.join(args.taxdump, "merged.dmp"))
    deleted = load_deleted(os.path.join(args.taxdump, "delnodes.dmp"))
    parents = load_parents(os.path.join(args.taxdump, "nodes.dmp"))

    resolved = {}             # raw taxid -> (current taxid or None, NA reason)
    lineages = {}             # current taxid -> [taxids]
    merged_used = {}          # old taxid -> new taxid
    na_taxids = {NA_DELETED: set(), NA_NOT_IN_DUMP: set()}
    for _, _, _, raw in rows:
        taxid, reason = parse_taxid(raw)
        if taxid is None or taxid in resolved:
            continue
        current = merged.get(taxid, taxid)
        if current != taxid:
            merged_used[taxid] = current
        if taxid in deleted or current in deleted:
            resolved[taxid] = (None, NA_DELETED)
            na_taxids[NA_DELETED].add(taxid)
        elif current not in parents:
            resolved[taxid] = (None, NA_NOT_IN_DUMP)
            na_taxids[NA_NOT_IN_DUMP].add(taxid)
        else:
            resolved[taxid] = (current, None)
            if current not in lineages:
                lineages[current] = lineage_taxids(current, parents)
    del parents

    wanted_taxids = {t for lineage in lineages.values() for t in lineage}
    wanted_names  = {name for _, _, expected, _ in rows for _, name in expected}
    sci_names, found_names = load_scientific_names(
        os.path.join(args.taxdump, "names.dmp"), wanted_taxids, wanted_names)
    lineage_names = {
        taxid: {normalise(sci_names[t]) for t in lineage if t in sci_names}
        for taxid, lineage in lineages.items()
    }

    # --- Match --------------------------------------------------------------
    results = {}  # row number -> (tax_match, rank_match)
    counts, rank_counts, na_counts = Counter(), Counter(), Counter()
    for r, sample, expected, raw in rows:
        taxid, reason = parse_taxid(raw)
        if taxid is not None:
            taxid, reason = resolved[taxid]
        if reason is None and not expected:
            reason = NA_NO_EXPECTED
        if reason is not None:
            results[r] = ("NA", "NA")
            counts["NA"] += 1
            na_counts[reason] += 1
            continue

        names = lineage_names[taxid]
        rank = next((rank for rank, name in expected if name in names), None)
        if rank is not None:
            results[r] = (True, rank)
            counts["TRUE"] += 1
            rank_counts[rank] += 1
        else:
            results[r] = (False, None)
            counts["FALSE"] += 1

    # --- Insert the new columns after the last BLAST column -----------------
    blast_cols = next((cols for name, cols in groups if name.upper() == BLAST_GROUP), [])
    anchor_cols = [positions[c] for c in blast_cols if c in positions]
    if not anchor_cols:
        anchor_cols = [positions[c] for c in (TITLE_COL, TAXID_COL) if c in positions]
    anchor = max(anchor_cols)

    ws.insert_cols(anchor + 1, amount=len(NEW_COLS))
    header_style = ws.cell(HEADER_ROW, anchor)._style
    for i, name in enumerate(NEW_COLS, start=1):
        cell = ws.cell(HEADER_ROW, anchor + i, value=name)
        cell._style = copy(header_style)
    for r, (tax_match, rank_match) in results.items():
        ws.cell(r, anchor + 1, value=tax_match)
        ws.cell(r, anchor + 2, value=rank_match)

    # Put the new group straight after BLAST (order only matters for reading)
    insert_at = next((i + 1 for i, (name, _) in enumerate(groups)
                      if name.upper() == BLAST_GROUP), len(groups))
    groups.insert(insert_at, (NEW_GROUP, list(NEW_COLS)))
    draw_group_headers(ws, groups)

    # --- Save (to a temporary file first, so a failure can't corrupt the input)
    out_dir = os.path.dirname(os.path.abspath(args.input))
    fd, tmp_path = tempfile.mkstemp(suffix=".xlsx", dir=out_dir)
    os.close(fd)
    try:
        wb.save(tmp_path)
        os.replace(tmp_path, args.input)
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)

    # --- Summary ------------------------------------------------------------
    print(f"  Rows: {len(rows)}")
    by_rank = ", ".join(f"{rank} {rank_counts[rank]}"
                        for rank in EXPECTED_RANKS if rank_counts[rank])
    print(f"  tax_match TRUE:  {counts['TRUE']}" + (f"  ({by_rank})" if by_rank else ""))
    print(f"  tax_match FALSE: {counts['FALSE']}")
    print(f"  tax_match NA:    {counts['NA']}")
    for reason in NA_REASONS:
        if na_counts[reason]:
            line = f"    {na_counts[reason]:>4}  {reason}"
            if na_taxids.get(reason):
                line += ": " + ", ".join(str(t) for t in sorted(na_taxids[reason]))
            print(line)
    if merged_used:
        print("  Merged taxids updated: " +
              ", ".join(f"{old} -> {new}" for old, new in sorted(merged_used.items())))

    # Expected names that are not a scientific name anywhere in NCBI can never
    # match (spelling, or a synonym / name NCBI doesn't use)
    not_found = {}
    for _, sample, expected, _ in rows:
        for rank, name in expected:
            if name not in found_names:
                not_found.setdefault((rank, name), []).append(sample)
    if not_found:
        print("  Expected names not found in the NCBI taxonomy (check spelling/synonyms):")
        for (rank, name), samples in sorted(not_found.items(),
                                            key=lambda kv: (EXPECTED_RANKS.index(kv[0][0]), kv[0][1])):
            print(f"    {rank}: {spelling[name]}  ({', '.join(samples)})")

    print(f"  Written to {args.input}")


if __name__ == "__main__":
    main()
