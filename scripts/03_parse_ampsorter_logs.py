#!/usr/bin/env python3
"""Parse amplicon_sorter wrapper logs into a summary table."""
import argparse
import csv
import glob
import os
import re
import sys

DEFAULT_PATTERN = "AS_array_*.log"
DEFAULT_BASENAME = "amplicon_sorter_summary.tsv"

FIELDS = ["array_id", "identifier", "dataset", "status", "reads_in_fastq", "reads_in_window",
          "used_reads", "subsampled", "ssg", "n_gene_groups", "n_consensus",
          "top_consensus_reads", "top_consensus_pct", "top2_pct", "total_assigned_pct",
          "nogroup_reads", "nogroup_pct", "flag"]


def parse(path):
    txt = open(path, errors="replace").read()
    r = {k: "" for k in FIELDS}
    m = re.search(r"AS_array_(\d+)", os.path.basename(path))
    r["array_id"] = m.group(1) if m else ""
    m = re.search(r"^Identifier:\s*(\S+)", txt, re.M);  r["identifier"] = m.group(1) if m else ""
    m = re.search(r"^Dataset:\s*(\S+)", txt, re.M);     r["dataset"] = m.group(1) if m else ""
    if "Processing complete!" in txt:
        r["status"] = "OK"
    elif "Number of usable sequences is lower than 5" in txt:
        r["status"] = "FAIL: <5 usable reads"
    elif "not created" in txt:
        r["status"] = "FAIL: no consensus"
    else:
        r["status"] = "FAIL: unknown"
    m = re.search(r"contains (\d+) reads\.", txt)
    if m: r["reads_in_fastq"] = int(m.group(1))
    # "Reading 10000 out of 45733 sequences between..."  or "reading all 583 out of 583"
    m = re.search(r"[Rr]eading (?:all )?(\d+) out of (\d+) sequences between", txt)
    if m:
        r["used_reads"], r["reads_in_window"] = int(m.group(1)), int(m.group(2))
        r["subsampled"] = "Y" if r["used_reads"] < r["reads_in_window"] else "N"
    m = re.search(r"Estimated ssg = (\d+)", txt)
    if m: r["ssg"] = int(m.group(1))
    # gene-level groups surviving the <5-read filter
    genes = re.findall(r"^  \S+_(\d+)\.group contains (\d+) sequences \(([\d.]+)%\)", txt, re.M)
    if genes: r["n_gene_groups"] = len(genes)
    # final consensus sequences: "--> ID_g_s.fasta contains N sequences (P% of total)"
    cons = [(int(a), float(b)) for a, b in
            re.findall(r"^--> \S+\.fasta contains (\d+) sequences \(([\d.]+)% of total\)", txt, re.M)]
    if cons:
        cons.sort(key=lambda x: -x[0])
        r["n_consensus"] = len(cons)
        r["top_consensus_reads"] = cons[0][0]
        r["top_consensus_pct"] = cons[0][1]
        r["top2_pct"] = round(sum(c[1] for c in cons[:2]), 2)
        r["total_assigned_pct"] = round(sum(c[1] for c in cons), 2)
    m = re.search(r"(\d+) sequences were not assigned in groups.*?\(([\d.]+)% of total\)", txt)
    if m:
        r["nogroup_reads"], r["nogroup_pct"] = int(m.group(1)), float(m.group(2))
    # interpretive flag
    flags = []
    if r["status"] != "OK":
        flags.append("failed")
    else:
        if r["n_consensus"] == 1:
            flags.append("single clean consensus")
        elif isinstance(r["top_consensus_pct"], float):
            if r["top_consensus_pct"] >= 80:
                flags.append("dominant + minor variants")
            elif r["top2_pct"] >= 75:
                flags.append("two co-dominant clusters")
            else:
                flags.append("MIXED / no dominant cluster")
        if isinstance(r["n_consensus"], int) and r["n_consensus"] >= 8:
            flags.append("high cluster count")
        if isinstance(r["total_assigned_pct"], float) and r["total_assigned_pct"] < 95:
            flags.append("low read recovery")
    r["flag"] = "; ".join(flags)
    return r


def resolve_out(out, indir):
    """-o may be a directory or an explicit .tsv path."""
    if out is None:
        return os.path.join(indir, DEFAULT_BASENAME)
    if out.endswith(os.sep) or os.path.isdir(out) or not out.lower().endswith((".tsv", ".txt", ".csv")):
        return os.path.join(out, DEFAULT_BASENAME)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-d", "--dir", default=".",
                    help=f"directory to scan for {DEFAULT_PATTERN} (default: cwd)")
    ap.add_argument("-o", "--out", default=None,
                    help=f"output directory or .tsv path (default: <dir>/{DEFAULT_BASENAME})")
    ap.add_argument("-p", "--pattern", default=DEFAULT_PATTERN,
                    help=f"log glob pattern (default: {DEFAULT_PATTERN})")
    ap.add_argument("-r", "--recursive", action="store_true",
                    help="search subdirectories too")
    ap.add_argument("-q", "--quiet", action="store_true", help="do not echo table to stdout")
    args = ap.parse_args()

    indir = os.path.abspath(args.dir)
    if not os.path.isdir(indir):
        sys.exit(f"Not a directory: {indir}")

    if args.recursive:
        paths = sorted(glob.glob(os.path.join(indir, "**", args.pattern), recursive=True))
    else:
        paths = sorted(glob.glob(os.path.join(indir, args.pattern)))
    if not paths:
        sys.exit(f"No logs matching {args.pattern} in {indir}")

    rows = [parse(p) for p in paths]
    rows.sort(key=lambda x: x["identifier"] or "zzz")

    out = resolve_out(args.out, indir)
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    if not args.quiet:
        print("\t".join(FIELDS))
        for r in rows:
            print("\t".join(str(r[k]) for k in FIELDS))

    n_ok = sum(1 for r in rows if r["status"] == "OK")
    print(f"\n{len(rows)} logs parsed: {n_ok} OK, {len(rows) - n_ok} failed -> {out}",
          file=sys.stderr)


if __name__ == "__main__":
    main()