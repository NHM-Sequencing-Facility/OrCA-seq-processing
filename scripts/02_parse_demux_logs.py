#!/usr/bin/env python3
"""Parse OrCA-seq step-02 demultiplexing logs (demux_and_clean_*.log) into summary tables.

Log structure: one log per barcode/dataset, containing
  Round 1  - one cutadapt run splitting the pychopped fastq by SP5 (5') index
  Round 2  - one cutadapt run per SP5 index, splitting by SP27 (3') index

Primary output is final-sample level: one row per SP5 x SP27 combination seen in
the log, with read counts taken from the per-adapter 'Trimmed: N times' line
inside each Round 2 sub-run. A second run-level table captures the Round 1
numbers, which cannot be reconstructed from the sample rows.

Reports only what the log states; no judgement about which index combinations
are expected. Cross-reference the sample sheet downstream if needed.
"""
import argparse
import csv
import glob
import os
import re
import sys

DEFAULT_PATTERN = "demux_and_clean_*.log"
SAMPLES_BASENAME = "demux_samples_summary.tsv"
RUNS_BASENAME = "demux_runs_summary.tsv"

SAMPLE_FIELDS = ["log_file", "dataset", "sp5_index", "sp27_index", "sample_id",
                 "reads", "rc_reads", "pct_of_sp5", "pct_of_run",
                 "sp5_reads", "sp5_assigned", "sp5_pct_assigned",
                 "sp5_unassigned", "note"]

RUN_FIELDS = ["log_file", "dataset", "input_fastq", "status",
              "r1_reads_in", "r1_assigned", "r1_pct_assigned",
              "r1_unassigned", "r1_pct_unassigned",
              "r1_bp_in", "r1_bp_written", "r1_pct_bp_written",
              "n_sp5_found", "n_sp5_processed", "n_sp5_no_reads",
              "n_sp27_seen", "n_combinations", "n_combinations_nonzero",
              "reads_in_combinations",
              "removed_unknown_files", "removed_invalid_files"]

# Sequence: ...; Type: regular 3'; Length: 57; Trimmed: 46162 times; Reverse-complemented: 4565 times
ADAPTER_RE = re.compile(
    r"^=== Adapter (\S+) ===\s*\n+\s*Sequence:.*?Trimmed:\s*(\d+) times"
    r"(?:; Reverse-complemented:\s*(\d+) times)?", re.M | re.S)
SUBRUN_SPLIT_RE = re.compile(r"^  Processing: (\S+)[ \t]*$", re.M)


def _int(s):
    return int(s.replace(",", ""))


def _index_num(name):
    """SP27_009 -> 9. Returns None if the name has no trailing number."""
    m = re.search(r"_(\d+)$", name)
    return int(m.group(1)) if m else None


def parse_summary(txt):
    """Pull the first cutadapt '=== Summary ===' block out of txt."""
    s = {}
    m = re.search(r"^Total reads processed:\s*([\d,]+)", txt, re.M)
    if m:
        s["reads_in"] = _int(m.group(1))
    m = re.search(r"^Reads with adapters:\s*([\d,]+)\s*\(([\d.]+)%\)", txt, re.M)
    if m:
        s["assigned"], s["pct_assigned"] = _int(m.group(1)), float(m.group(2))
    m = re.search(r"^Total basepairs processed:\s*([\d,]+) bp", txt, re.M)
    if m:
        s["bp_in"] = _int(m.group(1))
    m = re.search(r"^Total written \(filtered\):\s*([\d,]+) bp \(([\d.]+)%\)", txt, re.M)
    if m:
        s["bp_written"], s["pct_bp_written"] = _int(m.group(1)), float(m.group(2))
    return s


def parse_adapters(txt):
    """{adapter_name: (trimmed, reverse_complemented)} for one cutadapt run."""
    out = {}
    for name, trimmed, rc in ADAPTER_RE.findall(txt):
        out[name] = (int(trimmed), int(rc) if rc else 0)
    return out


def split_subruns(txt):
    """[(sp5_identifier, chunk_text), ...] for the Round 2 section."""
    parts = SUBRUN_SPLIT_RE.split(txt)
    # parts = [preamble, name1, chunk1, name2, chunk2, ...]
    return list(zip(parts[1::2], parts[2::2]))


def parse(path):
    txt = open(path, errors="replace").read()
    run = {k: "" for k in RUN_FIELDS}
    run["log_file"] = os.path.basename(path)

    m = re.search(r"^Dataset name:\s*(\S+)", txt, re.M)
    dataset = m.group(1) if m else ""
    run["dataset"] = dataset
    m = re.search(r"^Processing:\s*(\S+)", txt, re.M)
    if m:
        run["input_fastq"] = m.group(1)

    # --- status (as reported by the log) ----------------------------------
    err = re.search(r"^Error:\s*(.+)$", txt, re.M)
    if "Pipeline complete!" in txt:
        run["status"] = "OK"
    elif err:
        run["status"] = "FAIL: " + err.group(1).strip()[:60]
    else:
        run["status"] = "FAIL: incomplete/unknown"

    # --- split Round 1 from Round 2 ---------------------------------------
    m = re.search(r"^Round 2: Demultiplexing", txt, re.M)
    r1_txt, r2_txt = (txt[:m.start()], txt[m.start():]) if m else (txt, "")

    # --- Round 1 ----------------------------------------------------------
    s = parse_summary(r1_txt)
    run["r1_reads_in"] = s.get("reads_in", "")
    run["r1_assigned"] = s.get("assigned", "")
    run["r1_pct_assigned"] = s.get("pct_assigned", "")
    if "reads_in" in s and "assigned" in s:
        run["r1_unassigned"] = s["reads_in"] - s["assigned"]
        run["r1_pct_unassigned"] = round(100 - s["pct_assigned"], 2)
    run["r1_bp_in"] = s.get("bp_in", "")
    run["r1_bp_written"] = s.get("bp_written", "")
    run["r1_pct_bp_written"] = s.get("pct_bp_written", "")
    sp5_counts = {k: v[0] for k, v in parse_adapters(r1_txt).items()}

    m = re.search(r"^Found (\d+) identifiers", txt, re.M)
    if m:
        run["n_sp5_found"] = int(m.group(1))

    # --- Round 2 sub-runs -------------------------------------------------
    subruns = split_subruns(r2_txt)
    run["n_sp5_processed"] = len(subruns)

    per_sp5 = {}
    sp27_seen = set()
    for sp5, chunk in subruns:
        sub = parse_summary(chunk)
        adapters = parse_adapters(chunk)
        sp27_seen.update(adapters)
        per_sp5[sp5] = {
            "reads_in": sub.get("reads_in"),
            "assigned": sub.get("assigned"),
            "pct_assigned": sub.get("pct_assigned"),
            "adapters": adapters,
            "no_reads": "No reads processed!" in chunk,
        }
    run["n_sp5_no_reads"] = sum(1 for v in per_sp5.values() if v["no_reads"])
    run["n_sp27_seen"] = len(sp27_seen)

    # Every SP5 and every SP27 index named anywhere in the log, all combinations
    # including zeros.
    all_sp5 = sorted(set(sp5_counts) | set(per_sp5), key=lambda x: (_index_num(x) or 0, x))
    all_sp27 = sorted(sp27_seen, key=lambda x: (_index_num(x) or 0, x))

    r1_total = run["r1_reads_in"] if isinstance(run["r1_reads_in"], int) else None
    samples = []
    for sp5 in all_sp5:
        info = per_sp5.get(sp5)
        sp5_reads = (info["reads_in"] if info and info["reads_in"] is not None
                     else sp5_counts.get(sp5, ""))
        sp5_assigned = info["assigned"] if info and info["assigned"] is not None else ""
        sp5_pct_assigned = info["pct_assigned"] if info and info["pct_assigned"] is not None else ""
        sp5_unassigned = (sp5_reads - sp5_assigned
                          if isinstance(sp5_reads, int) and isinstance(sp5_assigned, int) else "")
        if info is None:
            note = "SP5 not processed in round 2"
        elif info["no_reads"]:
            note = "no reads processed for this SP5"
        else:
            note = ""
        for sp27 in all_sp27:
            reads, rc = (info["adapters"].get(sp27, (0, 0)) if info else (0, 0))
            samples.append({
                "log_file": run["log_file"],
                "dataset": dataset,
                "sp5_index": sp5,
                "sp27_index": sp27,
                # mirrors the pipeline's output filename: {SP27}_{SP5}_{dataset}
                "sample_id": f"{sp27}_{sp5}_{dataset}" if dataset else f"{sp27}_{sp5}",
                "reads": reads,
                "rc_reads": rc,
                "pct_of_sp5": (round(100 * reads / sp5_reads, 2)
                               if isinstance(sp5_reads, int) and sp5_reads else ""),
                "pct_of_run": (round(100 * reads / r1_total, 3) if r1_total else ""),
                "sp5_reads": sp5_reads,
                "sp5_assigned": sp5_assigned,
                "sp5_pct_assigned": sp5_pct_assigned,
                "sp5_unassigned": sp5_unassigned,
                "note": note,
            })

    run["n_combinations"] = len(samples)
    run["n_combinations_nonzero"] = sum(1 for r in samples if r["reads"])
    run["reads_in_combinations"] = sum(r["reads"] for r in samples)

    m = re.search(r"^Removed (\d+) files with 'unknown'", txt, re.M)
    if m:
        run["removed_unknown_files"] = int(m.group(1))
    m = re.search(r"^Removed (\d+) files with invalid SP27 indices", txt, re.M)
    if m:
        run["removed_invalid_files"] = int(m.group(1))

    return run, samples


def resolve_out(out, indir, basename):
    if out is None:
        return os.path.join(indir, basename)
    if out.endswith(os.sep) or os.path.isdir(out) or not out.lower().endswith((".tsv", ".txt", ".csv")):
        return os.path.join(out, basename)
    return out


def write_tsv(path, fields, rows):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-d", "--dir", default=".",
                    help=f"directory to scan for {DEFAULT_PATTERN} (default: cwd)")
    ap.add_argument("-o", "--out", default=None,
                    help=f"output directory or .tsv path (default: <dir>/{SAMPLES_BASENAME})")
    ap.add_argument("-p", "--pattern", default=DEFAULT_PATTERN,
                    help=f"log glob pattern (default: {DEFAULT_PATTERN})")
    ap.add_argument("-r", "--recursive", action="store_true", help="search subdirectories too")
    ap.add_argument("--no-runs", action="store_true", help="skip the run-level TSV")
    ap.add_argument("-q", "--quiet", action="store_true", help="do not echo tables to stdout")
    args = ap.parse_args()

    indir = os.path.abspath(args.dir)
    if not os.path.isdir(indir):
        sys.exit(f"Not a directory: {indir}")

    pat = os.path.join(indir, "**", args.pattern) if args.recursive else os.path.join(indir, args.pattern)
    paths = sorted(glob.glob(pat, recursive=args.recursive))
    if not paths:
        sys.exit(f"No logs matching {args.pattern} in {indir}")

    runs, samples = [], []
    for p in paths:
        run, srows = parse(p)
        runs.append(run)
        samples.extend(srows)

    runs.sort(key=lambda x: str(x["dataset"] or "zzz"))
    samples.sort(key=lambda x: (str(x["dataset"]), _index_num(x["sp5_index"]) or 0,
                                _index_num(x["sp27_index"]) or 0))

    samples_out = resolve_out(args.out, indir, SAMPLES_BASENAME)
    write_tsv(samples_out, SAMPLE_FIELDS, samples)
    runs_out = None
    if not args.no_runs:
        runs_out = os.path.join(os.path.dirname(os.path.abspath(samples_out)), RUNS_BASENAME)
        write_tsv(runs_out, RUN_FIELDS, runs)

    if not args.quiet:
        print("\t".join(SAMPLE_FIELDS))
        for r in samples:
            print("\t".join(str(r[k]) for k in SAMPLE_FIELDS))

    n_ok = sum(1 for r in runs if r["status"] == "OK")
    nz = sum(1 for r in samples if r["reads"])
    print(f"\n{len(runs)} logs parsed ({n_ok} OK), {len(samples)} combinations "
          f"({nz} with reads) -> {samples_out}", file=sys.stderr)
    if runs_out:
        print(f"run-level summary -> {runs_out}", file=sys.stderr)


if __name__ == "__main__":
    main()