#!/usr/bin/env python3
"""Parse primer-trimming wrapper logs (trim_<jobid>_<idx>.log) into a summary table.

Assumes Round 1 only (linked primer pair A); no Round 2 handling.
Counts are consensus sequences, not raw reads, hence seqs_* column naming
(cutadapt calls them "reads").
"""
import argparse
import csv
import glob
import os
import re
import sys

FIELDS = [
    "log_file", "job_id", "array_id", "array_total", "identifier", "amplicon",
    "status", "primer_fwd", "primer_rev",
    "seqs_in", "seqs_with_adapters", "pct_with_adapters",
    "seqs_untrimmed", "pct_untrimmed", "seqs_written", "pct_written",
    "bp_in", "bp_written", "pct_bp_written",
    "trimmed_5p", "trimmed_3p",
    "failsafe", "n_primers_checked", "final_seqs", "flag",
]

OK_FAILSAFE = "Failsafe passed. No residual primers detected in Round 1 output."


def _int(s):
    return int(s.replace(",", ""))


def _search(pat, txt, flags=re.M):
    return re.search(pat, txt, flags)


def parse(path):
    txt = open(path, errors="replace").read()
    r = {k: "" for k in FIELDS}
    r["log_file"] = os.path.basename(path)

    # --- provenance ---------------------------------------------------------
    m = re.match(r"trim_(\d+)_(\d+)\.log$", r["log_file"])
    if m:
        r["job_id"] = m.group(1)
    m = _search(r"^Array Task ID:\s*(\d+)\s*/\s*(\d+)", txt)
    if m:
        r["array_id"], r["array_total"] = int(m.group(1)), int(m.group(2))
    elif re.match(r"trim_\d+_(\d+)\.log$", r["log_file"]):
        # fall back to filename index if the log died before printing the header
        r["array_id"] = int(re.match(r"trim_\d+_(\d+)\.log$", r["log_file"]).group(1))

    # --- identity: prefer explicit lines, fall back to the Processing path ---
    m = _search(r"^Sample identifier:\s*(\S+)", txt)
    if m:
        r["identifier"] = m.group(1)
    m = _search(r"^Amplicon type:\s*(\S+)", txt)
    if m:
        r["amplicon"] = m.group(1)
    if not (r["identifier"] and r["amplicon"]):
        m = _search(r"^Processing:\s*(\S+)", txt)
        if m:
            parts = m.group(1).rstrip("/").split("/")
            if len(parts) >= 2:
                r["identifier"] = r["identifier"] or parts[-2]
                r["amplicon"] = r["amplicon"] or parts[-1]

    # --- status -------------------------------------------------------------
    err = _search(r"^Error:\s*(.+)$", txt)
    if "Processing completed for" in txt:
        r["status"] = "OK"
    elif err:
        msg = err.group(1).strip()
        if "No consensus file found" in msg:
            r["status"] = "FAIL: no consensus file"
        else:
            r["status"] = "FAIL: " + msg[:60]
    else:
        r["status"] = "FAIL: incomplete/unknown"

    # --- primers ------------------------------------------------------------
    m = _search(r"^\s*Round 1 Forward \(pair A\):\s*(\S+)", txt)
    if m:
        r["primer_fwd"] = m.group(1)
    m = _search(r"^\s*Round 1 Reverse \(pair A\):\s*(\S+)", txt)
    if m:
        r["primer_rev"] = m.group(1)

    # --- cutadapt summary ---------------------------------------------------
    m = _search(r"^Total reads processed:\s*([\d,]+)", txt)
    if m:
        r["seqs_in"] = _int(m.group(1))
    m = _search(r"^Reads with adapters:\s*([\d,]+)\s*\(([\d.]+)%\)", txt)
    if m:
        r["seqs_with_adapters"], r["pct_with_adapters"] = _int(m.group(1)), float(m.group(2))
    m = _search(r"^Reads discarded as untrimmed:\s*([\d,]+)\s*\(([\d.]+)%\)", txt)
    if m:
        r["seqs_untrimmed"], r["pct_untrimmed"] = _int(m.group(1)), float(m.group(2))
    m = _search(r"^Reads written \(passing filters\):\s*([\d,]+)\s*\(([\d.]+)%\)", txt)
    if m:
        r["seqs_written"], r["pct_written"] = _int(m.group(1)), float(m.group(2))
    m = _search(r"^Total basepairs processed:\s*([\d,]+) bp", txt)
    if m:
        r["bp_in"] = _int(m.group(1))
    m = _search(r"^Total written \(filtered\):\s*([\d,]+) bp \(([\d.]+)%\)", txt)
    if m:
        r["bp_written"], r["pct_bp_written"] = _int(m.group(1)), float(m.group(2))
    m = _search(r"5' trimmed:\s*(\d+) times; 3' trimmed:\s*(\d+) times", txt)
    if m:
        r["trimmed_5p"], r["trimmed_3p"] = int(m.group(1)), int(m.group(2))

    # --- failsafe -----------------------------------------------------------
    m = _search(r"^Checking for (\d+) primers", txt)
    if m:
        r["n_primers_checked"] = int(m.group(1))
    if OK_FAILSAFE in txt:
        r["failsafe"] = "passed"
    elif "--- Failsafe" in txt:
        r["failsafe"] = "NOT_PASSED"
    else:
        r["failsafe"] = "not_run"

    # Round 1 only, so the final cleaned file holds whatever cutadapt wrote.
    r["final_seqs"] = r["seqs_written"]

    # --- interpretive flag --------------------------------------------------
    flags = []
    if r["status"] != "OK":
        flags.append("failed")
    else:
        if r["failsafe"] == "NOT_PASSED":
            flags.append("RESIDUAL PRIMERS / failsafe not passed")
        elif r["failsafe"] == "not_run":
            flags.append("failsafe not run")
        if isinstance(r["seqs_written"], int) and r["seqs_written"] == 0:
            flags.append("no sequences survived")
        if isinstance(r["pct_untrimmed"], float) and r["pct_untrimmed"] > 0:
            flags.append(f"partial trimming ({r['pct_untrimmed']}% untrimmed)")
        if isinstance(r["pct_bp_written"], float) and r["pct_bp_written"] < 70:
            flags.append("low bp retention")
        if (isinstance(r["trimmed_5p"], int) and isinstance(r["trimmed_3p"], int)
                and r["trimmed_5p"] != r["trimmed_3p"]):
            flags.append("asymmetric 5'/3' trimming")
        if isinstance(r["seqs_in"], int) and 0 < r["seqs_in"] < 3:
            flags.append("very few input consensus seqs")
        if not flags:
            flags.append("clean")
    r["flag"] = "; ".join(flags)
    return r


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-d", "--dir", default=".", help="directory containing trim_*.log (default: cwd)")
    ap.add_argument("-p", "--pattern", default="trim_*.log", help="glob pattern (default: trim_*.log)")
    ap.add_argument("-o", "--out", default="trim_summary.tsv", help="output TSV path")
    ap.add_argument("-q", "--quiet", action="store_true", help="do not echo table to stdout")
    args = ap.parse_args()

    paths = sorted(glob.glob(os.path.join(args.dir, args.pattern)))
    if not paths:
        sys.exit(f"No logs matching {args.pattern} in {args.dir}")

    rows = [parse(p) for p in paths]
    rows.sort(key=lambda x: (str(x["identifier"] or "zzz"), str(x["amplicon"])))

    outdir = os.path.dirname(os.path.abspath(args.out))
    os.makedirs(outdir, exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    if not args.quiet:
        print("\t".join(FIELDS))
        for r in rows:
            print("\t".join(str(r[k]) for k in FIELDS))

    n_ok = sum(1 for r in rows if r["status"] == "OK")
    print(f"\n{len(rows)} logs parsed: {n_ok} OK, {len(rows) - n_ok} failed -> {args.out}",
          file=sys.stderr)


if __name__ == "__main__":
    main()