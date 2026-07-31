### Output TSV column definitions:

#log_file
Basename of the parsed log, e.g. trim_154688_82.log. The key back to the raw file.

#job_id
SLURM job ID, taken from the filename (trim_<job_id>_<idx>.log). Empty if the filename doesn't match that form.

#array_id
Array task index, from Array Task ID: N / 96 inside the log. Falls back to the filename index if the log died before printing that header.

#array_total
Total array size (the 96). Same for every row in a run; a sanity check that we have all the logs.

#identifier
Sample ID from Sample identifier:. For early-exit logs, recovered as the second-to-last component of the Processing: path.

#amplicon
Amplicon type from Amplicon type: (e.g. COI, coi). Recovered from the last path component on early exit. Note it preserves the log's case, so COI and coi won't group together.

#status
OK if Processing completed for is present; FAIL: no consensus file for the missing-input case; FAIL: <message> (first 60 chars) for any other Error: line; FAIL: incomplete/unknown if the log just stops (typically a timeout or OOM kill).

#primer_fwd
Round 1 forward primer, pair A, as parsed from the primer FASTA. Useful for confirming all samples used the same primer set.

#primer_rev
Round 1 reverse primer, pair A.

#n_primers_checked
Number of primers the seqkit failsafe searched for (Checking for N primers...) — 2 for a single linked pair.

#seqs_in
Consensus sequences in the input FASTA (Total reads processed).

#seqs_with_adapters
Sequences where at least one primer was matched.

#pct_with_adapters
The above as a percentage of seqs_in.

#seqs_untrimmed
Sequences with no primer match, diverted to the untrimmed file and discarded.

#pct_untrimmed
The above as a percentage. Anything above 0 means primer detection failed on some sequences.

#seqs_written
Sequences written to the trimmed output.

#pct_written
The above as a percentage of seqs_in.

#bp_in
Total basepairs in the input.

#bp_written
Total basepairs retained after trimming.

#pct_bp_written
Retention as a percentage. Expect roughly 90–95% for a ~700 bp COI amplicon with 51 bp of primer removed; much lower suggests over-trimming or unexpectedly short input.

#trimmed_5p
Times the 5' (forward) primer was trimmed.

#trimmed_3p
Times the 3' (reverse) primer was trimmed. Should equal trimmed_5p for full-length linked matches; a difference means some sequences were trimmed at one end only.

#final_seqs
Sequences in the final cleaned FASTA. With Round 1 only this is a copy of seqs_written — it exists so the column survives if Round 2 is added.

#failsafe
Passed only on the exact success string; NOT_PASSED if a Failsafe section ran but didn't report success; not_run if no Failsafe section exists (i.e. the log died first).

#flag
Semicolon-separated interpretive notes: clean (nothing of concern), failed, RESIDUAL PRIMERS / failsafe not passed, failsafe not run, no sequences survived, partial trimming (N% untrimmed), low bp retention (<70%), asymmetric 5'/3' trimming, very few input consensus seqs (<3).