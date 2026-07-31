### Output TSV column definitions:

#array_id
The SLURM array task number, taken from the log filename (AS_array_150333.log). Not biologically meaningful; it's how you get back to the specific job if you need to re-inspect stderr or resubmit.

#identifier
The sample identifier from the wrapper's Identifier: line. Encodes the dual-index combination and ONT barcode (SP27_001_SP5_004_barcode17 = SP27 index 001, SP5 index 004, native barcode 17). This is the join key to your sample sheet.

#dataset
The ONT barcode directory the sample came from, from the Dataset: line. Redundant with the tail of identifier here, but it survives if you rename samples.

#status
Did the run produce a consensus? OK if the log reached Processing complete! Otherwise a parsed failure reason: <5 usable reads if amplicon_sorter quit at the length-filter stage, or no consensus if it exited without writing consensusfile.fasta and the reason isn't in the wrapper log (it'll be in results.txt).

#reads_in_fastq
Total reads in the demultiplexed input, from contains N reads. Raw yield for that index combination, before any filtering.

#reads_in_window
How many of those fell between -min and -max, from Reading X out of Y sequences between Abp and Bbp. The gap between this and reads_in_fastq is your length-filter loss. This should be relatively small in every sample, which is expected for a clean amplicon.

#used_reads
How many reads amplicon_sorter actually clustered. Caps at maxreads (10,000). This is the denominator for every percentage in the table, which is the single most important thing: a cluster reported at 60% is 60% of the 10,000 subsample, not of the(e.g) 46,162 reads in the FASTQ.

#subsampled
Y when used_reads < reads_in_window, i.e. the cap was hit and the run saw only part of the library. A flag that cluster proportions are estimates.

#ssg (similar_species_groups)
The per-dataset identity threshold amplicon_sorter estimated for the initial species-level sub-grouping, from Estimated ssg = N. Higher means it started with finer splits. It's data-derived rather than user-set, so it varies (e.g. 90–99).

#n_gene_groups
How many gene-level groups (≥80% identity tier) survived the ≥5-read filter, counted from the X_N.group contains ... lines. In a single-locus COI run you expect 1. More than one usually means off-target amplification or something divergent enough to fail the 80% threshold, not necessarily a different gene.

#n_consensus
The number of consensus sequences actually written, counted from the --> ID_g_s.fasta contains ... lines. This is the headline result: how many distinct sequence variants amplicon_sorter thinks are in the sample after all merging at the 96% consensus threshold.

#top_consensus_reads / top_consensus_pct
Read support for the most abundant consensus, absolute and as a percentage of used_reads. For a single-specimen barcode you want this near 100%.

#top2_pct
The two most abundant consensuses combined. Useful for distinguishing "one taxon plus noise" from "two real things": if the top cluster is 60% but top-two is 95%, that's a two-component sample, whereas 60% / 79% means a long tail.

#total_assigned_pct
All consensus percentages summed, i.e. what fraction of the clustered reads ended up in any consensus. Below ~95% means a lot of reads were discarded as too divergent or in sub-threshold groups, which is a quality warning independent of how many clusters there are.

#nogroup_reads / nogroup_pct
Reads that failed to join any group at all, written to *_nogroup_unique.fasta. Usually a handful of chimeras or very low-quality reads. Note this is not the complement of total_assigned_pct: reads can also be lost in groups that fell below the minimum-size filter, so nogroup + assigned won't sum to 100%.

#flag
03_parse_ampsorter_logs.py summary judgement, not parsed from the log. Applied from n_consensus, top_consensus_pct and top2_pct using these thresholds: ≥80% top = dominant, ≥75% top-two = co-dominant, else mixed; ≥8 clusters = high cluster count; <95% assigned = low recovery.