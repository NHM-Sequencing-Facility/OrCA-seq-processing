# Processing OrCA-seq barcoding data
## *This is the code used to process COI/28S/18S amplicon data from the 2025 NERC/NHM Field workshop on the taxonomy and natural history of freshwater and limno-terrestrial meiofauna*

## Installation and setup
* You'll need this software:
  * pychopper v2.7.0 : `conda create -n pychopper bioconda::pychopper`
  * cutadapt v4.9.0: `conda create -n cutadapt bioconda::cutadapt`
  * Amplicon Sorter (no releases, and you will need to [download the python code to the `scripts/auxiliary_code` directory](https://github.com/avierstr/amplicon_sorter/blob/master/amplicon_sorter.py))

    `conda create -n amplicon_sorter && conda install bioconda::python-edlib biopython matplotlib` 
  * seqkit v2.9.0: `conda create -n seqkit bioconda::seqkit`
  * barrnap v0.9: `conda create -n barrnap bioconda::barrnap==0.9`
  * pybarrnap v0.5.1: `conda create -n pybarrnap -c conda-forge -c bioconda pybarrnap -c bioconda infernal`
  * seqtk v1.4-r122: `conda create -n seqtk bioconda::seqtk`
  * hmmer v3.1b2: `conda create -n hmmer bioconda::hmmer`
  * bedtools v2.31.1: `conda create -n bedtools bioconda::bedtools`

## How to structure
- Run the scripts in the parent directory to the cloned github repo (so the filepaths remain relative)
- This workflow is written for a SLURM HPC system (I used the CropDiversity cluster)
- To run this code, git clone this whole repo, keeping the adapter directory and the scripts directory in the relative places they're in so that the codes run smoothly (or, edit the code yourself only in your cloned remote repo!

## The workflow
![bioinformatic workflow for processing DNA barcoding, from raw data to reorienting, demuxing, consensus building, cleaning, getting gene sequences, and selecting the best hits](schematics/DNA_Barcoding_workflow.png)

  1. `01_pychopper.sh`: `sbatch $0 /path/to/dataset/raw_fastq_file`
     * To reorient and quality-score trim cDNA reads into the same orientation
     * This script takes in 2 extra files, which are hardcoded into the script:
       * one which has the orientation of adapters sequences listed (`M13_config_for_pychopper.txt`)
       * one which has the non-variable sequences for the adapters (`M13_seqs_for_pychopper.fa`)
  
  2. `02_cutadapt_loop.sh`: `sbatch $0 /path/to/dataset/pychopped/pychopped_sample1.fastq.gz`
       * A loop which first demultiplexes the reads based on the 5' SP5 primers, and then, one-by-one, demultiplexes each output from the SP5 demultiplexing by the SP27 3' primers. It also trims adapters on the fly.
       * This script takes in 2 extra files, which are hardcoded into the script: the 5' and 3' adapter sequences, `M13_amplicon_indices_forward.fa` and the reverse complement of the 3' sequences (since we reoriented before with pychopper, `M13_amplicon_indices_reverse_rc.fa`)
       * You may notice that we do bin reads into adapter combinations which do not exist! Though this is inefficient, it's not the end of the world, since we have sequenced deeply enough per individual that enough reads are binned into the correct adapter combinations
       * These non-existing adapter combos are removed after demuxing, along with the 'unknown' bins, since we don't want to analyse them later on. 
  
  3. `03_amplicon_sorter.sh`: `sbatch $0 -input /path/to/dataset/demuxed/2nd_round [-min <int> -max <int> -prefix <amplicon>]`
     * The workflow is more dynamic here, since the user may be analysing different amplicons for us. In _our_ wet-lab protocol:
       * For the rRNAs, we use 2 sets of primers to amplify overlapping segments of the nuclear rRNA cistron, so a sequence can be anywhere over 3Kb in length
       * For the COIs, we use 2 alternate options for a forward primer, and one option for a reverse primer, to amplify 'redundant' sequences of a COI barcode segment (the 2 alternate options for the forward primer allow for matching more taxa than just one), so a sequence can be between 300bp-900bp in length
     * The amplicon sorter step splits samples' reads by size and clusters them by sequence
     * The variables (min length, max length, prefix name for amplicon type, and input folder) are user-defined :)
  
  4. `04a_cleaning_primers.sh`: `sbatch $0 <amplicon_sorted_dir> <amplicon_type> --r1-primers <file> --r2-primers <file> [--run-round2]`
       * This script takes each clustered (amplicon-sorted) file and removes the primer sequences from the amplicon sequences, since these are synthetic
       * The user can define the primer sequences to remove based on the amplicon they are sequencing, but if more than one amplicon was sequenced in a run, the other primer sets can also be submitted as a back check
       * Sequences with lone primers, mismatched primers, or >1 primer of a kind are removed in a failsafe
       * Optionally, the user can choose to trim 'untrimmed' sequences with alternate primers
  
  5. `05a_pybarrnap_rDNA_extract.sh`: `sbatch $0 /path/to/dataset/primerless` and `05b_reorganise_COIs.sh`: `sbatch $0 /path/to/dataset/primerless`
       * 5a script uses pybarrnap version 0.5.1. It takes the assembled contigs and uses an covariance model based on Rfam(14.10) to extract sequences matching 28S and 18S rDNA profiles from our amplicon contigs.
       * 5b is a straightforward script copying over the cleaned, clustered/non-redundant primerless COIs from the primerless directory to a COI directory for clarity

## Further notes
The analysis can stop here, if you like. The scripts 06, 07 (both), 08, and 09, as well as those in the R analysis directory, can be used for command-line automated BLASTn searching and reorganising contigs for tree building or species delimitation. The scripts listed were all used in our paper. 
  
## Primer schematics
### rRNA amplification
![schematic for amplification of 18S/5.8S/partial 28S rRNA cistron, showing wiggly purple line representing DNA, yellow boxes where genes are (overlaid on the purple DNA), and red arrows labelled with primer sequence name at the positions where the primers sit to amplify genes, not to scale](schematics/18Setc.png)
![schematic for amplification of full 28S rRNA gene, showing wiggly purple line representing DNA, yellow boxes where genes are (overlaid on the purple DNA), and red arrows labelled with primer sequence name at the positions where the primers sit to amplify genes, not to scale](schematics/28S_full.png)

### COI amplification
![schematic for amplification of partial COI gene, showing wiggly purple line representing DNA, yellow boxes where COI gene is (overlaid on the purple DNA), and red arrows labelled with primer sequence name at the positions where the primers sit to amplify COI, not to scale](schematics/COI.png)
