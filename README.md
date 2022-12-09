# Analysis of Transposable Elements (TEs) expression from RNA-Seq data for project Boulard - Stork

* Until featureCounts output, analysis is performed using snakemake pipeline in snake-make/TE_RNASeq.Snakefile.
* Snakefile is run on the cluster using pipeline_wrapper.sh in src/sh
* Snakefile uses conda environments in env/conda and singularity containers built with recipes in singularity/recipes
* config file for Snakefile and for SLURM are in config/
* important for reproducibility: 
  - conda version 4.8.3
  - singularity version used to build containers is 3.5.3
  - snakemake version to run Snakefile is 5.9.1
* Using the Rdata files produced in Snakefile, I then continue the analysis in the Rmd files, which can run both on the cluster and locally on a personal computer (having at least 16Gb of RAM).

**Note on mapping rate**: the % of uniquely mapped reads using parameters as in STAR_align rule in TE_RNASeq.Snakefile (which are the ones suggested by Deborah Bourc'His in the review cited in the Snakefile), was unusually low (~50%), with an unusually high (~30%) % of short reads and too-many-mismatches reads. Therefore I tried to understand where this was coming from by making a STAR_align_debug rule where I use same parameters except for the fact that I allow splicing: this is done by both allowing spliced reads and by allowing mates to be far away as per default (instead of max 350 bp apart). With these modifications, the % of uniquely mapped reads increased to ~75-80%, therefore I do not think there is a problem in the STAR_align mapping, but simply that intron-containing fragments end up being trimmed by STAR (thus becoming too short) and not mapped anyway.
