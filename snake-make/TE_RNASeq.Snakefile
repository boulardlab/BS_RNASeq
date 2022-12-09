########################################################################################
# Pipeline for analysis of RNA-Seq data for analysis of Tranposable Elements Expression
# Date: 01/11/20
# Author: Sara Formichetti
# email: sara.formichetti@embl.it
# run: all specified in shell wrapper for the specific experiment's analysis
#######################################################################################

#####################
# Imports
#####################

import pandas
import os
import yaml

#####################
# Defining shell to use
#####################

shell.executable("/bin/bash")

#####################
# Create tmp dir if it does not exist
#####################
if not os.path.isdir('./tmp'):
    shell("mkdir tmp/")

#####################
# Functions
#####################

def read_samplesTable(samplesTable):
    samplesData = pandas.read_csv(samplesTable)
    # Verify column names
    if not {'SRR','sample','condition', 'group_or_time_point', 'biol_rep','tech_rep','run','library_layout','read_length'}.issubset(samplesData.columns.values):
	    raise KeyError("The samples file must contain the following named columns (name of the column[content format]): SRR[SRR ID i.e. string],sample[string with no spaces],condition[string with no spaces],group_or_time_point[string with no spaces],biol_rep[number],tech_rep[number],run[number],library_layout[SINGLE - |PAIRED - ],read_length[number]; please add this column and just write always 'none'(without quotes) under 'SRR' if they do not come from SRA and always '1'(without quotes) under 'condition/group_or_time_point/biol_rep/tech_rep' if there is no more than one group.")
    return samplesData

######################
# Config Variables
######################

# Checking that provided configfile contains all necessary parameters 

with open('config/TE_RNASeq_config.yaml') as f:
    my_config_dict = yaml.safe_load(f)

target_config_keys = ['experiment_folder', 'input_samples_table', 'fastq_input_dir', 'seq_output_dir', 'genome_dir', 'maj_chr_bed', 'annotation_dir', 'gene_annotation_bed', 'TE_annotation_bed', 'analysis_output_dir', 'tmp_dir', 'envs_folder']
missing_params = [ele for ele in target_config_keys if ele not in list(my_config_dict)] 
if len(missing_params) == 0:
    print("Configfile contains all necessary variables")
else :
    print("Configfile misses following parameters: ", str(missing_params))
    sys.exit()

# Getting value of variable containing samples' table file name

try:
  input_samples_table = config["input_samples_table"]
except KeyError:
    print("The parameter \"input_samples_table\" has not been defined in the configfile. The pipeline cannot start.")
    
#################################
# Reading input samples' table
#################################

inputSamplesData = read_samplesTable(input_samples_table)

########################
# Variables definition
########################

# Splitting the table into single or paired end experiments

index_single = inputSamplesData['library_layout'] == 'SINGLE - '
index_paired = inputSamplesData['library_layout'] == 'PAIRED - '
samplesData_single = inputSamplesData[index_single]
samplesData_paired = inputSamplesData[index_paired]

# Output files names

SINGLESAMPLES = samplesData_single['sample'].tolist()
PAIREDSAMPLES = samplesData_paired['sample'].tolist()

###############################################################################
# Rules
###############################################################################

# Defining singularity container to use 
#singularity: config["containers_folder"] + ADD SIF IMAGE 

######################
# Rule all
######################

rule all:
  input:
    expand("{path}qc/{sample}_R{mate}.fastqc.{ext}", path=config["seq_output_dir"], sample=PAIREDSAMPLES, ext=["html", "zip"], mate=["1","2"]),
    expand("{path}alignment/STAR/{sample}_amgm350.TE_only.bam", path=config["seq_output_dir"], sample=PAIREDSAMPLES),
    expand("{path}featureCounts/{sample}_TE_counts.txt", path=config["analysis_output_dir"], sample=PAIREDSAMPLES),
    expand("{path}featureCounts/{sample}_TE_counts_fl_only.txt", path=config["analysis_output_dir"], sample=PAIREDSAMPLES)

rule fastqc_my_fastq:
  input:
    expand("{path}{{sample}}_R{{mate}}.fastq.gz", path=config["fastq_input_dir"])
  params:
    output_dir=expand("{path}qc/", path=config["seq_output_dir"])
  output:
    expand("{path}qc/{{sample}}_R{{mate}}.fastqc.{ext}", path=config["seq_output_dir"], ext=["html", "zip"])
  threads: 6
  conda: config["envs_folder"] + "TE_RNASeq.yml" 
  shell:
    "fastqc -t 6 -o {params.output_dir} {input}"

rule sort_my_fastq:
  input:
    expand("{path}{{sample}}_R{{mate}}.fastq.gz", path=config["fastq_input_dir"])
  output:
    expand("{path}{{sample}}_R{{mate}}.sorted.fastq.gz", path=config["fastq_input_dir"])
  threads: 6
  shell:
    """
    zcat {input} | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" | gzip - > {output}
    """

rule trimgalore_my_fastq:
  input:
    R1=expand("{path}{{sample}}_R1.sorted.fastq.gz", path=config["fastq_input_dir"]),
    R2=expand("{path}{{sample}}_R2.sorted.fastq.gz", path=config["fastq_input_dir"])
  output:
    R1_trimmed_fq=expand("{path}trimming/{{sample}}_R1.sorted_val_1.fq.gz", path=config["seq_output_dir"]),
    R2_trimmed_fq=expand("{path}trimming/{{sample}}_R2.sorted_val_2.fq.gz", path=config["seq_output_dir"])
  params:
    output_dir=expand("{path}trimming/", path=config["seq_output_dir"]),
  threads: 6
  conda: config["envs_folder"] + "TE_RNASeq.yml"
  shell: 
    "trim_galore -j 6 -q 20 --stringency 1 -e 0.1 --length 20 --fastqc -o {params.output_dir} --paired {input.R1} {input.R2}" # trimgalore default parameters
 
rule STAR_align:
  input:
    R1=rules.trimgalore_my_fastq.output.R1_trimmed_fq,
    R2=rules.trimgalore_my_fastq.output.R2_trimmed_fq
  output:
    bam=expand("{path}alignment/STAR/{{sample}}_amgm350_Aligned.out.bam", path=config["seq_output_dir"])
  params:
    genome_folder=expand("{path}STAR/index_2.7.5c", path=config["genome_dir"]), # index built with STAR 2.7.c5 using GRCm38.primary_assembly.genome.fa.gz 
    output_prefix=expand("{path}alignment/STAR/{{sample}}_amgm350_", path=config["seq_output_dir"])
  threads: 6
  conda: config["envs_folder"] + "TE_RNASeq.yml"
  shell:
    "STAR --runThreadN {threads} --genomeDir {params.genome_folder} --readFilesCommand zcat --readFilesIn {input.R1} {input.R2} --outFileNamePrefix {params.output_prefix} --runMode alignReads --outSAMtype BAM Unsorted --outFilterMultimapNmax 5000 --outSAMmultNmax 1 --outFilterMismatchNmax 3 --outMultimapperOrder Random --winAnchorMultimapNmax 5000 --alignEndsType EndToEnd --alignIntronMax 1 --alignMatesGapMax 350 --seedSearchStartLmax 30 --alignTranscriptsPerReadNmax 30000 --alignWindowsPerReadNmax 30000 --alignTranscriptsPerWindowNmax 300 --seedPerReadNmax 3000 --seedPerWindowNmax 300 --seedNoneLociPerWindow 1000 --outReadsUnmapped Fastx --outSAMunmapped Within" ## parameters are from Teissandier et al. Mobile DNA (2019), alignment to mouse genome in random mode with STAR

#rule STAR_align_debug:
#  input:
#    R1=rules.trimgalore_my_fastq.output.R1_trimmed_fq,
#    R2=rules.trimgalore_my_fastq.output.R2_trimmed_fq
#  output:
#    bam=expand("{path}alignment/STAR_debug_3/{{sample}}_{{mm_rate}}_Aligned.out.bam", path=config["seq_output_dir"])
#  params:
#    genome_folder=expand("{path}STAR/index_2.7.5c", path=config["genome_dir"]), # index built with STAR 2.7.c5 using GRCm38.primary_assembly.genome.fa.gz 
#    output_prefix=expand("{path}alignment/STAR_debug_3/{{sample}}_{{mm_rate}}_", path=config["seq_output_dir"]),
#    mismatch_rate="{mm_rate}"
#  threads: 6
#  conda: config["envs_folder"] + "TE_RNASeq.yml"
#  shell:
#    "STAR --runThreadN {threads} --genomeDir {params.genome_folder} --readFilesCommand zcat --readFilesIn {input.R1} {input.R2} --outFileNamePrefix {params.output_prefix} --runMode alignReads --outSAMtype BAM Unsorted --outFilterMultimapNmax 5000 --outSAMmultNmax 1 --outFilterMismatchNmax 3 --outMultimapperOrder Random --winAnchorMultimapNmax 5000 --alignEndsType EndToEnd --alignIntronMax 1000000 --alignMatesGapMax 0 --seedSearchStartLmax 30 --alignTranscriptsPerReadNmax 30000 --alignWindowsPerReadNmax 30000 --alignTranscriptsPerWindowNmax 300 --seedPerReadNmax 3000 --seedPerWindowNmax 300 --seedNoneLociPerWindow 1000 --outReadsUnmapped Fastx --outSAMunmapped Within --outFilterScoreMinOverLread {params.mismatch_rate} --outFilterMatchNminOverLread {params.mismatch_rate}"

rule sort_annotation:
  input:
    TE_annotation=config["TE_annotation_bed"],
    gene_annotation=config["gene_annotation_bed"]
  output:
    sorted_TE_anno=config["TE_annotation_bed"].replace(".bed", ".sorted.bed"),
    sorted_gene_anno=config["gene_annotation_bed"].replace(".bed", ".sorted.bed")
  threads: 4
  shell:
    """
    sort -k1,1 -k2,2n {input.TE_annotation} > {output.sorted_TE_anno}
    sort -k1,1 -k2,2n {input.gene_annotation} > {output.sorted_gene_anno}
    """

rule filter_paired_bam:
  input:
    alignment=expand("{path}alignment/STAR/{{sample}}_amgm350_Aligned.out.bam", path=config["seq_output_dir"]),
    TE_annotation=rules.sort_annotation.output.sorted_TE_anno,
    gene_annotation=rules.sort_annotation.output.sorted_gene_anno
  output:
    bam_majchr=os.path.join(config["seq_output_dir"], "alignment/STAR/{sample}_amgm350_Aligned.majchr.sorted.bam"),
    TE_only_alignment=os.path.join(config["seq_output_dir"], "alignment/STAR/{sample}_amgm350.TE_only.bam")
  params:
    majchr_bed=config["maj_chr_bed"]
  threads: 4
  conda: config["envs_folder"] + "TE_RNASeq.yml"
  shell:
    """
    #first keeping only alignments on major chromosomes and sort by coordinates. Sorting is for speeding up intersectBed inside filter_paired_bam script
    samtools view -@ {threads} -L {params.majchr_bed} {input.alignment} -b | samtools sort -@ {threads} -T tmp -o {output.bam_majchr}
    #script which keep only reads from pairs where both mates are completely included into TE and not overlapping with gene bodies
    {config[script_dir]}sh/filter_paired_bam.sh {threads} {input.TE_annotation} {input.gene_annotation} {output.bam_majchr} {config[tmp_dir]}
    """

rule bed_to_saf:
  input:
    TE_annotation=rules.sort_annotation.output.sorted_TE_anno,
    saf_header_TE=os.path.join(config["annotation_dir"], "saf_header_TE")
  output:
    TE_SAF=config["TE_annotation_bed"].replace(".bed", ".SAF.noSimple.tsv"),
    TE_anno_fl_only=config["TE_annotation_bed"].replace(".bed", ".fl_only.bed"),
    TE_SAF_fl_only=config["TE_annotation_bed"].replace(".bed", ".fl_only.SAF.tsv")
  threads: 2
  shell:
    """
    #conversion to SAF format is necessary for featureCounts command in 'count_TE' rule; filtering out simple repeats and low complexity sequences is necessary in order to avoid that - in cases where these kind of repeats are inside other types of repeats - alignments are ambigously counted to both types of repeats by featureCounts
    cat {input.TE_annotation} | awk -v OFS="\t" '{{gsub(/\//, "\t", $7); print $4,$1,$2,$3,$6,$7}}' | cat {input.saf_header_TE} - | grep -E -w -v "Low_complexity|Simple_repeat" - > {output.TE_SAF}
    #creating a SAF file also for annotation containing only full-lengh RNA TE
    cat {input.TE_annotation} | awk '($7~/\/L1/ && $3-$2>5000) || ($4~/IAP/ && $3-$2>6000) || ($4~/MMERVK10C/ && $3-$2>4500){{print}}' > {output.TE_anno_fl_only}
    cat {output.TE_anno_fl_only} | awk -v OFS="\t" '{{gsub(/\//, "\t", $7); print $4,$1,$2,$3,$6,$7}}' | cat {input.saf_header_TE} - > {output.TE_SAF_fl_only}
    """

rule filter_TE_annotation:
  input:
    TE_annotation=config["TE_annotation_bed"],
    gene_annotation=config["gene_annotation_bed"]
  output:
    config["TE_annotation_bed"].replace(".bed", ".noGenes.noSimple.bed")
  params:
    maj_chr_bed=config["maj_chr_bed"]
    threads: 4
    conda: config["envs_folder"] + "TE_RNASeq.yml"
    shell:
    """
    bedtools intersect -a {input.TE_annotation} -b {params.maj_chr_bed} -u | bedtools intersect -a stdin -b {input.gene_annotation} -v | grep -E -w -v "Low_complexity|Simple_repeat" > {output}
    """

rule count_TE:
  input:
    alignment=rules.filter_paired_bam.output.TE_only_alignment,
    TE_SAF=rules.bed_to_saf.output.TE_SAF,
    TE_SAF_fl_only=rules.bed_to_saf.output.TE_SAF_fl_only
  output:
    TE_counts=os.path.join(config["analysis_output_dir"], "featureCounts/{sample}_TE_counts.txt"),
    TE_counts_fl_only=os.path.join(config["analysis_output_dir"], "featureCounts/{sample}_TE_counts_fl_only.txt")
  threads: 4
  conda: config["envs_folder"] + "TE_RNASeq.yml"
  shell:
    """
    featureCounts -F SAF -p -B -s 0 --fracOverlap 1 -M -T {threads} -a {input.TE_SAF} -o {output.TE_counts} {input.alignment}
    featureCounts -F SAF -p -B -s 0 --fracOverlap 1 -M -T {threads} -a {input.TE_SAF_fl_only} -o {output.TE_counts_fl_only} {input.alignment}

    """

onsuccess:
  if os.path.isdir('./tmp'):
    shell("rm -r tmp/")
    print("Workflow finished, no error")
  else:
    print("Workflow finished, no error")
  filename_re = 'fastqc'
  for filename in os.listdir(config["seq_output_dir"] + "trimming/"):
    if filename_re in filename:
      print("Moving trimmed fastqc files to qc directory since the latter is excluded from .gitignore (differently from the rest of 'data' subdirectories)")
      shell("mv {config[seq_output_dir]}trimming/*fastqc* {config[seq_output_dir]}qc/ && echo done")
  
