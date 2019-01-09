# DEEPEST-Fusion

This repository contains the JSON script written in Common Workflow Language (CWL) for "Data-Enriched Efficient PrEcise STatistical Fusion detection" (DEEPEST-Fusion), which is a statistical fusion detection algorithm particularly engineered for screening big RNA sequencing databases. 

# Software Requirements

- Anaconda
- SciPy 1.1.0
- Bowtie2 2.2.9
- Python 3.4 with Biopython/1.70 installed
- java
- R 3.5.1
- Trim Galore! 0.4.4
# DEEPEST-Fusion main script

All computational steps including any alignment step required for running DEEPEST-Fusion has been packaged in a single JSON file "DEEPEST-Fusion_pipeline.json". This script can be run on any local cluster using an input JSON file that provides the paths for reference index files and input RNA-seq files. To run the script, Rabix toolkit should be pre-installed on the local cluster.  
# Genome files

All Bowtie2 reference index files (genome, transcriptome, regular junctions, scrambled junctions, ribosome) have been pre processed and are available for human GRCh38 assembly. 

# Input file

All input parameters required for running DEEPEST-Fusion JSON script should be provided via an input JSON file "DEEPEST-Fusion_input.json". The following parameters should be set in the input JSON file:

- Bowtie2 index files for the reference genome
- Bowtie2 index files for the regular junctions
- Bowtie2 index files for the scrambled junctions
- Bowtie2 index files for ribosome
- Bowtie2 index files for the transcriptome
- Bowtie2 index files for indel junctions (for junctions with up to 5 symmetric indels at the splice site)
- Fastq files for the input RNA-Seq data
- Pickle files for genes/exons annotation
- Pickle file for the known fusions list (a list of known fusions constructed based on ChimerDB 3.0 curated list of known cancer fusions)
- Bowtie2 index files for known fusions

Note: Instructions for building all index and reference files are provided here: 
https://github.com/salzmanlab/DEEPEST-Fusion/tree/master/reference_files 

# Toolkit for executing DEEPEST-Fusion JSON script

For running DEEPEST-Fusion JSON script on a local cluster, Rabix should be installed first. Rabix is an open-source tool developed by Seven Bridges, that can be used to run a computational workflow written in Common Workflow Language (CWL) on a locul cluster. More information on how to install Rabix can be found in this GitHub repositiory: https://github.com/rabix/bunny  

# An example Batch script for submitting DEEPEST-Fusion jobs on a local cluster

An example batch script "DEEPEST-Fusion_submit_job.sbatch", based on job scheduler Slurm has been provided. In the batch script file the path to where Rabix has been installed, DEEPEST-Fusion pipeline JSON file (DEEPEST-Fusion_pipeline.json), and DEEPEST-Fusion input JSON file (DEEPEST-Fusion_input.json) should be provided. 

# Scripts for implementing Sequence Bloom Trees for an RNA-Seq database

All scripts needed for implementing Sequence Bloom Tree (SBT) filters for an RNA-Seq database can be found in this github repository: https://github.com/elehnertSBG/SBT-Apps. Mor information on the SBT algorithm and the order in which the SBT CWL scripts should be run can be found in the original SBT manual: https://www.cs.cmu.edu/~ckingsf/software/bloomtree/sbt-manual.pdf 
