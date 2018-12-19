# sMACHETE

This repository contains the JSON script written in Common Workflow Language (CWL) for "scalable Mismatched Alignment CHimEra Tracking Engines"  (sMACHETE), which is a precise fusion detection algorithm specially designed for screening massive RNA sequencing databases. 

# Software Requirements

- Anaconda
- SciPy 1.1.0
- Bowtie2 2.2.9
- Python 3.4 with Biopython/1.70 installed
- java
- R 3.5.1
- Trim Galore! 0.4.4
# sMACHETE main script

All computational steps including any alignment step required for running sMACHETE has been packaged in a single JSON file "sMACHETE_pipeline.json". This script can be run on any local cluster using an input JSON file that provides the paths for reference index files and input RNA-seq file. To run the script, Rabix toolkit should be pre-installed on the local cluster.  
# Genome files

All Bowtie2 reference index files (genome, transcriptome, regular junctions, scrambled junctions, ribosome) have been pre processed and are available for human GRCh38 assembly. 

# Input file

All input parameters required for running sMACHETE JSON script should be provided via an input JSON file "sMACHETE_input.json". The following parameters should be set in the input JSON file:

- Bowtie2 index files for the reference genome
- Bowtie2 index files for the regular junctions
- Bowtie2 index files for the scrambled junctions
- Bowtie2 index files for ribosome
- Bowtie2 index files for the transcriptome
- Bowtie2 index files for indel junctions (for junctions with up to 5 symmetric indels at the splice site)
- Fastq files for the input RNA-Seq data
- Pickle files for genes/exons annotation
- Pickle file for the known fusions list (a list of known fusions constructed based on ChimerDB 3.0 curated list of known cancer fusions)

# Executing sMACHETE JSON script

For running sMACHETE JSON script on a local cluster, Rabix should be installed first. Rabix is an open-source tool developed by Seven Bridges, that can be used to run a computational workflow written in Common Workflow Language (CWL) on a locul cluster. More information on how to install Rabix can be found in this GitHub repositiory: https://github.com/rabix/bunny  

